get.mod.fun <- function(mod.name) {
  x <- switch(mod.name, maxent.jar = dismo::maxent,
              maxnet = maxnet::maxnet,
              brt = dismo::gbm.step,
              bioclim = dismo::bioclim)
  return(x)
}

model.msgs <- function(tune.args, mod.name) {
  # tune.args checks and algorithm-specific set-up
  if(mod.name %in% c("maxent", "maxnet")) {
    if(!("rm" %in% names(tune.args)) | !("fc" %in% names(tune.args))) {
      stop("For Maxent, please specify both 'rm' and 'fc' settings. See ?tune.args for help.")
    }else{
      if(!is.numeric(tune.args[["rm"]])) {
        stop("Please input numeric values for 'rm' settings for Maxent.")
      }
      all.fc <- unlist(sapply(1:5, function(x) apply(combn(c("L","Q","H","P","T"), x), 2, function(y) paste(y, collapse = ""))))
      if(any(!tune.args[["fc"]] %in% all.fc)) {
        stop("Please input accepted values for 'fc' settings for Maxent.")
      }
    }
    
    if(mod.name == 'maxent') {
      if(is.null(getOption('dismo_rJavaLoaded'))) {
        # to avoid trouble on macs
        Sys.setenv(NOAWT=TRUE)
        if ( requireNamespace('rJava') ) {
          rJava::.jpackage('dismo')
          options(dismo_rJavaLoaded=TRUE)
        } else {
          stop('rJava cannot be loaded')
        }
      }
      mxe <- rJava::.jnew("meversion") 
      v <- try(rJava::.jcall(mxe, "S", "meversion"))
      algorithm.ver <- paste("maxent.jar v.", v, "from dismo package v.", packageVersion('dismo'))
    }
    
    if(mod.name == 'maxnet') {
      # construct user message with version info
      algorithm.ver <- paste("maxnet from maxnet package v.", packageVersion('maxnet'))
    }
  }
  
  if(mod.name == 'gbm.step') {
    if(!all("tree.complexity" %in% names(tune.args), "learning.rate" %in% names(tune.args), "bag.fraction" %in% names(tune.args))) {
      stop("BRT settings must include 'tree.complexity', 'learning.rate', and 'bag.fraction'.")
    }
    # construct user message with version info
    algorithm.ver <- paste("Boosted regression trees (BRTs) using the gbm.step() function from gbm package v.", packageVersion('gbm'), "and dismo package v.", packageVersion('dismo'))
  }
  
  if(mod.name == 'bioclim') {
    algorithm.ver <- paste("BIOCLIM from dismo v.", packageVersion('dismo'))
  }
  
  if(mod.name == "gam") {
    algorithm.ver <- paste("Generalized additive model (GAM) from gam package v.", packageVersion('gam'))
  }
  message(paste("*** Running ENMevaluate using", algorithm.ver, "***"))
}

#' @export
model.args <- function(tune.tbl.i, mod.name, occs.vals, bg.vals, other.args) {
  
  out <- list()
  # define data 
  d <- rbind(occs.vals, bg.vals)
  # define response
  p <- c(rep(1, nrow(occs.vals)), rep(0, nrow(bg.vals)))
  # maxent.jar
  if(mod.name == "maxent") {
    out$x <- d
    out$p <- p
    out$args <- c("noaddsamplestobackground", "noremoveDuplicates", "noautofeature")
    if(!grepl("L", tune.tbl.i$fc)) out$args <- c(out$args, "nolinear")
    if(!grepl("Q", tune.tbl.i$fc)) out$args <- c(out$args, "noquadratic")
    if(!grepl("H", tune.tbl.i$fc)) out$args <- c(out$args, "nohinge")
    if(!grepl("P", tune.tbl.i$fc)) out$args <- c(out$args, "noproduct")
    if(!grepl("T", tune.tbl.i$fc)) out$args <- c(out$args, "nothreshold")
    out$args <- c(out$args, paste0("betamultiplier=", tune.tbl.i$rm, sep=""))
  }
  # maxnet
  if(mod.name == "maxnet") {
    out$data <- d
    out$p <- p
    out$f <- maxnet::maxnet.formula(out$p, out$data, classes = tolower(tune.tbl.i$fc))
    out$regmult <- tune.tbl.i$rm
  }
  # BRT
  if(mod.name == "gbm.step") {
    out$tree.complexity <- tune.tbl.i$tree.complexity
    out$learning.rate <- tune.tbl.i$learning.rate
    out$bag.fraction <- tune.tbl.i$bag.fraction
    out$data <- cbind(p, d)
    out$gbm.x <- 2:ncol(out$data)
    out$gbm.y <- 1
    out$silent <- TRUE
  }
  # BIOCLIM
  if(mod.name == "bioclim") {
    out$x <- occs.vals 
  }
  
  # add other args
  out <- c(out, other.args)
  
  return(out)
}

# function to calculate AUC based on the model type
calcAUC <- function(occs.vals, bg.vals, mod, mod.name) {
  if(mod.name == "gbm.step") {
    auc <- dismo::evaluate(occs.vals, bg.vals, mod, n.trees = length(mod$trees))@auc
  }else{
    # currently, maxent.jar, maxnet, and BIOCLIM have no additional parameters
    auc <- dismo::evaluate(occs.vals, bg.vals, mod)@auc
  }
  
  return(auc)
}

# function to predict values to a raster based on the model type
rasterPred <- function(mod, envs, mod.name, doClamp) {
  if(mod.name == "maxent") {
    pred.args <- c("outputformat=raw", ifelse(doClamp == TRUE, "doclamp=true", "doclamp=false"))
    pred <- dismo::predict(mod, envs, args = pred.args, na.rm = TRUE)
  }else if(mod.name == "maxnet") {
    envs.n <- raster::nlayers(envs)
    envs.pts <- raster::rasterToPoints(envs)
    origNrow <- nrow(envs.pts)
    envs.pts <- na.omit(envs.pts)
    naOmitNrow <- nrow(envs.pts)
    rowDiff <- origNrow - naOmitNrow
    if(rowDiff > 0) {
      message(paste('\n', rowDiff, "grid cells found with at least one NA value: these cells were excluded from raster predictions."))
    }
    # mxnet.p <- maxnet::predict(mod, envs.pts, type=type, clamp=clamp)
    mxnet.p <- predict(mod, envs.pts, type = 'exponential', clamp = doClamp, na.rm = TRUE)
    envs.pts <- cbind(envs.pts, as.numeric(mxnet.p))
    pred <- raster::rasterFromXYZ(envs.pts[,c(1, 2, envs.n+3)], res=raster::res(envs))
  }else if(mod.name == "gbm.step") {
    pred <- raster::predict(envs, mod, type = "response", n.trees = mod$gbm.call$best.trees, na.rm = TRUE)
  }else if(mod.name == "bioclim") {
    # if no tails in other.args, defaults to NULL
    pred <- dismo::predict(mod, envs, tails = other.args$tails, na.rm = TRUE)
  }else{
    pred <- dismo::predict(mod, envs, na.rm = TRUE)
  }
  
  return(pred)
}

vectorPred <- function(mod, df, mod.name, doClamp) {
  if(mod.name == "maxent") {
    pred.args <- c("outputformat=raw", ifelse(doClamp == TRUE, "doclamp=true", "doclamp=false"))
    pred <- dismo::predict(mod, df, args = pred.args, na.rm = TRUE)
  }else if(mod.name == "maxnet") {
    pred <- dismo::predict(mod, df, type = 'exponential', clamp = doClamp, na.rm = TRUE)
  }else if(mod.name == "gbm.step") {
    pred <- dismo::predict(mod, df, type = "response", n.trees = mod$gbm.call$best.trees, na.rm = TRUE)
  }else if(mod.name == "bioclim") {
    # if no tails in other.args, defaults to NULL
    pred <- dismo::predict(mod, df, tails = other.args$tails, na.rm = TRUE)
  }else{
    pred <- dismo::predict(mod, df, na.rm = TRUE)
  }
  return(pred)
}

#' @export
# get total number of parameters
no.params <- function(mod, mod.name) {
  if(mod.name == 'maxnet') {
    return(length(mod$betas))
  }else if(mod.name == "maxent") {
    lambdas <- mod@lambdas[1:(length(mod@lambdas)-4)]
    countNonZeroParams <- function(x) if(strsplit(x, split=", ")[[1]][2] != '0.0') 1
    return(sum(unlist(sapply(lambdas, countNonZeroParams))))
  }else if(mod.name == "gbm.step") {
    # as no L1 regularization occurs, no parameters are dropped
    return(length(mod$var.names))
  }else if(mod.name == "bioclim") {
    return(length(mod@min))
  }
}

#' @title Calculate AICc from Maxent model prediction
#' @export
calc.aicc <- function(nparams, occs, preds) {
  # only functional for Maxent models currently
  out <- as.data.frame(matrix(nrow = length(nparams), ncol = 3, 
                              dimnames = list(NULL, c("AICc", "AICc.delta", "AICc.weights"))))
  AIC.valid <- nparams < nrow(occs)
  if(nlayers(preds) == 0) {
    warning("Cannot calculate AICc when rasterPreds = FALSE... returning NAs.")
  }else{
    vals <- raster::extract(preds, occs)
    probsum <- raster::cellStats(preds, sum)
    # The log-likelihood was incorrectly calculated (see next line) in ENMeval v.0.1.0 when working with >1 model at once.
    #   LL <- colSums(log(vals/probsum), na.rm=T)
    # The corrected calculation (since v.0.1.1) is:
    LL <- colSums(log(t(t(vals)/probsum)), na.rm=T)
    AICc <- (2*nparams - 2*LL) + (2*(nparams)*(nparams+1)/(nrow(occs)-nparams-1))
    AICc[AIC.valid==FALSE] <- NA
    AICc[is.infinite(AICc)] <- NA
    if(sum(is.na(AICc))==length(AICc)){
      warning("AICc not valid... returning NAs.")
    }else{
      out$AICc <- AICc
      out$AICc.delta <- (AICc - min(AICc, na.rm=TRUE))
      out$AICc.weights <- (exp(-0.5*out$AICc.delta))/(sum(exp(-0.5*out$AICc.delta), na.rm=TRUE))
    }    
  }
  return(out)
}

# define a corrected variance function
corrected.var <- function(x, nk){
  sum((x - mean(x))^2) * ((nk-1)/nk)
}

# repurposed from dismo::mess(), based on .messi3()
mess.vec <- function(p, v) {
  calc.mess <- function(p, v) {
    v <- stats::na.omit(v)
    f <- 100*findInterval(p, sort(v)) / length(v)
    minv <- min(v)
    maxv <- max(v)
    res <- 2*f 
    f[is.na(f)] <- -99
    i <- f>50 & f<100
    res[i] <- 200-res[i]
    
    i <- f==0 
    res[i] <- 100*(p[i]-minv)/(maxv-minv)
    i <- f==100
    res[i] <- 100*(maxv-p[i])/(maxv-minv)
    return(res)
  }
  
  x <- sapply(1:ncol(p), function(i) calc.mess(p[,i], v[,i]))
  rmess <- apply(x, 1, min, na.rm=TRUE)
  return(rmess)
}

#' @export

var.importance <- function(mod) {
  if(!'MaxEnt' %in% class(mod)){
    stop('Sorry, variable importance cannot currently be calculated with maxnet models (only maxent.jar)')
  } else {
    res <- mod@results
    pc <- res[grepl('contribution', rownames(res)),]
    pi <- res[grepl('permutation', rownames(res)),]
    varnames <- sapply(strsplit(names(pc), '.contribution'), function(x) x[1])
    df <- data.frame(variable=varnames, percent.contribution=pc, permutation.importance=pi, row.names=NULL)
    return(df)
  }
}
