get.mod.fun <- function(mod.name) {
  x <- switch(mod.name, maxent.jar = dismo::maxent,
              maxnet = maxnet::maxnet,
              brt = dismo::gbm.step,
              bioclim = dismo::bioclim)
  return(x)
}

mod.msgs <- function(tune.args, mod.name) {
  # tune.args checks and algorithm-specific set-up
  if(mod.name %in% c("maxent.jar", "maxnet")) {
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
    
    if(mod.name == 'maxent.jar') {
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
  
  if(mod.name == 'brt') {
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
  message(paste("*** Running ENMevaluate using", algorithm.ver, "***\n"))
}

#' @export
mod.args <- function(tune.tbl.i, mod.name, occs.vals, bg.vals, other.args) {
  
  out <- list()
  # define data 
  d <- rbind(occs.vals, bg.vals)
  # define response
  p <- c(rep(1, nrow(occs.vals)), rep(0, nrow(bg.vals)))
  # maxent.jar
  if(mod.name == "maxent.jar") {
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
  if(mod.name == "brt") {
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
mod.eval <- function(occs.vals, bg.vals, mod, mod.name, doClamp) {
  if(mod.name == "maxent.jar") {
    pred.args <- c("outputformat=raw", ifelse(doClamp == TRUE, "doclamp=true", "doclamp=false"))
    mod.eval <- dismo::evaluate(occs.vals, bg.vals, mod, args = pred.args)
  }else if(mod.name == "maxnet") {
    mod.eval <- dismo::evaluate(occs.vals, bg.vals, mod, type = 'exponential', clamp = doClamp)
  }else if(mod.name == "brt") {
    mod.eval <- dismo::evaluate(occs.vals, bg.vals, mod, n.trees = length(mod$trees))
  }else if(mod.name == "bioclim") {
    # if no tails in other.args, defaults to NULL
    mod.eval <- dismo::evaluate(occs.vals, bg.vals, mod, tails = other.args$tails)
  }else{
    # currently, maxent.jar, maxnet, and BIOCLIM have no additional parameters
    mod.eval <- dismo::evaluate(occs.vals, bg.vals, mod)
  }
  return(mod.eval)
}

# function to predict values to a raster based on the model type
mod.prediction <- function(mod, envs, mod.name, other.args, doClamp) {
  if(mod.name == "maxent.jar") {
    pred.args <- c("outputformat=raw", ifelse(doClamp == TRUE, "doclamp=true", "doclamp=false"))
    pred <- dismo::predict(mod, envs, args = pred.args, na.rm = TRUE)
  }else if(mod.name == "maxnet") {
    # if input envs is Raster*, return a RasterLayer of predicted values
    if(inherits(envs, "BasicRaster") == TRUE) {
      envs.n <- raster::nlayers(envs)
      envs.pts <- raster::rasterToPoints(envs)
      # mxnet.p <- maxnet::predict(mod, envs.pts, type=type, clamp=clamp)
      mxnet.p <- predict(mod, envs.pts, type = 'exponential', clamp = doClamp, na.rm = TRUE)
      envs.pts <- cbind(envs.pts, as.numeric(mxnet.p))
      pred <- raster::rasterFromXYZ(envs.pts[,c(1, 2, envs.n+3)], res=raster::res(envs))  
    }else{
      # otherwise, envs is data frame, so return data frame of predicted values
      pred <- dismo::predict(mod, envs, type = 'exponential', clamp = doClamp, na.rm = TRUE)
    }
  }else if(mod.name == "brt") {
    pred <- raster::predict(envs, mod, type = "response", n.trees = mod$gbm.call$best.trees, na.rm = TRUE)
  }else if(mod.name == "bioclim") {
    # if no tails in other.args, defaults to NULL
    pred <- dismo::predict(mod, envs, tails = other.args$tails, na.rm = TRUE)
  }else{
    pred <- dismo::predict(mod, envs, na.rm = TRUE)
  }
  
  return(pred)
}

#' @export
# get total number of parameters
mod.paramNum <- function(mod, mod.name) {
  if(mod.name == 'maxnet') {
    return(length(mod$betas))
  }else if(mod.name == "maxent.jar") {
    lambdas <- mod@lambdas[1:(length(mod@lambdas)-4)]
    countNonZeroParams <- function(x) if(strsplit(x, split=", ")[[1]][2] != '0.0') 1
    return(sum(unlist(sapply(lambdas, countNonZeroParams))))
  }else if(mod.name == "brt") {
    # as no L1 regularization occurs, no parameters are dropped
    return(length(mod$var.names))
  }else if(mod.name == "bioclim") {
    return(length(mod@min))
  }
}

#' @title Calculate AICc from Maxent model prediction
#' @description This function calculates AICc for Maxent models based on Warren 
#' and Seifert (2011).
#' @param nparams integer; number of model parameters (non-zero coefficients)
#' @param occs data frame; longitude (x) and latitude (y) of occurrence 
#' localities
#' @param preds Raster*; Maxent model predictions from \code{dismo::predict()}
#' @return data frame with four columns:
#' \code{AICc} is the Akaike Information Criterion corrected for small sample
#' sizes calculated as:
#' \deqn{ (2 * K - 2 * logLikelihood) + (2 * K) * (K+1) / (n - K - 1)}
#' where \emph{K} is the number of parameters in the model (i.e., number of
#' non-zero parameters in Maxent lambda file) and \emph{n} is the number of
#' occurrence localities.  The \emph{logLikelihood} is calculated as:
#' \deqn{ sum(log(vals / total))}
#' where \emph{vals} is a vector of Maxent raw values at occurrence localities
#' and \emph{total} is the sum of Maxent raw values across the entire study
#' area.
#' \code{delta.AICc} is the difference between the AICc of a given model and
#' the AICc of the model with the lowest AICc.
#' \code{w.AICc} is the Akaike weight (calculated as the relative likelihood of
#' a model (exp(-0.5 * \code{delta.AICc})) divided by the sum of the likelihood
#' values of all models included in a run.  These can be used for model
#' averaging (Burnham and Anderson 2002).
#' @aliases calc.aicc get.params
#' @details As motivated by Warren and Seifert (2011) and implemented in 
#' ENMTools (Warren  \emph{et al.} 2010), this function calculates the small 
#' sample size version of Akaike Information Criterion for ENMs (Akaike 1974).  
#' We use AICc (instead of AIC) regardless of sample size based on the 
#' recommendation of Burnham and Anderson (1998, 2004).  The number of 
#' parameters is determined by counting the number of non-zero parameters in 
#' the \code{maxent} lambda file.  See Warren \emph{et al.} (2014) for 
#' limitations of this approach, namely that the number of parameters is an 
#' estimate of the true degrees of freedom.  For Maxent ENMs, AICc is 
#' calculated by standardizing the raw output such that all cells in the study 
#' extent sum to 1.  The likelihood of the data for a given model is then 
#' calculated by taking the product of the raw output values for all grid cells 
#' that contain an occurrence locality (Warren and Seifert 2011).
#' @seealso \code{maxent} in the \pkg{dismo} package.
#' @note Returns all \code{NA}s if the number of parameters is larger than the
#' number of observations (occurrence localities).
#' @references Akaike, H. (1974) A new look at the statistical model
#' identification. \emph{IEEE Transactions on Automatic Control}, \bold{19}:
#' 716-723.
#' 
#' Burnham, K. P. and Anderson, D. R. (1998) Model selection and multimodel
#' inference: a practical information-theoretic approach. Springer, New York.
#' 
#' Burnham, K. P. and Anderson, D. R. (2004) Multimodel inference:
#' understanding AIC and BIC in model selection. \emph{Sociological Methods and
#' Research}, \bold{33}: 261-304.
#' 
#' Warren, D. L., Glor, R. E, and Turelli, M. (2010) ENMTools: a toolbox for
#' comparative studies of environmental niche models. \emph{Ecography},
#' \bold{33}: 607-611.
#' 
#' Warren, D. L. and Seifert, S. N. (2011) Ecological niche modeling in Maxent:
#' the importance of model complexity and the performance of model selection
#' criteria. \emph{Ecological Applications}, \bold{21}: 335-342.
#' 
#' Warren, D. L., Wright, A. N., Seifert, S. N., and Shaffer, H. B. (2014)
#' Incorporating model complexity and sampling bias into ecological niche
#' models of climate change risks faced by 90 California vertebrate species of
#' concern. \emph{Diversity and Distributions}, \bold{20}: 334-343.

#' @export
calc.aicc <- function(nparams, occs, preds) {
  # only functional for Maxent models currently
  out <- as.data.frame(matrix(nrow = length(nparams), ncol = 3, 
                              dimnames = list(NULL, c("AICc", "AICc.delta", "AICc.weights"))))
  AIC.valid <- nparams < nrow(occs)
  if(nlayers(preds) == 0) {
    warning("Cannot calculate AICc when skipRasters = TRUE... returning NAs.")
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

# var.importance <- function(mod) {
#   if(!'MaxEnt' %in% class(mod)){
#     stop('Sorry, variable importance cannot currently be calculated with maxnet models (only maxent.jar)')
#   } else {
#     res <- mod@results
#     pc <- res[grepl('contribution', rownames(res)),]
#     pi <- res[grepl('permutation', rownames(res)),]
#     varnames <- sapply(strsplit(names(pc), '.contribution'), function(x) x[1])
#     df <- data.frame(variable=varnames, percent.contribution=pc, permutation.importance=pi, row.names=NULL)
#     return(df)
#   }
# }

calc.niche.overlap <- function(preds, overlapStat){
  n <- raster::nlayers(preds)
  ov <- matrix(nrow = n, ncol = n)
  pb <- txtProgressBar(0, nlayers(preds) - 1, style = 3)
  for(i in 1:(n - 1)){
    setTxtProgressBar(pb, i)
    for(j in (i + 1):n){
      ov[j, i] <- dismo::nicheOverlap(preds[[i]], preds[[j]], stat = overlapStat)
    }
  }
  colnames(ov) <- names(preds)
  rownames(ov) <- names(preds)
  close(pb)
  return(ov)
}
