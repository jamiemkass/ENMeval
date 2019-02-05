#' @export
make.args <- function(tune.args.i, mod.fun.name, occs.vals, bg.vals, other.args) {
  
  out <- list()
  
  if(mod.fun.name == "maxent") {
    out$args <- c("noaddsamplestobackground", "noremoveDuplicates", "noautofeature")
    
    if(!grepl("L", tune.args.i$fc)) out$args <- c(out$args, "nolinear")
    if(!grepl("Q", tune.args.i$fc)) out$args <- c(out$args, "noquadratic")
    if(!grepl("H", tune.args.i$fc)) out$args <- c(out$args, "nohinge")
    if(!grepl("P", tune.args.i$fc)) out$args <- c(out$args, "noproduct")
    if(!grepl("T", tune.args.i$fc)) out$args <- c(out$args, "nothreshold")
    out$args <- c(out$args, paste0("betamultiplier=", tune.args.i$rm, sep=""))
    out$x <- as.data.frame(rbind(occs.vals, bg.vals))
    out$p <- c(rep(1, nrow(occs.vals)), rep(0, nrow(bg.vals)))
  }
  
  if(mod.fun.name == "maxnet") {
    out$data <- rbind(occs.vals, bg.vals)
    out$p <- c(rep(1, nrow(occs.vals)), rep(0, nrow(bg.vals)))
    out$f <- maxnet::maxnet.formula(out$p, out$data, classes = tolower(tune.args.i$fc))
    out$regmult <- tune.args.i$rm
  }
  
  # need logic for BRTs
  
  # add other args
  out <- c(out, other.args)
  
  return(out)
}

#' @export
# get total number of parameters
getNoParams <- function(m, mod.fun.name) {
  if(mod.fun.name == 'maxnet') {
    return(length(m$betas))
  }
  if(mod.fun.name == "maxent") {
    lambdas <- m@lambdas[1:(length(m@lambdas)-4)]
    countNonZeroParams <- function(x) {
      if(strsplit(x, split=", ")[[1]][2] != '0.0') 1
    }
    return(sum(unlist(sapply(lambdas, countNonZeroParams))))
  }
}

#' @export
calc.aicc <- function(nparam, occs, preds) {
  AIC.valid <- nparam < nrow(occs)
  if (nlayers(preds) == 0) {
    res <- data.frame(cbind(AICc = NA, AICc.delta = NA, AICc.weights = NA, params = nparam))
    warning("Cannot calculate AICc when rasterPreds = FALSE... returning NA's.")
  } else {
    vals <- raster::extract(preds, occs)
    probsum <- raster::cellStats(preds, sum)
    # The log-likelihood was incorrectly calculated (see next line) in ENMeval v.1.0.0 when working with >1 model at once.
    #   LL <- colSums(log(vals/probsum), na.rm=T)
    # The corrected calculation (since v.0.1.1) is:
    LL <- colSums(log(t(t(vals)/probsum)), na.rm=T)
    AICc <- (2*nparam - 2*LL) + (2*(nparam)*(nparam+1)/(nrow(occs)-nparam-1))
    AICc[AIC.valid==FALSE] <- NA
    AICc[is.infinite(AICc)] <- NA
    if(sum(is.na(AICc))==length(AICc)){
      warning("AICc not valid... returning NA's.")
      res <- data.frame(cbind(AICc, AICc.delta=NA, AICc.weights=NA, params = nparam))
    } else {
      AICc.delta <- (AICc - min(AICc, na.rm=TRUE))
      AICc.weights <- (exp(-0.5*AICc.delta))/(sum(exp(-0.5*AICc.delta), na.rm=TRUE))
      res <- data.frame(AICc, AICc.delta, AICc.weights, params = nparam)
      rownames(res) <- NULL
    }    
  }
  rownames(res) <- NULL
  return(res)
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

#' @export

# function to make a raster prediction from a maxnet object
maxnet.predictRaster <- function(mod, env, type, clamp) {
  env.n <- raster::nlayers(env)
  env.pts <- raster::rasterToPoints(env)
  origNrow <- nrow(env.pts)
  env.pts <- na.omit(env.pts)
  naOmitNrow <- nrow(env.pts)
  rowDiff <- origNrow - naOmitNrow
  if (rowDiff > 0) {
    message(paste('\n', rowDiff, "grid cells found with at least one NA value: these cells were excluded from raster predictions."))
  }
  # mxnet.p <- maxnet::predict(mod, env.pts, type=type, clamp=clamp)
  mxnet.p <- predict(mod, env.pts, type=type, clamp=clamp)
  env.pts <- cbind(env.pts, as.numeric(mxnet.p))
  mxnet.p <- raster::rasterFromXYZ(env.pts[,c(1, 2, env.n+3)], res=raster::res(env))
  return(mxnet.p)
}

#' @export
maxentJARversion <- function() {
  if (is.null(getOption('dismo_rJavaLoaded'))) {
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
  return(v)
}
