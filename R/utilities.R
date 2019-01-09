#' @export
calc.aicc <- function(nparam, occ, predictive.maps) {
  AIC.valid <- nparam < nrow(occ)
  if (nlayers(predictive.maps) == 0) {
    res <- data.frame(cbind(AICc=NA, delta.AICc=NA, w.AIC=NA, parameters=nparam))
    warning("Cannot calculate AICc when rasterPreds = FALSE... returning NA's.")
  } else {
    vals <- extract(predictive.maps, occ)
    probsum <- cellStats(predictive.maps, sum)
    # The log-likelihood was incorrectly calculated (see next line) in ENMeval v.1.0.0 when working with >1 model at once.
    #   LL <- colSums(log(vals/probsum), na.rm=T)
    # The corrected calculation (since v.0.1.1) is:
    LL <- colSums(log(t(t(vals)/probsum)), na.rm=T)
    AICc <- (2*nparam - 2*LL) + (2*(nparam)*(nparam+1)/(nrow(occ)-nparam-1))
    AICc[AIC.valid==FALSE] <- NA
    AICc[is.infinite(AICc)] <- NA
    if(sum(is.na(AICc))==length(AICc)){
      warning("AICc not valid... returning NA's.")
      res <- data.frame(cbind(AICc, delta.AICc=NA, w.AIC=NA, parameters=nparam))
    } else {
      delta.AICc <- (AICc - min(AICc, na.rm=TRUE))
      w.AIC <- (exp(-0.5*delta.AICc))/(sum(exp(-0.5*delta.AICc), na.rm=TRUE))
      res <- data.frame(AICc, delta.AICc, w.AIC, parameters=nparam)
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
