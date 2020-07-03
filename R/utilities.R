#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#' @export
# remove.env.na <- function(d) {
#   d.envs <- d[,3:ncol(d)]
#   ind.NA <- unique(which(is.na(d.envs), arr.ind = TRUE)[,1])
#   names(ind.NA) <- NULL
#   d.envs.NA <- d.envs[ind.NA,]
#   ind.NA.occs <- ind.NA[which(d.envs.NA$pb == 1)]
#   ind.NA.bg <- ind.NA[which(d.envs.NA$pb == 0)]
#   occs.msg <- paste0("occurrences: ", paste(ind.NA.occs, collapse = ","))
#   bg.msg <- paste0("background: ", paste(ind.NA.bg, collapse = ","))
#   if(length(ind.NA.occs) > 0 | length(ind.NA.bg) > 0) {
#     envNA.msg <- dplyr::case_when(length(ind.NA.occs) > 0 & length(ind.NA.bg) > 0 ~ paste(occs.msg, bg.msg, sep = ", "),
#                             length(ind.NA.occs) > 0 ~ occs.msg,
#                             length(ind.NA.bg) > 0 ~ bg.msg)
#     message(paste0("* Records found with NA for at least one predictor variable with the following row numbers: (", envNA.msg, "). Removing from analysis..."))
#     d.naRem <- d[-c(ind.NA.occs, ind.NA.bg),]
#     return(d.naRem)    
#   }
#   return(d)
# }

#' @title Calculate AICc from Maxent model prediction
#' @description This function calculates AICc for Maxent models based on Warren 
#' and Seifert (2011).
#' @param p.occs data frame of raw (maxent.jar) or exponential (maxnet) predictions for
#' the occurrence localities based on one or more models
#' @param nparam integer; number of model parameters (non-zero coefficients)
#' @param p Raster* of raw (maxent.jar) or exponential (maxnet) model predictions;
#' if NULL, AICc will be calculated based on the background points, which already
#' have predictions that sum to 1 and thus need no correction -- this assumes that
#' the background points represent a good sample of the study extent 
#' @return data frame with three columns:
#' \code{AICc} is the Akaike Information Criterion corrected for small sample
#' sizes calculated as:
#' \deqn{ (2 * K - 2 * logLikelihood) + (2 * K) * (K+1) / (n - K - 1)}
#' where \emph{K} is the number of parameters in the model (i.e., number of
#' non-zero parameters) and \emph{n} is the number of
#' occurrence localities.  The \emph{logLikelihood} is calculated as:
#' \deqn{ sum(log(vals))}
#' where \emph{vals} is a vector of Maxent raw/expontential values at occurrence localities
#' and the sum of these values across the study extent is equal to 1.
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
#' the \code{maxent} lambda file (\code{m@lambdas} for maxent.jar and \code{m$betas} for maxnet.  
#' See Warren \emph{et al.} (2014) for limitations of this approach, namely that the number of parameters is an 
#' estimate of the true degrees of freedom.  For Maxent ENMs, AICc is 
#' calculated by standardizing the raw output such that all cells in the study 
#' extent sum to 1. The likelihood of the data for a given model is then 
#' calculated by taking the product of the raw output values (or the sum of their logs, as is implemented here)
#' for all grid cells that contain an occurrence locality (Warren and Seifert 2011).
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
aic.maxent <- function(p.occs, nparams, p = NULL) {
  # differential behavior for summing if p is Raster* or data frame
  if(!is.null(p)) {
    p.sum <- raster::cellStats(p, sum)  
    # if total does not sum to 1, standardize so that the sum is 1
    for(i in 1:nlayers(p)) if(p.sum[i] != 1) p.occs[,i] <- p.occs[,i] / p.sum[i]
  }
  # if more parameters than data points, determine AIC to be invalid:
  # this avoids considering overly complex models at all
  n.occs <- nrow(p.occs)
  AIC.valid <- nparams < n.occs
  # calculate log likelihood
  LL <- colSums(log(p.occs), na.rm = TRUE)
  AICc <- (2 * nparams - 2 * LL) + (2 * (nparams) * (nparams + 1) / (n.occs - nparams - 1))
  # if determined invalid or if infinite, make AICc NA
  AICc <- sapply(1:length(AICc), function(x) ifelse(AIC.valid[x] == FALSE | is.infinite(AICc[x]), NA, AICc[x]))
  # make output table
  out <- data.frame(AICc = AICc, delta.AICc = AICc - min(AICc, na.rm=TRUE))
  out$w.AIC <- exp(-0.5 * out$delta.AICc) / sum(exp(-0.5 * out$delta.AICc), na.rm=TRUE)
  return(out)
}

# Define a corrected variance function
#' Calculate variance corrected for non-independence of \emph{k}-fold iterations
#'
#' `corrected.var` calculates variance corrected for non-independence of \emph{k}-fold iterations.  See Appendix of Shcheglovitova & Anderson (2013) and other references (Miller 1974; Parr 1985; Shao and Wu 1989) for additional details. 
#' 
#' This function calculates variance that is corrected for the non-independence of \emph{k} cross-validation iterations.  Following Shao and Wu (1989): 
#' 
#' \deqn{Sum Of Squares * ((n-1)/n)} 
#' 
#' where \emph{n} = the number of \emph{k}-fold iterations.
#'
#' @param x A numeric vector.
#' @param nk Number of \emph{k}-fold iterations.
#' @return A numeric value of the corrected variance.
#' @author Robert Muscarella <bob.muscarella@gmail.com>
#' @references 
#'   Miller, R. G. (1974) The jackknife - a review. \emph{Biometrika}, \bold{61}: 1-15.
#'   
#'   Parr, W. C. (1985) Jackknifing differentiable statistical functionals. \emph{Journal of the Royal Statistics Society, Series B}, \bold{47}: 56-66.
#'   
#'   Shao J. and Wu, C. F. J. (1989) A general theory for jackknife variance estimation. \emph{Annals of Statistics}, \bold{17}: 1176-1197.
#'   
#'   Shcheglovitova, M. and Anderson, R. P. (2013) Estimating optimal complexity for ecological niche models: a jackknife approach for species with small sample sizes. \emph{Ecological Modelling}, \bold{269}: 9-17.
#'   
#' @export
corrected.var <- function(x, nk){
  sum((x - mean(x))^2) * ((nk-1)/nk)
}

# function to calculate the 10 percentile threshold from training predictions
#' @export
calc.10p.trainThresh <- function(pred.train) {
  n <- length(pred.train)
  if(n < 10) {
    pct90.train <- floor(n * 0.9)
  }else{
    pct90.train <- ceiling(n * 0.9)
  }
  pct10.train.thr <- rev(sort(pred.train))[pct90.train]
  return(pct10.train.thr)
}


#' @title Compute multivariate environmental similarity surfaces (MESS)
#' @description Compute multivariate environmental similarity surfaces (MESS) (i.e., Elith \emph{et al.} 2010).
#' @details Repurposed from dismo::mess(), based on .messi3()
#' @param p Raster object
#' @param v Matrix with reference values
#' @return 
#' A RasterBrick with layers corresponding to the input layers and an additional layer with the mess values (if full=TRUE and nlayers(x) > 1) or a RasterLayer with the MESS values (if full=FALSE).
#' @references 
#' Elith J., M. Kearney M., and Phillips, S. (2010) The art of modelling range-shifting species. \emph{Methods in Ecology and Evolution}, \bold{1}: 330-342.
#' @author 
#' Based on \pkg{dismo}::\code{mess}
#' Jean-Pierre Rossi <jean-pierre.rossi@supagro.inra.fr>, Robert Hijmans, Paulo van Breugel

# mess.vec <- function(p, v) {
#   calc.mess <- function(p, v) {
#     v <- stats::na.omit(v)
#     f <- 100*findInterval(p, sort(v)) / length(v)
#     minv <- min(v)
#     maxv <- max(v)
#     res <- 2*f 
#     f[is.na(f)] <- -99
#     i <- f>50 & f<100
#     res[i] <- 200-res[i]
#     
#     i <- f==0 
#     res[i] <- 100*(p[i]-minv)/(maxv-minv)
#     i <- f==100
#     res[i] <- 100*(maxv-p[i])/(maxv-minv)
#     return(res)
#   }
#   
#   x <- sapply(1:ncol(p), function(i) calc.mess(p[,i], v[,i]))
#   rmess <- apply(x, 1, min, na.rm=TRUE)
#   return(rmess)
# }
# 
# calc.mess.kstats <- function(occs.train.vals, bg.train.vals, occs.test.vals, bg.test.vals) {
#   p <- rbind(occs.train.vals, bg.train.vals)
#   v <- rbind(occs.test.vals, bg.test.vals)
#   cat.j <- which(sapply(occs.train.vals, is.factor) == 1)
#   if(length(cat.j) > 0) {
#     p <- p[,-cat.j]
#     v <- v[,-cat.j]
#   }
#   mss <- mess.vec(p, v)
#   mess.quant <- quantile(mss)
#   names(mess.quant) <- paste0("mess.", gsub("%", "p", names(mess.quant)))
#   return(mess.quant)
# }


#' @title Calculate Similarity of ENMs in Geographic Space
#' 
#' @description Compute pairwise "niche overlap" in geographic space for Maxent predictions. The value ranges from 0 (no overlap) to 1 (identical predictions).  The function uses the \code{nicheOverlap} function of the \pkg{dismo} package (Hijmans \emph{et al.} 2011).
#' 
#' @aliases calc.niche.overlap
#' @usage 
#' calc.niche.overlap(predictive.maps, overlapStat = "D", maxent.args)
#' @param envs A rasterStack of at least 2 Maxent predictive raster layers.
#' @param overlapStat The statistic calculated by the \code{nicheOverlap} function of the \pkg{dismo} package.  Defaults to Schoeners \emph{D} (Schoener 1968) but can also accept \code{"I"} to calculate the \emph{I} similarity statistic from Warren \emph{et al.} (2008).
#' @return 
#' A matrix with the lower triangle giving values of pairwise "niche overlap" in geographic space.  Row and column names are given by the \code{\link{make.args}} argument when run by the \code{\link{ENMevaluate}} function.
#' @references 
#' Hijmans, R. J., Phillips, S., Leathwick, J. and Elith, J. (2011) dismo package for R. Available online at: \url{https://cran.r-project.org/package=dismo}.
#' Schoener, T. W. (1968) The \emph{Anolis} lizards of Bimini: resource partitioning in a complex fauna. \emph{Ecology}, \bold{49}: 704-726.
#' Warren, D. L., Glor, R. E., Turelli, M. and Funk, D. (2008) Environmental niche equivalency versus conservatism: quantitative approaches to niche evolution. \emph{Evolution}, \bold{62}: 2868-2883.
#' @author 
#' Based on \pkg{dismo}::\code{nicheOverlap}, which is based on \pkg{SDMTools}::\code{Istat}
#' Robert Muscarella <bob.muscarella@gmail.com>
#' @seealso 
#' `nicheOverlap` in the \pkg{dismo} package

#' @export
calc.niche.overlap <- function(envs, overlapStat, quiet=FALSE){
  n <- raster::nlayers(envs)
  ov <- matrix(nrow = n, ncol = n)
  if(quiet != TRUE) pb <- txtProgressBar(0, n - 1, style = 3)
  for(i in 1:(n - 1)){
    if(quiet != TRUE) setTxtProgressBar(pb, i)
    for(j in (i + 1):n){
      ov[j, i] <- dismo::nicheOverlap(envs[[i]], envs[[j]], stat = overlapStat)
    }
  }
  colnames(ov) <- names(envs)
  rownames(ov) <- names(envs)
  if(quiet != TRUE) close(pb)
  return(ov)
}

# function to look up the corresponding ENMdetails abject
#' @export
lookup.enm <- function(mod.name) {
  x <- switch(mod.name, 
              maxent.jar = enm.maxent.jar,
              maxnet = enm.maxnet,
              brt = enm.brt,
              bioclim = enm.bioclim)
  return(x)
}

# Modified version of dismo::mess
# This version ignores raster cells with NA for every variable, thus avoiding
# the generation of -Inf values and the corresponding warnings
.messi3 <- function(p,v) {
  # seems 2-3 times faster than messi2
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
  res
}

#' @export
mess <- function(x, v, full=FALSE, filename='', ...) {
  
  stopifnot(NCOL(v) == nlayers(x))
  out <- raster(x)
  nl <- nlayers(x)
  filename <- trim(filename)
  nms <- paste(names(x), '_mess', sep='')
  
  if (canProcessInMemory(x)) {
    x <- getValues(x)
    if (nl == 1) {
      rmess <- .messi3(x, v)
      names(out) <- 'mess'
      out <- setValues(out, rmess)
    } else {
      x <- sapply(1:ncol(x), function(i) .messi3(x[,i], v[,i]))
      rmess <- apply(x, 1, function(x) ifelse(!all(is.na(x)), min(x, na.rm=TRUE), NA))
      if (full) {
        out <- brick(out, nl=nl+1)
        names(out) <- c(nms, "mess")
        out <- setValues(out, cbind(x, rmess))
      } else {
        names(out) <- 'mess'
        out <- setValues(out, rmess)
      }
    }	
    if (filename != '') {
      out <- writeRaster(out, filename, ...)
    }
    return(out)
    
  } else {
    
    if (nl == 1) {
      
      names(out) <- "mess"
      tr <- blockSize(out)
      pb <- pbCreate(tr$n, ...)	
      out <- writeStart(out, filename, ...)
      for (i in 1:tr$n) {
        vv <- getValues(x, row=tr$row[i], nrows=tr$nrows[i])
        vv <- .messi3(vv, v)
        out <- writeValues(out, vv, tr$row[i])
        pbStep(pb) 
      }
      
    } else {
      
      if (full) {
        out <- brick(out, nl=nl+1)
        names(out) <- c(nms, "mess")
        tr <- blockSize(out)
        pb <- pbCreate(tr$n, ...)	
        out <- writeStart(out, filename, ...)
        for (i in 1:tr$n) {
          vv <- getValues(x, row=tr$row[i], nrows=tr$nrows[i])
          vv <- sapply(1:ncol(v), function(i) .messi3(vv[,i], v[,i]))
          m <- apply(vv, 1, min, na.rm=TRUE)
          out <- writeValues(out, cbind(vv, m), tr$row[i])
          pbStep(pb) 
        }
        
      } else {
        
        names(out) <- "mess"
        tr <- blockSize(out)
        pb <- pbCreate(tr$n, ...)	
        out <- writeStart(out, filename, ...)
        for (i in 1:tr$n) {
          vv <- getValues(x, row=tr$row[i], nrows=tr$nrows[i])
          vv <- sapply(1:ncol(v), function(i) .messi3(vv[,i], v[,i]))
          m <- apply(vv, 1, min, na.rm=TRUE)
          out <- writeValues(out, m, tr$row[i])
          pbStep(pb) 
        }
      }
    }
    out <- writeStop(out)
    pbClose(pb) 
  }	
  out
}

# function to look up the version of maxent.jar
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
