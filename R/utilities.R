#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

msg <- function(x, quiet) if(quiet == FALSE) message(x)

#' @export
remove.env.na <- function(d, quiet = FALSE) {
  d.envs <- d[,3:ncol(d)]
  ind.NA <- unique(which(is.na(d.envs), arr.ind = TRUE)[,1])
  names(ind.NA) <- NULL
  d.envs.NA <- d.envs[ind.NA,]
  ind.NA.occs <- ind.NA[which(d.envs.NA$pb == 1)]
  ind.NA.bg <- ind.NA[which(d.envs.NA$pb == 0)]
  occs.msg <- paste0("occurrences: ", paste(ind.NA.occs, collapse = ","))
  bg.msg <- paste0("background: ", paste(ind.NA.bg, collapse = ","))
  if(length(ind.NA.occs) > 0 | length(ind.NA.bg) > 0) {
    envNA.msg <- dplyr::case_when(length(ind.NA.occs) > 0 & length(ind.NA.bg) > 0 ~ paste(occs.msg, bg.msg, sep = ", "),
                            length(ind.NA.occs) > 0 ~ occs.msg,
                            length(ind.NA.bg) > 0 ~ bg.msg)
    msg(paste0("* Records found with NA for at least one predictor variable with the following row numbers: (", envNA.msg, "). Removing from analysis...\n"), quiet)
    d.naRem <- d[-c(ind.NA.occs, ind.NA.bg),]
    return(d.naRem)    
  }
  return(d)
}

#' @export
maxnet.predictRaster <- function(mod, envs, doClamp, type, other.args) {
  if(inherits(envs, "BasicRaster") == TRUE) {
    envs.n <- raster::nlayers(envs)
    envs.pts <- na.omit(raster::rasterToPoints(envs))
    mxnet.p <- predict(mod, envs.pts, type = type, clamp = doClamp, na.rm = TRUE, other.args)
    p.vals <- cbind(envs.pts[,1:2], as.numeric(mxnet.p))
    p.ras <- raster::rasterFromXYZ(p.vals, res=raster::res(envs))
  }else{
    # otherwise, envs is data frame, so return data frame of predicted values
    p.ras <- dismo::predict(mod, envs, type = type, clamp = doClamp, na.rm = TRUE, other.args)
  }
  return(p.ras)
}

#' @title Calculate AICc from Maxent model prediction
#' @description This function calculates AICc for Maxent models based on Warren 
#' and Seifert (2011).
#' @param occs data frame; longitude (x) and latitude (y) of occurrence 
#' localities
#' @param nparam integer; number of model parameters (non-zero coefficients)
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
calc.aicc <- function(occs.preds, nparams, preds) {
  # only functional for Maxent models currently
  out <- as.data.frame(matrix(nrow = length(nparams), ncol = 3, 
                              dimnames = list(NULL, c("AICc", "delta.AICc", "w.AIC"))))
  AIC.valid <- nparams < nrow(occs.preds)
  if(raster::nlayers(preds) == 0) {
    warning("Cannot calculate AICc without prediction rasters... returning NAs.", immediate. = TRUE)
  }else{
    probsum <- raster::cellStats(preds, sum)
    # The log-likelihood was incorrectly calculated (see next line) in ENMeval v.0.1.0 when working with >1 model at once.
    #   LL <- colSums(log(vals/probsum), na.rm=T)
    # The corrected calculation (since v.0.1.1) is:
    LL <- colSums(log(t(t(occs.preds)/probsum)), na.rm=T)
    AICc <- (2*nparams - 2*LL) + (2*(nparams)*(nparams+1)/(nrow(occs.preds)-nparams-1))
    AICc[AIC.valid==FALSE] <- NA
    AICc[is.infinite(AICc)] <- NA
    if(sum(is.na(AICc))==length(AICc)){
      warning("AICc not valid... returning NAs.")
    }else{
      out$AICc <- AICc
      out$delta.AICc <- (AICc - min(AICc, na.rm=TRUE))
      out$w.AIC <- (exp(-0.5*out$delta.AICc))/(sum(exp(-0.5*out$delta.AICc), na.rm=TRUE))
    }    
  }
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

calc.mess.kstats <- function(occs.train.vals, bg.train.vals, occs.test.vals, bg.test.vals) {
  p <- rbind(occs.train.vals, bg.train.vals)
  v <- rbind(occs.test.vals, bg.test.vals)
  cat.j <- which(sapply(occs.train.vals, is.factor) == 1)
  if(length(cat.j) > 0) {
    p <- p[,-cat.j]
    v <- v[,-cat.j]
  }
  mss <- mess.vec(p, v)
  mess.quant <- quantile(mss)
  names(mess.quant) <- paste0("mess.", gsub("%", "p", names(mess.quant)))
  return(mess.quant)
}

# #' @export

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


#' @title Calculate Similarity of ENMs in Geographic Space
#' 
#' @description Compute pairwise "niche overlap" in geographic space for Maxent predictions. The value ranges from 0 (no overlap) to 1 (identical predictions).  The function uses the \code{nicheOverlap} function of the \pkg{dismo} package (Hijmans \emph{et al.} 2011).
#' 
#' @aliases calc.niche.overlap
#' @usage 
#' calc.niche.overlap(predictive.maps, overlapStat = "D", maxent.args)
#' @param preds A rasterStack of at least 2 Maxent predictive raster layers.
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
calc.niche.overlap <- function(preds, overlapStat){
  n <- raster::nlayers(preds)
  ov <- matrix(nrow = n, ncol = n)
  pb <- txtProgressBar(0, n - 1, style = 3)
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


#' @export
lookup.enm <- function(mod.name) {
  x <- switch(mod.name, 
              maxent.jar = enm.maxent.jar,
              maxnet = enm.maxnet,
              brt = enm.brt,
              bioclim = enm.bioclim)
  return(x)
}


#' @title An object of class `ENMevaluation`.
#' @description An example results file based on a call of `ENMevaluate` (see example).
#' @details The dataset is based on the simulated dataset and call of \code{\link{ENMevaluate}} shown in the example section below.
#' @format An object of class `ENMevaluation`.
#' @source Simulated data from `ENMevaluate`.
#' @examples
#' require(raster)
#' ### Simulated data environmental covariates
#' set.seed(1)
#' r1 <- raster(matrix(nrow=50, ncol=50, data=runif(10000, 0, 25)))
#' r2 <- raster(matrix(nrow=50, ncol=50, data=rep(1:100, each=100), byrow=TRUE))
#' r3 <- raster(matrix(nrow=50, ncol=50, data=rep(1:100, each=100)))
#' r4 <- raster(matrix(nrow=50, ncol=50, data=c(rep(1,1000),rep(2,500)),byrow=TRUE))
#' values(r4) <- as.factor(values(r4))
#' env <- stack(r1,r2,r3,r4)
#' 
#' ### Simulate occurrence localities
#' nocc <- 50
#' x <- (rpois(nocc, 2) + abs(rnorm(nocc)))/11
#' y <- runif(nocc, 0, .99)
#' occ <- cbind(x,y)
#' \dontrun{
#' enmeval_results <- ENMevaluate(occ, env, n.bg=500, 
#'                                partitions="block", 
#'                                categoricals="layer.4",
#'                                mod.name='maxnet', 
#'                                tune.args=list(fc = c("L","LQ","LQH","LQHP","LQHPT"), 
#'                                               rm = 1:4),
#'                                overlap=T, overlapStat="D")
#' }
#' 
#' data(enmeval_results)
#' enmeval_results
#' 
#' ### See table of evaluation metrics
#' enmeval_results@results
#' 
#' ### Plot prediction with lowest AICc
#' plot(enmeval_results@predictions[[which (enmeval_results@results$delta.AICc == 0) ]])
#' points(enmeval_results@occ.pts, pch=21, bg= enmeval_results@occ.grp)
#' 
#' ### Niche overlap statistics between model predictions
#' enmeval_results@overlap
"enmeval_results"

