#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#' @title Convert old ENMevaluation objects to new ones
#' @description Converts ENMevaluation objects made with version <=0.3.1 to
#' new ones made with version >=2.0.0.
#' @param e ENMevaluation object: the old object to convert
#' @param envs RasterStack: the original predictor variables used to generate
#' the old ENMevaluation object (these are used to make the new occs and bg slots
#' which contain the predictor variable values)
#' @note If bin.output was set to TRUE, \code{`e@results`} will be equivalent to 
#' the new results.partitions slot. Some slots are unable to be filled in because
#' previous versions of ENMeval did not record them in ENMevaluation objects:
#' variable.importance, partition.settings, other.settings, doClamp (set to TRUE
#' arbitrarily to avoid errors, but may actually have been FALSE), clamp.directions,
#' taxon.name, and rmm.
#' @importFrom rlang .data
#' @export
ENMevaluation_convert <- function(e, envs) {
  .data <- NULL
  alg <- ifelse(grepl("Maxent", e@algorithm), "maxent.jar", "maxnet")
  ts <- dplyr::distinct(e@results, fc = .data$features, rm) %>% as.data.frame()
  targs <- apply(ts, 1, function(x) paste(names(x), x, collapse = "_", sep = "."))
  ts <- cbind(ts, tune.args = targs)
  rs <- cbind(ts, e@results[,-1:-3])
  names(rs)[-1:-3] <- c("auc.train", "auc.val.avg", "auc.val.sd", "auc.diff.avg", "auc.diff.sd", 
                        "or.10p.avg", "or.10p.sd", "or.mtp.avg", "or.mtp.sd", "AICc", 
                        "delta.AICc", "w.AIC", "ncoef")
  occs <- e@occ.pts %>% dplyr::rename(lon = .data$LON, lat = .data$LAT) %>% as.data.frame()
  occs <- cbind(occs, raster::extract(envs, occs))
  bg <- e@bg.pts %>% dplyr::rename(lon = .data$LON, lat = .data$LAT) %>% as.data.frame()
  bg <- cbind(bg, raster::extract(envs, bg))
  ms <- e@models
  names(ms) <- rs$tune.args
  e_new <- ENMevaluation(algorithm = alg, tune.settings = as.data.frame(ts),
                         results = rs, results.partitions = data.frame(),
                         predictions = e@predictions, models = ms, 
                         variable.importance = list(),
                         partition.method = e@partition.method, partition.settings = list(),
                         other.settings = list(), doClamp = TRUE, clamp.directions = list(), 
                         taxon.name = "", occs = occs, occs.testing = data.frame(), 
                         occs.grp = factor(e@occ.grp), bg = bg, bg.grp = factor(e@bg.grp),
                         rmm = list())
  return(e_new)
}

#' @title Find NA cells in a RasterStack
#' @description Finds cells that are NA for at least one raster in a RasterStack.
#' @param envs RasterStack
#' 
rasStackNAs <- function(envs) {
  envs.z <- raster::values(envs)
  envs.naMismatch <- sum(apply(envs.z, 1, function(x) !all(is.na(x)) & !all(!is.na(x))))
  if(envs.naMismatch > 0) {
    message(paste0("* Found ", envs.naMismatch, " raster cells that were NA for one or more, but not all, predictor variables. Converting these cells to NA for all predictor variables."))
    envs.names <- names(envs)
    envs <- raster::stack(raster::calc(envs, fun = function(x) if(sum(is.na(x)) > 0) x * NA else x))
    names(envs) <- envs.names
  }
  return(envs)
}

#' @title Clamp predictor variables
#' @author Stephen J. Phillips, Jamie M. Kass, Gonzalo Pinilla-Buitrago
#' @param orig.vals RasterStack / matrix / data frame: environmental predictor variables (must be in same geographic 
#' projection as occurrence data), or predictor variables values for the original records
#' @param ref.vals matrix / data frame: predictor variable values for the reference records
#' (not including coordinates), used to determine the minimums and maximums -- 
#' this should ideally be the occurrences + background (can be made with raster::extract())
#' @param left character vector: names of variables to get a minimum clamp; can be "none" to turn
#' off minimum clamping
#' @param right character vector: names of variables to get a maximum clamp, can be "none" to turn
#' off maximum clamping
#' @param categoricals character vector: name or names of categorical environmental variables
#' @description This function restricts the values of one or more predictor variable rasters
#' to stay within the bounds of the input occurrence and background data (argument "ref.vals").
#' This is termed "clamping", and is mainly used to avoid making extreme extrapolations
#' when making model predictions to environmental conditions outside the range of the
#' occurrence / background data used to train the model. Clamping can be done on variables of
#' choice on one or both tails of their distributions (i.e., arguments "left" and "right" for
#' minimum and maximum clamps, respectively). If "left" and/or "right" are not specified and 
#' left at the default NULL, the function will clamp all variables for that tail (thus, the 
#' function default is to clamp all variables on both sides). To turn off clamping for one side, 
#' enter "none" for either "left" or "right".
#' 
#' Categorical variables need to be declared with the argument "categoricals". These variables
#' are excluded from the clamping analysis, but are put back into the RasterStack that is returned.
#' @return The clamped Raster* object.
#' @export

clamp.vars <- function(orig.vals, ref.vals, left = NULL, right = NULL, categoricals = NULL) {
  # find if orig.vals is a raster or not
  isRas <- inherits(orig.vals, "BasicRaster") == TRUE
  # error if "none" is included in left or right alongside variable names
  if((("none" %in% left) & length(left) > 1) | (("none" %in% right) & length(right) > 1)) {
    stop('To turn clamping off, specify the argument left, right or both of them to "none".')
  }
  # error if both sides are set to "none"
  if(!is.null(left) & !is.null(right)) {
    if(("none" %in% left) & ("none" %in% right)) {
      warning('Both left and right were set to "none", so clamping was not performed.')
      return(orig.vals)
    }
  }
  # convert all RasterStack grid cells to NA that are NA for any one raster
  # if NAs exist in an input data frame, stop the function and request user to remove manually
  if(isRas == TRUE) {
    orig.vals <- rasStackNAs(orig.vals)
  }else{
    orig.na <- sum(rowSums(is.na(orig.vals)))
    ref.na <- sum(rowSums(is.na(ref.vals)))
    if(orig.na > 0) stop("There are one or more NAs in the orig.vals table. Please remove them and rerun.")
    if(ref.na > 0) stop("There are one or more NAs in the ref.vals table. Please remove them and rerun.")  
  }
  
  # remove categorical variables from clamping analysis
  if(!is.null(categoricals)) {
    if(inherits(orig.vals, "BasicRaster") == TRUE) {
      p <- orig.vals[[-which(names(orig.vals) %in% categoricals)]]
    }else{
      p <- orig.vals[,-which(colnames(orig.vals) %in% categoricals)]
    }
    ref.vals <- ref.vals[,-which(colnames(ref.vals) %in% categoricals)]
  }else{
    p <- orig.vals
  }
  # get mins and maxs of input variable values for occs and bg
  minmaxes <- data.frame(min = apply(ref.vals, 2, min, na.rm = TRUE),
                         max = apply(ref.vals, 2, max, na.rm = TRUE))
  # function to clamp the values of the input raster
  adjustRas <- function(pp, toadjust, mm) {
    raster::stack(lapply(slot(pp, "layers"), function(oldlayer) {
      layername <- names(oldlayer)
      if (!(layername %in% toadjust)) return(oldlayer)
      newlayer <- if(mm == "min") max(oldlayer, minmaxes[layername, "min"]) else min(oldlayer, minmaxes[layername, "max"])
      names(newlayer) <- layername
      return(newlayer)
    }))
  }
  adjustDF <- function(pp, toadjust, mm) {
    inds <- which(names(pp) %in% toadjust)
    for(i in 1:ncol(pp)) {
      if(i %in% inds) {
        if(mm == "min") pp[,i] <- ifelse(pp[,i] < minmaxes[i,"min"], minmaxes[i,"min"], pp[,i]) 
        if(mm == "max") pp[,i] <- ifelse(pp[,i] > minmaxes[i,"max"], minmaxes[i,"max"], pp[,i]) 
      }
    }
    return(pp)
  }
  # default to all variables if left and/or right are NULL
  if(is.null(left)) left <- names(p)
  if(is.null(right)) right <- names(p)
  
  f <- ifelse(inherits(orig.vals, "BasicRaster") == TRUE, adjustRas, adjustDF)
  # clamp both sides unless left or right is "none"
  if("none" %in% left) {
    out <- f(p, right, "max")
  }else if("none" %in% right) {
    out <- f(p, left, "min")
  }else{
    out <- f(f(p, left, "min"), right, "max")  
  }
  # add back the categorical variable to the stack
  if(!is.null(categoricals)) {
    if(inherits(orig.vals, "BasicRaster") == TRUE) {
      out <- raster::addLayer(out, orig.vals[[categoricals]])
      out <- out[[names(orig.vals)]]
    }else{
      for(i in 1:length(categoricals)) {
        out <- cbind(out, orig.vals[,categoricals[i]])
        names(out)[ncol(out)] <- categoricals[i]
      }
      out <- out[,names(orig.vals)]
    }
  }
  return(out)
}


#' @title Calculate AICc from Maxent model prediction
#' @description This function calculates AICc for Maxent models based on Warren 
#' and Seifert (2011).
#' @param p.occs data frame: raw (maxent.jar) or exponential (maxnet) predictions for
#' the occurrence localities based on one or more models
#' @param ncoefs numeric: number of non-zero model coefficients
#' @param p RasterStack: raw (maxent.jar) or exponential (maxnet) model predictions;
#' if NULL, AICc will be calculated based on the background points, which already
#' have predictions that sum to 1 and thus need no correction -- this assumes that
#' the background points represent a good sample of the study extent 
#' @return data frame with three columns:
#' \code{AICc} is the Akaike Information Criterion corrected for small sample
#' sizes calculated as:
#' \deqn{ (2 * K - 2 * logLikelihood) + (2 * K) * (K+1) / (n - K - 1)}
#' where \emph{K} is the number of non-zero coefficients in the model and \emph{n} is the number of
#' occurrence localities.  The \emph{logLikelihood} is calculated as:
#' \deqn{ sum(log(vals))}
#' where \emph{vals} is a vector of Maxent raw/exponential values at occurrence localities
#' and the sum of these values across the study extent is equal to 1.
#' \code{delta.AICc} is the difference between the AICc of a given model and
#' the AICc of the model with the lowest AICc.
#' \code{w.AICc} is the Akaike weight (calculated as the relative likelihood of
#' a model (exp(-0.5 * \code{delta.AICc})) divided by the sum of the likelihood
#' values of all models included in a run.  These can be used for model
#' averaging (Burnham and Anderson 2002).
#' @aliases calc.aicc
#' @details As motivated by Warren and Seifert (2011) and implemented in ENMTools (Warren 
#' \emph{et al.} 2010), this function calculates the small sample size version of Akaike 
#' Information Criterion for ENMs (Akaike 1974). We use AICc (instead of AIC) regardless of 
#' sample size based on the recommendation of Burnham and Anderson (1998, 2004).  The number of 
#' coefficients is determined by counting the number of non-zero coefficients in the 
#' \code{maxent} lambda file (\code{m@lambdas} for maxent.jar and \code{m$betas} for maxnet.  
#' See Warren \emph{et al.} (2014) for limitations of this approach, namely that the number of 
#' non-zero coefficients is an estimate of the true degrees of freedom. For Maxent ENMs, AICc 
#' is calculated by first standardizing the raw predictions such that all cells in the study 
#' extent sum to 1, then extracting the occurrence record predictions. The predictions of the
#' study extent may not sum to 1 if the background does not cover every grid cell -- as the 
#' background predictions sum to 1 by definition, extra predictions for grid cells not in 
#' the training data will add to this sum. When no raster data is provided, the raw predictions 
#' of the occurrence records are used to calculate AICc without standardization, with the 
#' assumption that the background records have adequately represented the occurrence records. 
#' The standardization is not necessary here because the background predictions sum to 1 
#' already, and the occurrence data is a subset of the background. This will not be true if 
#' the background does not adequately represent the occurrence records, in which case the 
#' occurrences are not a subset of the background and the raster approach should be used 
#' instead. The likelihood of the data for a given model is then calculated by taking the 
#' product of the raw occurrence predictions (Warren and Seifert 2011), or the sum of their 
#' logs, as is implemented here.
#' 
#' @seealso \code{maxent} in the \pkg{dismo} package.
#' 
#' @note Returns all \code{NA}s if the number of non-zero coefficients is larger than the
#' number of observations (occurrence localities).
#' 
#' @references 
#' Akaike, H. (1974) A new look at the statistical model identification. \emph{IEEE Transactions on Automatic Control}, \bold{19}: 716-723. \url{https://doi.org/10.1109/TAC.1974.1100705}
#' 
#' Burnham, K. P. and Anderson, D. R. (1998) Model selection and multimodel inference: a practical information-theoretic approach. Springer, New York.
#' 
#' Burnham, K. P. and Anderson, D. R. (2004) Multimodel inference: understanding AIC and BIC in model selection. \emph{Sociological Methods and Research}, \bold{33}: 261-304. \url{https://doi.org/10.1177/0049124104268644}
#' 
#' Warren, D. L., Glor, R. E, and Turelli, M. (2010) ENMTools: a toolbox for comparative studies of environmental niche models. \emph{Ecography}, \bold{33}: 607-611. \url{https://doi.org/10.1111/j.1600-0587.2009.06142.x}
#' 
#' Warren, D. L., & Seifert, S. N. (2011). Ecological niche modeling in Maxent: the importance of model complexity and the performance of model selection criteria. \emph{Ecological Applications}, \bold{21}: 335-342. \url{https://doi.org/10.1890/10-1171.1}
#' 
#' Warren, D. L., Wright, A. N., Seifert, S. N., and Shaffer, H. B. (2014). Incorporating model complexity and sampling bias into ecological niche models of climate change risks faced by 90 California vertebrate species of concern. \emph{Diversity and Distributions}, \bold{20}: 334-343. \url{https://doi.org/10.1111/ddi.12160}
#' 
#' @export
aic.maxent <- function(p.occs, ncoefs, p = NULL) {
  # differential behavior for summing if p is Raster* or data frame
  if(!is.null(p)) {
    p.sum <- raster::cellStats(p, sum)  
    # if total does not sum to 1 (this happens when the background does not fully
    # cover the study extent), standardize the study extent predictions so that they sum to 1
    # and use these corrected occurrence predictions as likelihoods
    # (dividing the occurrence predictions by the sum of the study extent predictions
    # achieves the above)
    for(i in 1:raster::nlayers(p)) if(p.sum[i] != 1) p.occs[,i] <- p.occs[,i] / p.sum[i]
  }
  # if more model coefficients than data points, determine AIC to be invalid:
  # this avoids considering overly complex models at all
  n.occs <- nrow(p.occs)
  AIC.valid <- ncoefs < n.occs
  for(i in 1:length(AIC.valid)) {
    if(AIC.valid[i] == FALSE) {
      message(paste("Warning: model", names(AIC.valid)[i], "has more non-zero coefficients (ncoef) than occurrence records for training, so AIC cannot be calculated."))
    }
  }
  # calculate log likelihood
  LL <- colSums(log(p.occs), na.rm = TRUE)
  AICc <- (2 * ncoefs - 2 * LL) + (2 * (ncoefs) * (ncoefs + 1) / (n.occs - ncoefs - 1))
  # if determined invalid or if infinite, make AICc NA
  AICc <- sapply(1:length(AICc), function(x) ifelse(AIC.valid[x] == FALSE | is.infinite(AICc[x]), NA, AICc[x]))
  # make output table
  out <- data.frame(AICc = AICc, delta.AICc = AICc - min(AICc, na.rm=TRUE), row.names = NULL)
  out$w.AIC <- exp(-0.5 * out$delta.AICc) / sum(exp(-0.5 * out$delta.AICc), na.rm=TRUE)
  return(out)
}

#' @title Corrected variance function
#' @description Calculate variance corrected for non-independence of \emph{k}-fold iterations
#'
#' @details `corrected.var` calculates variance corrected for non-independence of \emph{k}-fold iterations.  
#' See Appendix of Shcheglovitova & Anderson (2013) and other references (Miller 1974; Parr 1985; Shao and Wu 1989) for additional details. 
#' This function calculates variance that is corrected for the non-independence of \emph{k} cross-validation iterations.  
#' Following Shao and Wu (1989): 
#' 
#' \deqn{Sum Of Squares * ((n-1)/n)} 
#' 
#' where \emph{n} = the number of \emph{k}-fold iterations.
#'
#' @param x numeric vector: input values
#' @param nk numeric: number of \emph{k}-fold iterations
#' @return A numeric value of the corrected variance.
#' @author Robert Muscarella <bob.muscarella@gmail.com>
#' @references 
#' Miller, R. G. (1974) The jackknife - a review. \emph{Biometrika}, \bold{61}: 1-15. \url{https://doi.org/10.1093/biomet/61.1.1}
#'   
#' Parr, W. C. (1985) Jackknifing differentiable statistical functionals. \emph{Journal of the Royal Statistics Society, Series B}, \bold{47}: 56-66. \url{https://doi.org/10.1111/j.2517-6161.1985.tb01330.x}
#'   
#' Shao J. and Wu, C. F. J. (1989) A general theory for jackknife variance estimation. \emph{Annals of Statistics}, \bold{17}: 1176-1197. \url{https://doi.org/10.1214/aos/1176347263}
#'   
#' Shcheglovitova, M. and Anderson, R. P. (2013) Estimating optimal complexity for ecological niche models: a jackknife approach for species with small sample sizes. \emph{Ecological Modelling}, \bold{269}: 9-17. \url{https://doi.org/10.1016/j.ecolmodel.2013.08.011}
#' @export
corrected.var <- function(x, nk){
  sum((x - mean(x))^2) * ((nk-1)/nk)
}

#' @title Calculate 10 percentile threshold
#' @description Function to calculate the 10 percentile threshold from training predictions
#' @param pred.train numeric vector: training occurrence predictions
#' 
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

#' @title Calculate Similarity of ENMs in Geographic Space
#' @description Compute pairwise "niche overlap" in geographic space for Maxent predictions. The value ranges from 0 (no overlap) to 1 (identical predictions).  The function uses the \code{nicheOverlap} function of the \pkg{dismo} package (Hijmans \emph{et al.} 2011).
#' @param predictors RasterStack: at least 2 Maxent raster predictions
#' @param overlapStat character: either "D" or "I", the statistic calculated by the \code{nicheOverlap} function of the \pkg{dismo} package (default: "D")
#' @param quiet boolean: if TRUE, silence all function messages (but not errors)
#' @details "D" refers to Schoeners \emph{D} (Schoener 1968), while "I" refers to the \emph{I} similarity statistic from Warren \emph{et al.} (2008).
#' @return 
#' A matrix with the lower triangle giving values of pairwise "niche overlap" in geographic space.  Row and column names correspond to the results table output by \code{\link{ENMevaluate}()}.
#' @references 
#' Hijmans, R. J., Phillips, S., Leathwick, J. & Elith, J. (2011) dismo package for R. Available online at: \url{https://cran.r-project.org/package=dismo}.
#' 
#' Schoener, T. W. (1968) The \emph{Anolis} lizards of Bimini: resource partitioning in a complex fauna. \emph{Ecology}, \bold{49}: 704-726. \url{https://doi.org/10.2307/1935534}
#' 
#' Warren, D. L., Glor, R. E., Turelli, M. & Funk, D. (2008) Environmental niche equivalency versus conservatism: quantitative approaches to niche evolution. \emph{Evolution}, \bold{62}: 2868-2883. \url{https://doi.org/10.1111/j.1558-5646.2008.00482.x}
#' 
#' @author 
#' Based on \pkg{dismo}::\code{nicheOverlap}, which is based on \pkg{SDMTools}::\code{Istat}
#' Robert Muscarella <bob.muscarella@gmail.com>
#' @seealso 
#' `nicheOverlap` in the \pkg{dismo} package

#' @export
calc.niche.overlap <- function(predictors, overlapStat, quiet=FALSE){
  n <- raster::nlayers(predictors)
  ov <- matrix(nrow = n, ncol = n)
  if(quiet != TRUE) pb <- txtProgressBar(0, n - 1, style = 3)
  for(i in 1:(n - 1)){
    if(quiet != TRUE) setTxtProgressBar(pb, i)
    for(j in (i + 1):n){
      ov[j, i] <- dismo::nicheOverlap(predictors[[i]], predictors[[j]], stat = overlapStat)
    }
  }
  colnames(ov) <- names(predictors)
  rownames(ov) <- names(predictors)
  if(quiet != TRUE) close(pb)
  return(ov)
}

#' @title Look up ENMdetails abject
#' @description Internal function to look up ENMdetails objects.
#' @param algorithm character: algorithm name (must be implemented as ENMdetails object)
#' 
lookup.enm <- function(algorithm) {
  x <- switch(algorithm, 
              maxent.jar = enm.maxent.jar,
              maxnet = enm.maxnet,
              # randomForest = enm.randomForest,
              # boostedRegressionTrees = enm.boostedRegressionTrees,
              bioclim = enm.bioclim
  )
  return(x)
}


#' @title Look up version of maxent.jar
#' @description Internal function to look up the version of the maxent.jar being used.
#' 
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
