#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#' @title User-specified evaluation function for model validation
#' @name user.eval
#' @usage my.user.eval <- myFunc(vars)
#' ENMevaluate(..., user.eval = my.user.eval)
#' @description This is a custom function for specifying performance metrics not included in ENMeval.
#' The function must be first defined and then input as the argument "user.eval" into ENMevaluate(). 
#' This function has a single argument called "vars", which is a list that includes different data 
#' that can be used to calculate the metric. Below is the list of data included in vars, which can
#' be accessed with $, as in vars$occs.train.z. See the vignette for a worked example.
#' @param enm ENMdetails object
#' @param occs.train.z data frame: predictor variable values for training occurrences
#' @param occs.val.z data frame: predictor variable values for validation occurrences
#' @param bg.train.z data frame: predictor variable values for training background
#' @param bg.val.z data frame: predictor variable values for validation background
#' @param mod.k Model object for current partition (k)
#' @param nk numeric: number of folds (i.e., partitions)
#' @param other.settings named list: other settings specified in ENMevaluate()
#' @param partitions character: name of the partition method (e.g., "block")
#' @param occs.train.pred numeric: predictions made by mod.k for training occurrences
#' @param occs.val.pred numeric: predictions made by mod.k for validation occurrences
#' @param bg.train.pred numeric: predictions made by mod.k for training background
#' @param bg.val.pred numeric: predictions made by mod.k for validation background
NULL

#' @title Partition settings
#' @name partition.settings
#' @usage ps <- list(kfolds = 5)
#' ps <- list(orientation = "lat_lat")
#' ps <- list(aggregation.factor = c(4,4))
#' ENMevaluate(..., partition.settings = ps)
#' @description This is a named list used to specify certain settings for partitioning schema.
#' It is inserted as an argument to ENMevaluate(). Some partitioning schema require a setting to
#' be specified. It is helpful to use the evalplot.grps() plotting function to visualize differences in settings.
#' @param orientation character: one of "lat_lon", "lon_lat", "lat_lat", or "lon_lon" (required for block partition)
#' @param aggregation.factor numeric vector: one or two numbers specifying
#' the factor with which to aggregate the envs raster to assign partitions (required for the checkerboard partitions)
#' @param kfolds Required for the random partition. This specifies the number of random partition groups, or folds, to make.
#' @details For the block partition, the orientation specifications are abbreviations for "latitude" and "longitude", and 
#' they determine the order and orientations with which the block partitioning function creates the partition groups. For example,
#' "lat_lon" will split the occurrence localities first by latitude, then by longitude.
#' 
#' For the checkerboard partitions, the aggregation factor specifies how much to aggregate the existing cells in the envs raster
#' to make new spatial partitions. For example, checkerboard1 with an aggregation factor value of 2 will make the grid cells 
#' 4 times larger and then assign occurrence and background records to partition groups based on which cell they are in. 
#' The checkerboard2 partition is hierarchical, so cells are first aggregated to define groups like checkerboard1, but a 
#' second aggregation is then made to separate the resulting 2 bins into 4 bins. For checkerboard2, two different numbers can be used
#' to specify the two levels of the hierarchy, or if a single number is inserted, that value will be used for both levels.
NULL

#' @title Other settings
#' @name other.settings
#' @usage # an example of specifying settings different from default
#' os <- list(pred.type = "logistic", 
#'   abs.auc.diff = FALSE, 
#'   validation.bg = "partition")
#' ENMevaluate(..., other.settings = os)
#' @description This is a named list used to specify extra settings for the analysis.
#' It is inserted as an argument to ENMevaluate(). All of these settings have internal defaults,
#' so if they are not specified the analysis will be run with default settings.
#' @param abs.auc.diff boolean: if TRUE, take absolute value of AUCdiff (default: TRUE)
#' @param validation.bg character: either "full" to calculate AUC and CBI with respect to the full background, or
#' "partition" to calculate them with respect to the validation partition background (default: "full")
#' @param pred.type character: specifies which prediction type should be used to generate maxnet or maxent.jar prediction rasters (default: "cloglog")
#' @param other.args named list: any additional model arguments not specified for tuning
NULL

#' @title Clamp predictor variables
#' @author Stephen J. Phillips, Jamie M. Kass, Gonzalo Pinilla-Buitrago
#' @param predictors RasterStack: environmental predictor variables (must be in same geographic projection as occurrence data)
#' @param p.z matrix / data frame: predictor variable values for the reference records
#' (not including coordinates), used to determine the minimums and maximums -- 
#' this should ideally be the occurrences + background (can be made with raster::extract())
#' @param left character vector: names of variables to get a minimum clamp; can be "none" to turn
#' off minimum clamping
#' @param right character vector: names of variables to get a maximum clamp, can be "none" to turn
#' off maximum clamping
#' @param categoricals character vector: name or names of categorical environmental variables
#' @description This function restricts the values of one or more predictor variable rasters
#' to stay within the bounds of the input occurrence and background data (argument "p.z").
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

clamp.vars <- function(predictors, p.z, left = NULL, right = NULL, categoricals = NULL) {
  if((("none" %in% left) & length(left) > 1) | (("none" %in% right) & length(right) > 1)) {
    stop('To turn clamping off, specify the argument left, right or both of them to "none".')
  }
  
  if(!is.null(left) & !is.null(right)) {
    if(("none" %in% left) & ("none" %in% right)) {
      warning('Both left and right were set to "none", so clamping was not performed.')
      return(predictors)
    }
  }
  # remove categorical variables from clamping analysis
  if(!is.null(categoricals)) {
    p <- predictors[[-which(names(predictors) == categoricals)]]
    p.z <- p.z[,-which(colnames(p.z) == categoricals)]
  }else{
    p <- predictors
  }
  # get mins and maxs of input variable values for occs and bg
  minmaxes <- data.frame(min = apply(p.z, 2, min, na.rm = TRUE),
                         max = apply(p.z, 2, max, na.rm = TRUE))
  # function to clamp the values of the input raster
  adjust <- function(pp, toadjust, mm) {
    raster::stack(lapply(slot(pp, "layers"), function(oldlayer) {
      layername <- names(oldlayer)
      if (!(layername %in% toadjust)) return(oldlayer)
      newlayer <- if(mm == "min") max(oldlayer, minmaxes[layername, "min"]) else min(oldlayer, minmaxes[layername, "max"])
      names(newlayer) <- layername
      return(newlayer)
    }))}
  # default to all variables if left and/or right are NULL
  if(is.null(left)) left <- names(p)
  if(is.null(right)) right <- names(p)
  
  # clamp both sides unless left or right is "none"
  if("none" %in% left) {
    out <- adjust(p, right, "max")
  }else if("none" %in% right) {
    out <- adjust(p, left, "min")
  }else{
    out <- adjust(adjust(p, left, "min"), right, "max")  
  }
  # add back the categorical variable to the stack
  if(!is.null(categoricals)) out <- raster::addLayer(out, predictors[[categoricals]])
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
#' coefficients is determined by counting the number of non-zero coefficients in 
#' the \code{maxent} lambda file (\code{m@lambdas} for maxent.jar and \code{m$betas} for maxnet.  
#' See Warren \emph{et al.} (2014) for limitations of this approach, namely that the number of non-zero coefficients is an 
#' estimate of the true degrees of freedom.  For Maxent ENMs, AICc is 
#' calculated by standardizing the raw output such that all cells in the study 
#' extent sum to 1. The likelihood of the data for a given model is then 
#' calculated by taking the product of the raw output values (or the sum of their logs, as is implemented here)
#' for all grid cells that contain an occurrence locality (Warren and Seifert 2011).
#' @seealso \code{maxent} in the \pkg{dismo} package.
#' @note Returns all \code{NA}s if the number of non-zero coefficients is larger than the
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
aic.maxent <- function(p.occs, ncoefs, p = NULL) {
  # differential behavior for summing if p is Raster* or data frame
  if(!is.null(p)) {
    p.sum <- raster::cellStats(p, sum)  
    # if total does not sum to 1, standardize so that the sum is 1
    for(i in 1:raster::nlayers(p)) if(p.sum[i] != 1) p.occs[,i] <- p.occs[,i] / p.sum[i]
  }
  # if more model coefficients than data points, determine AIC to be invalid:
  # this avoids considering overly complex models at all
  n.occs <- nrow(p.occs)
  AIC.valid <- ncoefs < n.occs
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
#' @param x numeric vector: input values
#' @param nk numeric: number of \emph{k}-fold iterations
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

#' @title Calculate Similarity of ENMs in Geographic Space
#' 
#' @description Compute pairwise "niche overlap" in geographic space for Maxent predictions. The value ranges from 0 (no overlap) to 1 (identical predictions).  The function uses the \code{nicheOverlap} function of the \pkg{dismo} package (Hijmans \emph{et al.} 2011).
#' 
#' @aliases calc.niche.overlap
#' @usage 
#' calc.niche.overlap(predictive.maps, overlapStat = "D", maxent.args)
#' @param predictors RasterStack: at least 2 Maxent raster predictions
#' @param overlapStat character: either "D" or "I", the statistic calculated by the \code{nicheOverlap} function of the \pkg{dismo} package (default: "D")
#' @details "D" refers to Schoeners \emph{D} (Schoener 1968), while "I" refers to the \emph{I} similarity statistic from Warren \emph{et al.} (2008).
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

# function to look up the corresponding ENMdetails abject
#' @export
lookup.enm <- function(algorithm) {
  x <- switch(algorithm, 
              maxent.jar = enm.maxent.jar,
              maxnet = enm.maxnet,
              randomForest = enm.randomForest,
              boostedRegressionTrees = enm.boostedRegressionTrees,
              bioclim = enm.bioclim
  )
  return(x)
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
