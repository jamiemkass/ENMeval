#' Tuning and evaluation of ecological niche models
#' 
#' @param occs matrix or data frame with two columns for longitude and latitude 
#' of occurrence localities, in that order
#' @param envs Raster* object of environmental variables (must be in 
#' same geographic projection as occurrence data)
#' @param bg matrix or data frame with two columns for longitude and latitude of 
#' background (or pseudo-absence) localities, in that order; if NULL, points 
#' will be randomly sampled across \code{envs} with the number specified by 
#' parameter \code{n.bg}
#' @param occs.vals matrix or data frame of environmental values corresponding
#' to occurrence localities, intended to be input when environmental rasters
#' are not used (\code{envs} is NULL) 
#' @param bg.vals matrix or data frame of environmental values corresponding
#' to background (or pseudo-absence) localities, intended to be input when 
#' environmental rasters are not used (\code{envs} is NULL) 
#' @param mod.fun function of chosen model
#' @param tune.args named list of model settings to be tuned
#' @param other.args named list of any additional model arguments not specified 
#' for tuning
#' @param categoricals character vector of names of categorical 
#' environmental variables
#' @param partitions character of name of partitioning technique (see
#' \code{?partitions})
#' @param occs.folds numeric vector of partition group (fold) for each
#' occurrence locality, intended for user-defined partitions
#' @param bg.folds numeric vector of partition group (fold) for each background 
#' (or pseudo-absence) locality, intended for user-defined partitions
#' @param occs.ind matrix or data frame with two columns for longitude and latitude 
#' of occurrence localities, in that order, intended for independent evaluation;
#' when \code{partitions = "independent"}; these occurrences will be used only 
#' for evaluation, and not for model training, and thus no cross validation will 
#' be done
#' @param kfolds numeric for number of partition groups (folds), only for random
#' k-fold partitioning
#' @param aggregation.factor numeric vector with length 2 that specifies the
#' factors for aggregating \code{envs} in order to perform checkerboard
#' partitioning
#' @param n.bg numeric for number of random background (or pseudo-absence) points
#' to sample; necessary if \code{bg} is NULL
#' @param overlap boolean (TRUE or FALSE); if TRUE, calculate niche overlap 
#' statistics
#' @param doClamp boolean (TRUE or FALSE); if TRUE, clamp model responses; only
#' applicable for Maxent models
#' @param skipRasters boolean (TRUE or FALSE); if TRUE, skip raster predictions
#' @param abs.auc.diff boolean (TRUE or FALSE); if TRUE, take absolute value of
#' AUCdiff; default is TRUE
#' @param parallel boolean (TRUE or FALSE); if TRUE, run with parallel processing
#' @param numCores boolean (TRUE or FALSE); if TRUE, use specifed number of cores
#' for parallel processing
#' @param updateProgress boolean (TRUE or FALSE); if TRUE, use shiny progress
#' bar; only for use in shiny apps
#'
#' @return 
#'
#' @examples
#'
#' @importFrom magrittr %>%
#' @export 

ENMevaluate <- function(occs, envs = NULL, bg = NULL, occs.vals = NULL, bg.vals = NULL, 
                        mod.fun, tune.args, other.args = NULL, categoricals = NULL, 
                        partitions = NULL, occs.folds = NULL, bg.folds = NULL, occs.ind = NULL, 
                        kfolds = NA, aggregation.factor = c(2, 2), n.bg = 10000, overlap = FALSE,   
                        doClamp = TRUE, skipRasters = FALSE, abs.auc.diff = TRUE, parallel = FALSE, 
                        numCores = NULL, updateProgress = FALSE) {
  
  # record start time
  start.time <- proc.time()
  # get model function's name
  mod.name <- as.character(substitute(mod.fun))[3]
  
  ########### #
  # CHECKS ####
  ########### #
  
  ## general parameter checks
  if(!(partitions %in% c("jackknife", "randomkfold", "block", "checkerboard1", 
                         "checkerboard2", "user", "independent", "none"))) {
    stop("Please make sure partition method is one of the available options.")
  }
  
  # print model-specific message
  model.msgs(tune.args, mod.name)
  
  # make table for all tuning parameter combinations
  tune.tbl <- expand.grid(tune.args, stringsAsFactors = FALSE)
  
  ## data checks and formatting
  # if occs is combined occurrence and background with environmental
  # predictor values (SWD format)
  if(is.null(envs)) {
    if(is.null(occs.vals)) {
      stop("If inputting data without rasters (SWD), please specify both 
           occs.vals and bg.vals.")
    }
    occs.vals <- as.data.frame(occs.vals)
    if(is.null(bg.vals)) {
      stop("If inputting data without rasters (SWD), please specify both 
           occs.vals and bg.vals.")
    }
    bg.vals <- as.data.frame(bg.vals)
    if(is.null(bg)) {
      stop("If inputting data without rasters (SWD), please specify bg in 
           addition to occs.")
    }
    message("Data without rasters were input (SWD format), so no raster 
            predictions will be generated and AICc cannot be calculated.")
  }else{
    # if no background points specified, generate random ones
    if(is.null(bg)) {
      bg <- dismo::randomPoints(envs, n = n.bg)
    }  
  }
  
  # make sure occs and bg are data frames with identical column names
  occs <- data.frame(longitude = occs[,1], latitude = occs[,2])
  bg <- data.frame(longitude = bg[,1], latitude = bg[,2])
  
  # print message for selected cross-validation method
  # for occs.ind settings, partitions should be NULL
  
  if(partitions == "jackknife") {
    folds <- get.jackknife(occs, bg)
    message("Doing model evaluations with k-1 jackknife cross validation...")
  }
  if(partitions == "randomkfold") {
    folds <- get.randomkfold(occs, bg, kfolds)
    message(paste0("Doing model evaluations with random", kfolds, "-fold cross validation..."))
  }
  if(partitions == "block") {
    folds <- get.block(occs, bg)
    message("Doing model evaluations with spatial block (4-fold) cross validation...")
  }
  if(partitions == "checkerboard1") {
    if(is.null(envs)) stop("Cannot use checkerboard partitioning if envs is NULL.")
    folds <- get.checkerboard1(occs, envs, bg, aggregation.factor)
    message("Doing model evaluations with checkerboard (2-fold) cross validation...")
  }
  if(partitions == "checkerboard2") {
    if(is.null(envs)) stop("Cannot use checkerboard partitioning if envs is NULL.")
    folds <- get.checkerboard2(occs, envs, bg, aggregation.factor)
    message("Doing model evaluations with hierarchical checkerboard (4-fold) cross validation...")
  }
  if(partitions == "user") {
    folds <- list(occs.folds = occs.folds, bg.folds = bg.folds)
    userk <- length(unique(occs.folds))
    message(paste0("Doing model evaluations with user-defined ", userk, "-fold cross validation..."))
  }
  if(partitions == "independent") {
    if(is.null(occs.ind)) stop("Cannot use independent partitioning if occs.ind is NULL.")
    folds <- NULL
    message("Doing model evaluations with independent testing data...")
  }
  if(partitions == "none") {
    folds <- NULL
    message("Skipping model evaluations (only calculating AICc)...")
  }
  
  # unpack occs.folds and bg.folds
  occs.folds <- folds$occs.folds
  bg.folds <- folds$bg.folds
  
  ############# #
  # ANALYSIS ####
  ############# #
  
  # extract predictor variable values at coordinates for occs and bg
  if(!is.null(envs)) {
    occs.vals <- as.data.frame(raster::extract(envs, occs))
    bg.vals <- as.data.frame(raster::extract(envs, bg))  
  }
  
  # remove rows from occs, occs.vals, occs.folds with NA for any predictor variable
  occs.vals.na <- sum(rowSums(is.na(occs.vals)) > 0)
  bg.vals.na <- sum(rowSums(is.na(bg.vals)) > 0)
  if(occs.vals.na > 0) {
    i <- !apply(occs.vals, 1, anyNA)
    occs <- occs[i,]
    occs.folds <- occs.folds[i]
    occs.vals <- occs.vals[i,]
    message(paste("There were", occs.vals.na, "occurrence records with NA for at least one predictor variable. Removed these from analysis, resulting in", nrow(occs.vals), "occurrence records."))
  }
  # do the same for associated bg variables
  if(bg.vals.na > 0) {
    i <- !apply(bg.vals, 1, anyNA)
    occs <- occs[i,]
    bg.folds <- bg.folds[i]
    bg.vals <- bg.vals[i,]
    message(paste("There were", bg.vals.na, "background records with NA for at least one predictor variable. Removed these from analysis, resulting in", nrow(bg.vals), "background records."))
  }
  
  # convert fields for categorical data to factor class
  if(!is.null(categoricals)) {
    for(i in 1:length(categoricals)) {
      occs.vals[, categoricals[i]] <- as.factor(occs.vals[, categoricals[i]])
      bg.vals[, categoricals[i]] <- as.factor(bg.vals[, categoricals[i]])
    }
  }
  
  if(parallel == TRUE) {
    results <- tune.enm.parallel(occs.vals, bg.vals, occs.folds, bg.folds, envs, 
                                 mod.fun, mod.name, partitions, tune.tbl, other.args, categoricals, 
                                 occs.ind, doClamp, skipRasters, abs.auc.diff)
  } 
  else {
    results <- tune.enm(occs.vals, bg.vals, occs.folds, bg.folds, envs, 
                        mod.fun, mod.name, partitions, tune.tbl, other.args, categoricals, 
                        occs.ind, doClamp, skipRasters, abs.auc.diff, updateProgress)
  }
  
  res <- collateResults(results, tune.tbl, envs, mod.name, partitions, skipRasters)
  if(is.null(occs.folds)) occs.folds <- 0
  if(is.null(bg.folds)) bg.folds <- 0
  e <- ENMevaluation(algorithm = mod.name, 
                     stats = res$stats, kstats = res$kstats,
                     predictions = res$preds, models = res$mods, 
                     partition.method = partitions,
                     occs = occs, occs.folds = occs.folds,
                     bg = bg, bg.folds = bg.folds)
  
  
  # if niche overlap selected, calculate and add the resulting matrix to results
  if (overlap == TRUE) {
    if (length(args) > 1) {
      if(nlayers(results@predictions) > 1) {
        message("Calculating niche overlap")
        overlap.mat <- calc.niche.overlap(results@predictions, "D")
        results@overlap <- overlap.mat
      }
      else {
        message("Cannot calculate niche overlap without rasterPreds")
      }
    }
    else {
      message("Need >1 settings to do niche overlap")
    }
  }
  
  # calculate time expended and print message
  timed <- proc.time() - start.time
  t.min <- floor(timed[3] / 60)
  t.sec <- timed[3] - (t.min * 60)
  message(paste("ENMeval completed in", t.min, "minutes", round(t.sec, 1), "seconds."))
  
  return(e)
}
