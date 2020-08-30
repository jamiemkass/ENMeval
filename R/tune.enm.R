#' @title Iterate tuning of ENMs
#' @description Internal functions to tune and summarize results for ecological niche models (ENMs) iteratively across a range of user-specified tuning settings. 
#' Function \code{tune.parallel()} tunes ENMs with parallelization. Function \code{cv.enm()} calculates training and validation evaluation statistics for one set of specified tuning parameters.
#' @aliases tune.parallel tune.regular cv.enm
#' @param d data frame from \code{ENMevaluate()} with occurrence and background coordinates (or coordinates plus predictor variable values) and partition group values
#' @param envs Raster* object of environmental variables (must be in same geographic projection as occurrence data)
#' @param enm Object of class \link{ENMdetails}.
#' @param partitions character of name of partitioning technique (see \code{?partitions})
#' @param tune.tbl Data frame of tuning parameter combinations.
#' @param other.settings list of settings from \code{ENMevaluate()} containing other.args, clamp, pred.type, abs.auc.diff, cbi.cv
#' @param numCores boolean (TRUE or FALSE); if TRUE, use specifed number of cores for parallel processing
#' @param parallelType character of either "doParallel" or "doSNOW" to conduct parallelization

#' @name tune.enm
NULL

#' @rdname tune.enm
tune.train <- function(enm, occs.z, bg.z, mod.full, mod.full.pred, envs, other.settings, partitions, quiet) {
  # get model predictions for training data
  occs.pred <- enm@predict(mod.full, occs.z, other.settings)
  bg.pred <- enm@predict(mod.full, bg.z, other.settings)
  
  # training AUC
  e <- dismo::evaluate(occs.pred, bg.pred)
  auc.train <- e@auc
  
  # calculate training CBI with background, not raster -- this is in order
  # to calculate CBI in a standardized way across training and validation data for 
  # spatial partitions (ENMeval does not mask rasters to partition areas and hence 
  # does not have this information)
  # combine occs and bg predictions as input for "fit" because the interval
  # is determined from "fit" only, and if occs have higher predictions than
  # bg, the interval will be cut short
  cbi.train <- ecospat::ecospat.boyce(c(bg.pred, occs.pred), occs.pred, PEplot = FALSE)$Spearman.cor
  
  out.df <- data.frame(auc.train = auc.train, cbi.train = cbi.train)
  return(out.df)
}

#' @rdname tune.enm
#' @description Validation CBI is calculated here with background values, not raster data, in order
# to standardize the methodology for both training and validation data for 
# spatial partitions, as ENMeval does not mask rasters to partition areas and hence 
# does not have partitioned raster data. Further, predictions for occurrence and background localities 
# are combined as input for the parameter "fit" in \code{ecospat::ecospat_boyce()} because the interval
# is determined from "fit" only, and if test occurrences all have higher predictions than the background, 
# the interval will be cut short.
tune.validate <- function(enm, occs.train.z, occs.val.z, bg.train.z, bg.val.z, mod.k, nk, other.settings, partitions, user.eval, quiet) {
  # get model predictions for training and validation data for partition k
  occs.train.pred <- enm@predict(mod.k, occs.train.z, other.settings)
  occs.val.pred <- enm@predict(mod.k, occs.val.z, other.settings)
  bg.train.pred <- enm@predict(mod.k, bg.train.z, other.settings)
  # if the background was partitioned, get predictions for this subset; if not, it will be NULL
  if(nrow(bg.val.z) > 0) {
    bg.val.pred <- enm@predict(mod.k, bg.val.z, other.settings)
  }else{
    bg.val.pred <- NULL
    # if(quiet != TRUE) message("\nNOTE: Background points were not partitioned for this analysis, so model validation will proceed on the full background.")
  }
  
  # if validation.bg == "full", calculate training and validation AUC and CBI based on the full background 
  # (training + validation background for all statistics) 
  # see Radosavljevic & Anderson 2014 for an example of calculating validation AUC with spatial partitions over a shared background
  if(other.settings$validation.bg == "full") {
    e.train <- dismo::evaluate(occs.train.pred, c(bg.train.pred, bg.val.pred))
    auc.train <- e.train@auc
    e.val <- dismo::evaluate(occs.val.pred, c(bg.train.pred, bg.val.pred))
    auc.val <- e.val@auc
    # calculate AUC diff as training AUC minus validation AUC with a shared background
    auc.diff <- auc.train - auc.val
    # calculate CBI based on the full background (do not calculate for jackknife partitions)
    if(partitions != "jackknife") {
      cbi.val <- ecospat::ecospat.boyce(c(bg.train.pred, bg.val.pred, occs.val.pred), occs.val.pred, PEplot = FALSE)$Spearman.cor
    }else{
      cbi.val <- NA
    }
  # if validation.bg == "partition", calculate training and validation AUC and CBI based on the partitioned backgrounds only 
  # (training background for training statistics and validation background for validation statistics) 
  }else if(other.settings$validation.bg == "partition") {
    e.train <- dismo::evaluate(occs.train.pred, bg.train.pred)
    auc.train <- e.train@auc
    e.val <- dismo::evaluate(occs.val.pred, bg.val.pred)
    auc.val <- e.val@auc
    # calculate AUC diff as training AUC minus validation AUC with different backgrounds
    auc.diff <- auc.train - auc.val
    # calculate CBI based on the validation background only (do not calculate for jackknife partitions)
    if(partitions != "jackknife") {
      cbi.val <- ecospat::ecospat.boyce(c(bg.val.pred, occs.val.pred), occs.val.pred, PEplot = FALSE)$Spearman.cor
    }else{
      cbi.val <- NA
    }
  }
  
  # get absolute value of AUC diff if requested
  if(other.settings$abs.auc.diff == TRUE) auc.diff <- abs(auc.diff)
  
  ## omission rates
  # get minimum training presence threshold (expected no omission)
  min.train.thr <- min(occs.train.pred)
  or.mtp <- mean(occs.val.pred < min.train.thr)
  # get 10 percentile training presence threshold (expected 0.1 omission)
  pct10.train.thr <- calc.10p.trainThresh(occs.train.pred)
  or.10p <- mean(occs.val.pred < pct10.train.thr)
  
  # perform user-specified validation statistics if available
  if(is.function(user.eval)) {
    user.eval.out <- user.eval(enm, occs.train.z, occs.val.z, bg.train.z, bg.val.z, mod.k, nk, other.settings, partitions)  
  }else{
    user.eval.out <- NULL
  }
  
  # gather all evaluation statistics for k
  out.df <- data.frame(auc.val = auc.val, auc.diff = auc.diff, cbi.val = cbi.val, or.mtp = or.mtp, or.10p = or.10p)
  if(!is.null(user.eval.out)) out.df <- cbind(out.df, user.eval.out)
  
  return(out.df)
}

#' @rdname tune.enm
tune.parallel <- function(d, envs, enm, partitions, tune.tbl, other.settings, user.val.grps, numCores, parallelType, user.eval, quiet) {
  # set up parallel processing functionality
  allCores <- parallel::detectCores()
  if (is.null(numCores)) {
    numCores <- allCores
  }
  cl <- parallel::makeCluster(numCores)
  n <- ifelse(!is.null(tune.tbl), nrow(tune.tbl), 1)
  if(quiet != TRUE) pb <- txtProgressBar(0, n, style = 3)
  if(quiet != TRUE) progress <- function(n) setTxtProgressBar(pb, n)  
  
  if(parallelType == "doParallel") {
    doParallel::registerDoParallel(cl)
    opts <- NULL
  } else if(parallelType == "doSNOW") {
    doSNOW::registerDoSNOW(cl)
    if(quiet != TRUE) opts <- list(progress=progress) else opts <- NULL
  }
  numCoresUsed <- foreach::getDoParWorkers()
  if(quiet != TRUE) message(paste0("\nOf ", allCores, " total cores using ", numCoresUsed, "..."))
  if(quiet != TRUE) message(paste0("Running in parallel using ", parallelType, "..."))
  
  results <- foreach::foreach(i = 1:n, .options.snow = opts, .export = "cv.enm") %dopar% {
    cv.enm(d, envs, enm, partitions, tune.tbl[i,], other.settings, user.val.grps, user.eval, quiet)
  }
  if(quiet != TRUE) close(pb)
  parallel::stopCluster(cl)
  return(results)
}

#' @rdname tune.enm
tune.regular <- function(d, envs, enm, partitions, tune.tbl, other.settings, user.val.grps, updateProgress, user.eval, quiet) {
  results <- list()
  n <- ifelse(!is.null(tune.tbl), nrow(tune.tbl), 1)
  
  # set up the console progress bar
  if(quiet != TRUE) pb <- txtProgressBar(0, n, style = 3)
  
  for(i in 1:n) {
    # and (optionally) the shiny progress bar (updateProgress)
    if(n > 1) {
      if(is.function(updateProgress)) {
        text <- paste0('Running ', paste(as.character(tune.tbl[i,]), collapse = ""), '...')
        updateProgress(detail = text)
      }
      if(quiet != TRUE) setTxtProgressBar(pb, i)
    }
    # set the current tune settings
    tune.i <- tune.tbl[i,]
    results[[i]] <- cv.enm(d, envs, enm, partitions, tune.i, other.settings, user.val.grps, user.eval, quiet)
  }
  if(quiet != TRUE) close(pb)
  return(results)
}

#' @param tune.i vector of single set of tuning parameters
#' @rdname tune.enm
cv.enm <- function(d, envs, enm, partitions, tune.i, other.settings, user.val.grps, user.eval, quiet) {
  envs.names <- names(d[, 3:(ncol(d)-2)])
  # unpack predictor variable values for occs and bg
  occs.xy <- d %>% dplyr::filter(pb == 1) %>% dplyr::select(1:2)
  occs.z <- d %>% dplyr::filter(pb == 1) %>% dplyr::select(all_of(envs.names))
  bg.xy <- d %>% dplyr::filter(pb == 0) %>% dplyr::select(1:2)
  bg.z <- d %>% dplyr::filter(pb == 0) %>% dplyr::select(all_of(envs.names))
  # build the full model from all the data
  mod.full.args <- enm@args(occs.z, bg.z, tune.i, other.settings)
  mod.full <- do.call(enm@fun, mod.full.args)
  if(is.null(mod.full)) stop('Training model is NULL. Consider changing the tuning parameters.')
  # make full model prediction as raster using raster envs (if raster envs exists) 
  # or full model prediction table using the occs and bg values (if raster envs does not exist)
  if(!is.null(envs)) {
    pred.envs <- envs
  }else{
    pred.envs <- d %>% dplyr::select(all_of(envs.names))
  }
  # if using BIOCLIM, put the current tail specification into other.settings
  # for use in doing evaluations and making predictions
  if(enm@name == "bioclim") {
    other.settings$tails <- tune.i
  }
  mod.full.pred <- enm@predict(mod.full, pred.envs, other.settings)
  # get evaluation statistics for training data
  train <- tune.train(enm, occs.z, bg.z, mod.full, mod.full.pred, envs, other.settings, partitions, quiet)
  # make training stats table
  tune.args.col <- paste(names(tune.i), tune.i, collapse = "_", sep = ":")
  train.stats.df <- data.frame(tune.args = tune.args.col, stringsAsFactors = FALSE) %>% cbind(train)
  
  # define number of grp (the value of "k") for occurrences
  # k is 1 for partition "testing"
  # k is 0 for partitions "none" and "user"
  occGrps <- unique(d[d$pb == 1, "grp"])
  if(length(occGrps) == 1 & 0 %in% occGrps) {
    nk <- 0
  }else{
    nk <- length(occGrps)
  }
  
  # if no partitions, return results without cv.stats
  if(nk == 0) {
    cv.res <- list(mod.full = mod.full, mod.full.pred = mod.full.pred, train.stats = train.stats.df, cv.stats = NULL)
    return(cv.res)
  }
  
  # list to contain cv statistics
  cv.stats <- list()
  
  for(k in 1:nk) {
    # assign partitions for training and validation occurrence data and for background data
    occs.train.xy <- d %>% dplyr::filter(pb == 1, grp != k) %>% dplyr::select(1:2)
    occs.train.z <- d %>% dplyr::filter(pb == 1, grp != k) %>% dplyr::select(all_of(envs.names))
    bg.train.z <- d %>% dplyr::filter(pb == 0, grp != k) %>% dplyr::select(all_of(envs.names))
    if(is.null(user.val.grps)) {
      occs.val.xy <- d %>% dplyr::filter(pb == 1, grp == k) %>% dplyr::select(1:2)
      occs.val.z <- d %>% dplyr::filter(pb == 1, grp == k) %>% dplyr::select(all_of(envs.names))
      bg.val.z <- d %>% dplyr::filter(pb == 0, grp == k) %>% dplyr::select(all_of(envs.names))
    }else{
      # assign partitions for training and validation occurrence data and for background data based on user data
      occs.val.xy <- user.val.grps %>% dplyr::filter(grp == k) %>% dplyr::select(1:2)
      occs.val.z <- user.val.grps %>% dplyr::filter(grp == k) %>% dplyr::select(all_of(envs.names))
      bg.val.z <- d %>% dplyr::filter(pb == 0, grp == k) %>% dplyr::select(envs.names)
    }
    
    # if no cross validation (nk = 1), define the model used for evaluation (mod.k) 
    # as the full model (mod.full) to avoid having to refit the same model
    if(nk != 1) {
      # define model arguments for current model k
      mod.k.args <- enm@args(occs.train.z, bg.train.z, tune.i, other.settings)
      # run the current model k
      mod.k <- tryCatch({
        do.call(enm@fun, mod.k.args)  
      }, error = function(cond) {
        if(quiet != TRUE) message(paste0("\n", cond, "\n"))
        # Choose a return value in case of error
        return(NULL)
      })
    }else{
      mod.k <- mod.full
    }
    
    # if model is NULL for some reason, continue but report to user
    if(is.null(mod.k)) {
      if(quiet != TRUE) message(paste0("\nThe model for settings ", paste(names(tune.i), tune.i, collapse = ", "), " for partition ", 
                                       k, " failed (resulted in NULL). Consider changing partitions. Cross validation averages will ignore this model."))
      next
    }
    
    validate <- tune.validate(enm, occs.train.z, occs.val.z, bg.train.z, bg.val.z, mod.k, nk, other.settings, partitions, user.eval, quiet)
    
    # put into list as one-row data frame for easy binding
    cv.stats[[k]] <- data.frame(tune.args = tune.args.col, fold = k, stringsAsFactors = FALSE) %>% cbind(validate)
  } 
  
  cv.stats.df <- dplyr::bind_rows(cv.stats)
  
  cv.res <- list(mod.full = mod.full, mod.full.pred = mod.full.pred, 
                 train.stats = train.stats.df, cv.stats = cv.stats.df)
  
  return(cv.res)
}
