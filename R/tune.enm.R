#' @title Iterate tuning of ENMs
#' @description Internal functions to tune and summarize results for ecological niche models (ENMs) iteratively across a range of user-specified tuning settings. 
#' See \link{ENMevaluate} for descriptions of shared arguments.
#' Function \code{tune.parallel()} tunes ENMs with parallelization. Function \code{cv.enm()} calculates training and validation evaluation statistics for one set of specified tuning parameters.
#' @aliases tune.parallel tune.regular cv.enm
#' @param enm \link{ENMdetails} object
#' @param occs.z data.frame: the envs values for the coordinates at the full dataset occurrence records
#' @param bg.z data.frame: the envs values for the coordinates at the full dataset background records
#' @param mod.full model object: the model trained on the full dataset
#' @param tune.tbl data frame: all combinations of tuning parameters
#' @param tune.tbl.i vector: single set of tuning parameters
#' @param partitions character: name of partitioning technique (see \code{?partitions})
#' @param occs.train.z data.frame: the envs values for the coordinates at the training occurrence records
#' @param occs.val.z data.frame: the envs values for the coordinates at the validation occurrence records
#' @param occs.testing.z data.frame: when fully withheld testing data is provided, the envs values for the 
#' coordinates at the testing occurrence records
#' @param bg.train.z data.frame: the envs values for the coordinates at the training background records
#' @param bg.val.z data.frame: the envs values for the coordinates at the validation background records
#' @param mod.k model object: the model trained on the training dataset that becomes evaluated on the validation data
#' @param nk numeric: the number of folds (i.e., partitions) -- will be equal to \code{kfolds} for random partitions
#' @param d data frame: data frame from \code{ENMevaluate()} with occurrence and background coordinates (or coordinates plus predictor variable values) and partition group values
#' @description Validation CBI is calculated here with background values, not raster data, in order
#' to standardize the methodology for both training and validation data for spatial partitions, as ENMeval
#' does not mask rasters to partition areas and hence does not have partitioned raster data. Further, 
#' predictions for occurrence and background localities are combined as input for the parameter "fit" in 
#' \code{ecospat::ecospat_boyce()} because the interval is determined from "fit" only, and if test occurrences 
#' all have higher predictions than the background, the interval will be cut short.
#' @inheritParams ENMevaluate
#' @name tune.enm
NULL

#' @rdname tune.enm
tune.train <- function(enm, occs.z, bg.z, mod.full, tune.tbl.i, other.settings, partitions, quiet) {
  # get model predictions for training data
  occs.pred <- enm@predict(mod.full, occs.z, other.settings)
  bg.pred <- enm@predict(mod.full, bg.z, other.settings)
  
  # training AUC
  e <- predicts::pa_evaluate(occs.pred, bg.pred)
  auc.train <- e@stats$auc
  
  # calculate training CBI with background, not raster -- this is in order
  # to calculate CBI in a standardized way across training and validation data for 
  # spatial partitions (ENMeval does not mask rasters to partition areas and hence 
  # does not have this information)
  # combine occs and bg predictions as input for "fit" because the interval
  # is determined from "fit" only, and if occs have higher predictions than
  # bg, the interval will be cut short
  if(other.settings$ecospat.use == TRUE) {
    cbi.train <- ecospat::ecospat.boyce(c(bg.pred, occs.pred), occs.pred, PEplot = FALSE)$cor
  }else{
    cbi.train <- NA
  }
  
  out.df <- data.frame(auc.train = auc.train, cbi.train = cbi.train)
  return(out.df)
}

#' @rdname tune.enm
tune.validate <- function(enm, occs.train.z, occs.val.z, bg.train.z, bg.val.z, 
                          mod.k, nk, tune.tbl.i, other.settings, partitions, 
                          user.eval, quiet) {
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
    e.train <- predicts::pa_evaluate(occs.train.pred, c(bg.train.pred, bg.val.pred))
    auc.train <- e.train@stats$auc
    e.val <- predicts::pa_evaluate(occs.val.pred, c(bg.train.pred, bg.val.pred))
    auc.val <- e.val@stats$auc
    # calculate AUC diff as training AUC minus validation AUC with a shared background
    auc.diff <- auc.train - auc.val
    # calculate CBI based on the full background (do not calculate for jackknife partitions)
    if(other.settings$ecospat.use == TRUE) {
      if(partitions != "jackknife") {
        cbi.val <- ecospat::ecospat.boyce(c(bg.train.pred, bg.val.pred, occs.val.pred), occs.val.pred, PEplot = FALSE)$cor  
      }else{
        cbi.val <- NA
      }
    }else{
      cbi.val <- NA
    }
    # if validation.bg == "partition", calculate training and validation AUC and CBI based on the partitioned backgrounds only 
    # (training background for training statistics and validation background for validation statistics) 
  }else if(other.settings$validation.bg == "partition") {
    e.train <- predicts::pa_evaluate(occs.train.pred, bg.train.pred)
    auc.train <- e.train@stats$auc
    e.val <- predicts::pa_evaluate(occs.val.pred, bg.val.pred)
    auc.val <- e.val@stats$auc
    # calculate AUC diff as training AUC minus validation AUC with different backgrounds
    auc.diff <- auc.train - auc.val
    # calculate CBI based on the validation background only (do not calculate for jackknife partitions)
    if(other.settings$ecospat.use == TRUE) {
      if(partitions != "jackknife") {
        cbi.val <- ecospat::ecospat.boyce(c(bg.val.pred, occs.val.pred), occs.val.pred, PEplot = FALSE)$cor
      }else{
        cbi.val <- NA
      }
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
    vars <- list(enm, occs.train.z, occs.val.z, bg.train.z, bg.val.z, mod.k, nk, 
                 other.settings, partitions, occs.train.pred, occs.val.pred,
                 bg.train.pred, bg.val.pred)
    names(vars) <- c("enm", "occs.train.z", "occs.val.z", "bg.train.z", "bg.val.z", "mod.k", "nk", 
                     "other.settings", "partitions", "occs.train.pred", "occs.val.pred",
                     "bg.train.pred", "bg.val.pred")
    user.eval.out <- user.eval(vars)
  }else{
    user.eval.out <- NULL
  }
  
  # gather all evaluation statistics for k
  out.df <- data.frame(auc.val = auc.val, auc.diff = auc.diff, cbi.val = cbi.val, or.mtp = or.mtp, or.10p = or.10p)
  if(!is.null(user.eval.out)) out.df <- cbind(out.df, user.eval.out)
  
  return(out.df)
}

#' @rdname tune.enm
tune <- function(d, enm, partitions, tune.tbl, doClamp, other.settings, 
                 partition.settings, user.val.grps, occs.testing.z, 
                 numCores, parallel, parallelType, user.eval, algorithm, 
                 updateProgress, quiet) {
  if(parallel == TRUE) {
    # set up parallel processing functionality
    allCores <- parallel::detectCores()
    if (is.null(numCores)) {
      numCores <- allCores
    }
    cl <- parallel::makeCluster(numCores, setup_strategy = "sequential")
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
      cv.enm(d, enm, partitions, tune.tbl[i,], doClamp, other.settings, 
             partition.settings, user.val.grps, occs.testing.z, user.eval, 
             algorithm, quiet)
    }
    parallel::stopCluster(cl)
  }else{
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
      tune.tbl.i <- tune.tbl[i,]
      results[[i]] <- cv.enm(d, enm, partitions, tune.tbl.i, doClamp,
                             other.settings, partition.settings, user.val.grps, 
                             occs.testing.z, user.eval, algorithm, quiet)
    }
  }
  if(quiet != TRUE) close(pb)
  return(results)
}

#' @rdname tune.enm
cv.enm <- function(d, enm, partitions, tune.tbl.i, doClamp, other.settings, 
                   partition.settings, user.val.grps, occs.testing.z, 
                   user.eval, algorithm, quiet) {
  envs.names <- names(d[, 3:(ncol(d)-2)])
  # unpack predictor variable values for occs and bg
  occs.xy <- d |> dplyr::filter(pb == 1) |> dplyr::select(1:2)
  occs.z <- d |> dplyr::filter(pb == 1) |> dplyr::select(dplyr::all_of(envs.names))
  bg.xy <- d |> dplyr::filter(pb == 0) |> dplyr::select(1:2)
  bg.z <- d |> dplyr::filter(pb == 0) |> dplyr::select(dplyr::all_of(envs.names))
  
  # define number of grp (the value of "k") for occurrences
  nk <- length(unique(d[d$pb == 1, "grp"]))
  
  # build the full model from all the data
  # assign arguments
  mod.full.args <- enm@args(occs.z, bg.z, tune.tbl.i, other.settings)
  # run training model with specified arguments
  mod.full <- do.call(enm@fun, mod.full.args)
  
  if(is.null(mod.full)) stop('Training model is NULL. Consider changing the tuning parameters or inputting more background points.')
  
  # as bioclim can be tuned with different "tails" settings that affect not the 
  # model but the prediction, these settings need to be moved from tune.tbl.i
  # to other.settings as eval.predict() does not accept tune.args
  if(algorithm == "bioclim") other.settings$tails <- tune.tbl.i$tails
  
  mod.full.pred <- enm@predict(mod.full, rbind(occs.z, bg.z), other.settings)
  # get evaluation statistics for training data
  train <- tune.train(enm, occs.z, bg.z, mod.full, tune.tbl.i, other.settings, partitions, quiet)
  # make training stats table
  tune.args.col <- paste(names(tune.tbl.i), tune.tbl.i, collapse = "_", sep = ".")
  train.stats.df <- data.frame(tune.args = tune.args.col, stringsAsFactors = FALSE) |> cbind(train)
  
  # if no partitions, return results without cv.stats
  if(partitions == "none") {
    cv.res <- list(mod.full = mod.full, mod.full.pred = mod.full.pred, train.stats = train.stats.df, cv.stats = NULL)
    return(cv.res)
  }
  
  if(partitions == "testing") {
    bg.val.z <- data.frame()
    occs.testing.zEnvs <- occs.testing.z |> dplyr::select(dplyr::all_of(envs.names))
    if(doClamp == TRUE) {
      occs.testing.zEnvs <- clamp.vars(orig.vals = occs.testing.zEnvs, ref.vals = rbind(occs.z, bg.z), 
                                       left = other.settings$clamp.directions$left, right = other.settings$clamp.directions$right, 
                                       categoricals = other.settings$categoricals)
    }
    validate <- tune.validate(enm, occs.z, occs.testing.zEnvs, bg.z, bg.val.z, mod.full, 0, tune.tbl.i, other.settings, partitions, user.eval, quiet)
    test.stats.df <- data.frame(tune.args = tune.args.col, fold = 0, stringsAsFactors = FALSE) |> cbind(validate)
    cv.res <- list(mod.full = mod.full, mod.full.pred = mod.full.pred, train.stats = train.stats.df, cv.stats = test.stats.df) 
    return(cv.res)
  }
  
  # list to contain cv statistics
  cv.stats <- list()
  
  for(k in 1:nk) {
    # assign partitions for training and validation occurrence data and for background data
    occs.train.z <- d |> dplyr::filter(pb == 1, grp != k) |> dplyr::select(dplyr::all_of(envs.names))
    bg.train.z <- d |> dplyr::filter(pb == 0, grp != k) |> dplyr::select(dplyr::all_of(envs.names))
    if(is.null(user.val.grps)) {
      occs.val.z <- d |> dplyr::filter(pb == 1, grp == k) |> dplyr::select(dplyr::all_of(envs.names))
      bg.val.z <- d |> dplyr::filter(pb == 0, grp == k) |> dplyr::select(dplyr::all_of(envs.names))
    }else{
      # assign partitions for training and validation occurrence data and for background data based on user data
      occs.val.z <- user.val.grps |> dplyr::filter(grp == k) |> dplyr::select(dplyr::all_of(envs.names))
      bg.val.z <- d |> dplyr::filter(pb == 0, grp == k) |> dplyr::select(envs.names)
    }
    
    # if doClamp is on, make sure that the validation data for each validation model is also clamped
    # this means for each partition, making sure no values in validation data are more extreme than those in training data
    if(doClamp == TRUE) {
      val.z <- clamp.vars(orig.vals = rbind(occs.val.z, bg.val.z), 
                          ref.vals = rbind(occs.train.z, bg.train.z), 
                          left = other.settings$clamp.directions$left, 
                          right = other.settings$clamp.directions$right, 
                          categoricals = other.settings$categoricals)
      occs.val.z <- val.z[1:nrow(occs.val.z),]
      if(nrow(bg.val.z) > 0) bg.val.z <- val.z[(nrow(occs.val.z)+1):nrow(bg.val.z),]  
    }
    
    # define model arguments for current model k
    mod.k.args <- enm@args(occs.train.z, bg.train.z, tune.tbl.i, other.settings)
    # run model k with specified arguments
    mod.k <- tryCatch({
      do.call(enm@fun, mod.k.args)  
    }, error = function(cond) {
      if(quiet != TRUE) message(paste0("\n", cond, "\n"))
      # Choose a return value in case of error
      return(NULL)
    })  
    
    # if model is NULL for some reason, continue but report to user
    if(is.null(mod.k)) {
      if(quiet != TRUE) message(paste0("\nThe results were NULL for model with settings ", 
                                       paste(names(tune.tbl.i), 
                                             tune.tbl.i, 
                                             collapse = ", "), 
                                       " for partition ", k, 
                                       ". These settings will have NA results."))
      validate <- data.frame(auc.val = NA, auc.diff = NA, cbi.val = NA, 
                             or.mtp = NA, or.10p = NA)
    }else{
      validate <- tune.validate(enm, occs.train.z, occs.val.z, bg.train.z, 
                                bg.val.z, mod.k, nk, tune.tbl.i, other.settings, 
                                partitions, user.eval, quiet)  
    }
    
    # put into list as one-row data frame for easy binding
    cv.stats[[k]] <- data.frame(tune.args = tune.args.col, fold = k, 
                                stringsAsFactors = FALSE) |> cbind(validate)
  } 
  
  cv.stats.df <- dplyr::bind_rows(cv.stats)
  
  cv.res <- list(mod.full = mod.full, mod.full.pred = mod.full.pred, 
                 train.stats = train.stats.df, cv.stats = cv.stats.df)
  
  return(cv.res)
}
