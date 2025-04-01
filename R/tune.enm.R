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
tune.train <- function(pred.occs, pred.bg) {
  # training AUC
  e <- predicts::pa_evaluate(pred.occs, pred.bg)
  auc.train <- e@stats$auc
  
  # calculate training CBI with background, not raster -- this is in order
  # to calculate CBI in a standardized way across training and validation data for 
  # spatial partitions (ENMeval does not mask rasters to partition areas and hence 
  # does not have this information)
  # combine occs and bg predictions as input for "fit" because the interval
  # is determined from "fit" only, and if occs have higher predictions than
  # bg, the interval will be cut short
  if(other.settings$ecospat.use == TRUE) {
    cbi.train <- ecospat::ecospat.boyce(c(pred.bg, pred.occs), pred.occs, PEplot = FALSE)$cor
  }else{
    cbi.train <- NA
  }
  
  out.df <- data.frame(auc.train = auc.train, cbi.train = cbi.train)
  return(out.df)
}

#' @rdname tune.enm
tune.validate <- function(pred.occs.train, pred.occs.val, pred.bg.train, 
                          pred.bg.val, other.settings, user.eval) {
  # if validation.bg == "full", calculate training and validation AUC and CBI based on the full background 
  # (training + validation background for all statistics) 
  # see Radosavljevic & Anderson 2014 for an example of calculating validation AUC with spatial partitions over a shared background
  if(other.settings$validation.bg == "full") {
    e.train <- predicts::pa_evaluate(pred.occs.train, c(pred.bg.train, pred.bg.val))
    auc.train <- e.train@stats$auc
    e.val <- predicts::pa_evaluate(pred.occs.val, c(pred.bg.train, pred.bg.val))
    auc.val <- e.val@stats$auc
    # calculate AUC diff as training AUC minus validation AUC with a shared background
    auc.diff <- auc.train - auc.val
    # calculate CBI based on the full background (do not calculate for jackknife partitions)
    if(other.settings$ecospat.use == TRUE & other.settings$cbi.cv == TRUE) {
      cbi.val <- ecospat::ecospat.boyce(c(pred.bg.train, pred.bg.val, pred.occs.val), pred.occs.val, PEplot = FALSE)$cor  
    }else{
      cbi.val <- NA
    }
    # if validation.bg == "partition", calculate training and validation AUC and CBI based on the partitioned backgrounds only 
    # (training background for training statistics and validation background for validation statistics) 
  }else if(other.settings$validation.bg == "partition") {
    e.train <- predicts::pa_evaluate(pred.occs.train, pred.bg.train)
    auc.train <- e.train@stats$auc
    e.val <- predicts::pa_evaluate(pred.occs.val, pred.bg.val)
    auc.val <- e.val@stats$auc
    # calculate AUC diff as training AUC minus validation AUC with different backgrounds
    auc.diff <- auc.train - auc.val
    # calculate CBI based on the validation background only (do not calculate for jackknife partitions)
    if(other.settings$ecospat.use == TRUE & other.settings$cbi.cv == TRUE) {
      cbi.val <- ecospat::ecospat.boyce(c(pred.bg.val, pred.occs.val), pred.occs.val, PEplot = FALSE)$cor
    }else{
      cbi.val <- NA
    }
  }
  
  # get absolute value of AUC diff if requested
  if(other.settings$abs.auc.diff == TRUE) auc.diff <- abs(auc.diff)
  
  ## omission rates
  # get minimum training presence threshold (expected no omission)
  min.train.thr <- min(pred.occs.train)
  or.mtp <- mean(pred.occs.val < min.train.thr)
  # get 10 percentile training presence threshold (expected 0.1 omission)
  pct10.train.thr <- quantile(pred.occs.train, 0.1) |> as.numeric()
  or.10p <- mean(pred.occs.val < pct10.train.thr)
  
  # perform user-specified validation statistics if available
  if(is.function(user.eval)) {
    vars <- list(pred.occs.train, pred.occs.val, pred.bg.train, 
                 pred.bg.val, other.settings)
    names(vars) <- c("pred.occs.train", "pred.occs.val", "pred.bg.train", 
                     "pred.bg.val", "other.settings")
    user.eval.out <- user.eval(vars)
  }else{
    user.eval.out <- NULL
  }
  
  # gather all evaluation statistics for k
  out.df <- data.frame(auc.val = auc.val, auc.diff = auc.diff, 
                       cbi.val = cbi.val, or.mtp = or.mtp, or.10p = or.10p)
  # if user.eval was used, add it to table
  if(!is.null(user.eval.out)) out.df <- cbind(out.df, user.eval.out)
  
  return(out.df)
}

#' @rdname tune.enm
tuneParallel <- function(occs.z, bg.z, grps, enm, partitions, tune.tbl, doClamp, 
                         other.settings, partition.settings, user.val.grps, 
                         occs.testing.z, numCores, user.eval, algorithm, updateProgress, quiet) {
  # set up parallel processing functionality
  allCores <- parallel::detectCores()
  if (is.null(numCores)) {
    numCores <- allCores
  }
  cl <- parallel::makeCluster(numCores, setup_strategy = "sequential")
  parallel::clusterExport(cl, c("d", "enm", "partitions", "tune.tbl", 
                                "doClamp", "other.settings", 
                                "partition.settings", "user.val.grps", 
                                "occs.testing.z", "user.eval", "algorithm", 
                                "quiet", "cv.enm", "tune.train", "clamp.vars", 
                                "tune.validate"),
                          envir = environment())
  n <- ifelse(!is.null(tune.tbl), nrow(tune.tbl), 1)
  
  msg(paste0("\nOf ", allCores, " total cores using ", numCores, "..."))
  msg("Running in parallel ...")
  
  par.cv.enm <- function(i) {
    cv.enm(occs.z, bg.z, grps, enm, partitions, tune.tbl[i,], doClamp, 
           other.settings, partition.settings, user.val.grps, occs.testing.z, 
           user.eval, algorithm, quiet)
  }
  
  results <- parallel::parLapply(cl, 1:n, par.cv.enm)
  return(results)
}

#' @rdname tune.enm
tune <- function(occs.z, bg.z, grps, enm, partitions, tune.tbl, doClamp, 
                 other.settings, partition.settings, user.val.grps, 
                 occs.testing.z, numCores, user.eval, algorithm, updateProgress, quiet) {
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
    results[[i]] <- cv.enm(occs.z, bg.z, grps, enm, partitions, tune.tbl.i, 
                           doClamp,other.settings, partition.settings, 
                           user.val.grps, occs.testing.z, user.eval, algorithm, quiet)
  }
  if(quiet != TRUE) close(pb)
  return(results)
}


#' @rdname tune.enm
cv.enm <- function(occs.z, bg.z, grps, enm, partitions, tune.tbl.i, doClamp, 
                   other.settings, partition.settings, user.val.grps, 
                   occs.testing.z, user.eval, algorithm, quiet) {
  envs.names <- names(occs.z)
  # unpack predictor variable values for occs and bg
  # occs.xy <- d |> dplyr::filter(pb == 1) |> dplyr::select(1:2)
  # occs.z <- d |> dplyr::filter(pb == 1) |> dplyr::select(dplyr::all_of(envs.names))
  # bg.xy <- d |> dplyr::filter(pb == 0) |> dplyr::select(1:2)
  # bg.z <- d |> dplyr::filter(pb == 0) |> dplyr::select(dplyr::all_of(envs.names))
  
  # define number of grp (the value of "k") for occurrences
  nk <- length(unique(grps$occs.grp))
  
  ## build the full model from all the data
  # assign arguments
  mod.full.args <- enm@args(occs.z, bg.z, tune.tbl.i, other.settings)
  # run training model with specified arguments
  mod.full <- do.call(enm@fun, mod.full.args)
  
  if(is.null(mod.full)) stop('Training model is NULL. Consider changing the tuning parameters or inputting more background points.')
  
  # as BIOCLIM can be tuned with different "tails" settings that affect not the 
  # model but the prediction, these settings need to be moved from tune.tbl.i
  # to other.settings as eval.predict() does not accept tune.args
  if(algorithm == "bioclim") other.settings$tails <- tune.tbl.i$tails
  # get model predictions for occs and bg
  pred.occs <- enm@predict(mod.full, occs.z, other.settings)
  pred.bg <- enm@predict(mod.full, bg.z, other.settings)
  pred.full <- c(pred.occs, pred.bg)
  # get evaluation statistics for training data
  train <- tune.train(pred.occs, pred.bg)
  # make training stats table
  tune.args.col <- paste(names(tune.tbl.i), tune.tbl.i, collapse = "_", sep = ".")
  train.stats.df <- data.frame(tune.args = tune.args.col, stringsAsFactors = FALSE) |> cbind(train)
  
  # if no partitions, return results without cv.stats
  if(partitions == "none") {
    cv.res <- list(mod.full = mod.full, pred.full = pred.full, 
                   train.stats = train.stats.df, cv.stats = NULL)
    return(cv.res)
  }
  
  # calculate validation stats on user-specified testing dataset
  if(partitions == "testing") {
    if(doClamp == TRUE) {
      occs.testing.z <- clamp.vars(orig.vals = occs.testing.z, 
                                   ref.vals = rbind(occs.z, bg.z), 
                                   left = other.settings$clamp.directions$left, 
                                   right = other.settings$clamp.directions$right, 
                                   categoricals = other.settings$categoricals)
    }
    # make model prediction of occs.training
    pred.occs.testing <- enm@predict(mod.full, occs.testing.z, other.settings)
    # calculate validation stats
    validate <- tune.validate(pred.occs, pred.occs.testing, pred.bg, 
                              bg.val.z = data.frame(), other.settings, 
                              user.eval)
    test.stats.df <- data.frame(tune.args = tune.args.col, fold = 0, 
                                stringsAsFactors = FALSE) |> cbind(validate)
    cv.res <- list(mod.full = mod.full, pred.full = pred.full, 
                   train.stats = train.stats.df, cv.stats = test.stats.df) 
    return(cv.res)
  }
  
  # list to contain cv statistics
  cv.stats <- list()
  
  for(k in 1:nk) {
    # if doClamp is on, make sure that the validation data for each validation model is also clamped
    # this means for each partition, making sure no values in validation data are more extreme than those in training data
    occs.train.z <- occs.z[-which(grps$occs.grp == k),]
    occs.val.z <- occs.z[which(grps$occs.grp == k),]
    bg.train.z <- bg.z[-which(grps$bg.grp == k),]
    bg.val.z <- bg.z[which(grps$bg.grp == k),]
    # if bg not partitioned for cross-validation, this will be an empty df
    # if so, make bg.val the whole bg
    if(nrow(bg.val.z) == 0) {
      bg.val.z <- bg.z
    }
    
    # if clamping is on, clamp the validation data by the training data
    if(doClamp == TRUE) {
      occs.val.z <- clamp.vars(orig.vals = occs.val.z, 
                               ref.vals = rbind(occs.train.z, bg.train.z), 
                               left = other.settings$clamp.directions$left, 
                               right = other.settings$clamp.directions$right, 
                               categoricals = other.settings$categoricals)
      bg.val.z <- clamp.vars(orig.vals = bg.val.z, 
                               ref.vals = rbind(occs.train.z, bg.train.z), 
                               left = other.settings$clamp.directions$left, 
                               right = other.settings$clamp.directions$right, 
                               categoricals = other.settings$categoricals)
    }
    
    ## NOT SURE ABOUT HERE
    if(!is.null(user.val.grps)) {
      # assign partitions for training and validation occurrence data and for background data based on user data
      occs.val.z <- user.val.grps |> dplyr::filter(grp == k) |> dplyr::select(dplyr::all_of(envs.names))
      bg.val.z <- d |> dplyr::filter(pb == 0, grp == k) |> dplyr::select(envs.names)
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
    
    # get model predictions for training and validation data for partition k
    pred.occs.train <- enm@predict(mod.k, occs.train.z, other.settings)
    pred.occs.val <- enm@predict(mod.k, occs.val.z, other.settings)
    pred.bg.train <- enm@predict(mod.k, bg.train.z, other.settings)
    # if the background was partitioned, get predictions for this subset; if not, it will be NULL
    if(nrow(bg.val.z) > 0) {
      pred.bg.val <- enm@predict(mod.k, bg.val.z, other.settings)
    }else{
      bg.val.pred <- NULL
    }
    
    # if model ran without problems, make validation table
    if(!is.null(mod.k)) {
      validate <- tune.validate(pred.occs.train, pred.occs.val, pred.bg.train, 
                                pred.bg.val, other.settings, user.eval)  
    # if model is NULL for some reason, continue but report to user
    }else{
      mod.name <- paste(names(tune.tbl.i), tune.tbl.i, collapse = ", ")
      msg(paste0("The results were NULL for model with settings ", mod.name,
                                       " for partition ", k, 
                                       ". These settings will have NA results."))
      validate <- data.frame(auc.val = NA, auc.diff = NA, cbi.val = NA, 
                             or.mtp = NA, or.10p = NA)
    }
    
    # put into list as one-row data frame for easy binding
    cv.stats[[k]] <- data.frame(tune.args = tune.args.col, fold = k, 
                                stringsAsFactors = FALSE) |> cbind(validate)
  } 
  
  cv.stats.df <- dplyr::bind_rows(cv.stats)
  
  cv.res <- list(mod.full = mod.full, pred.full = pred.full, 
                 train.stats = train.stats.df, cv.stats = cv.stats.df)
  
  return(cv.res)
}
