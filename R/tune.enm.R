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
tune.parallel <- function(d, envs, enm, partitions, tune.tbl, other.settings, user.val.grps, numCores, parallelType, quiet) {
  # set up parallel processing functionality
  allCores <- parallel::detectCores()
  if (is.null(numCores)) {
    numCores <- allCores
  }
  cl <- parallel::makeCluster(numCores)
  n <- ifelse(nrow(tune.tbl) > 0, nrow(tune.tbl), 1)
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
  if(quiet != TRUE) message(paste0("Of ", allCores, " total cores using ", numCoresUsed, "..."))
  if(quiet != TRUE) message(paste0("Running in parallel using ", parallelType, "..."))
  
  results <- foreach::foreach(i = 1:n, .packages = enm.pkgs(enm), .options.snow = opts, .export = "cv.enm") %dopar% {
    cv.enm(d, envs, enm, partitions, tune.tbl[i,], other.settings, user.val.grps, quiet)
  }
  if(quiet != TRUE) close(pb)
  parallel::stopCluster(cl)
  return(results)
}

#' @rdname tune.enm
tune.regular <- function(d, envs, enm, partitions, tune.tbl, other.settings, user.val.grps, updateProgress, quiet) {
  results <- list()
  n <- ifelse(nrow(tune.tbl) > 0, nrow(tune.tbl), 1)
  
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
    results[[i]] <- cv.enm(d, envs, enm, partitions, tune.i, other.settings, user.val.grps, quiet)
  }
  if(quiet != TRUE) close(pb)
  return(results)
}

#' @param tune.i vector of single set of tuning parameters

#' @rdname tune.enm
cv.enm <- function(d, envs, enm, partitions, tune.i, other.settings, user.val.grps, quiet) {
  envs.names <- names(d[, 3:(ncol(d)-2)])
  # unpack predictor variable values for occs and bg
  occs.xy <- d %>% dplyr::filter(pb == 1) %>% dplyr::select(1:2)
  occs.z <- d %>% dplyr::filter(pb == 1) %>% dplyr::select(all_of(envs.names))
  bg.xy <- d %>% dplyr::filter(pb == 0) %>% dplyr::select(1:2)
  bg.z <- d %>% dplyr::filter(pb == 0) %>% dplyr::select(all_of(envs.names))
  # build the full model from all the data
  mod.full.args <- enm@args(occs.z, bg.z, tune.i, other.settings)
  mod.full <- do.call(enm@fun, mod.full.args)
  if(is.null(mod.full)) stop("Training model is NULL. Consider changing the tuning parameters.")
  # make full model prediction as raster using raster envs (if raster envs exists) 
  # or full model prediction table using the occs and bg values (if raster envs does not exist)
  if(!is.null(envs)) {
    pred.envs <- envs
  }else{
    pred.envs <- d %>% dplyr::select(all_of(envs.names))
  }
  mod.full.pred <- enm@pred(mod.full, pred.envs, other.settings)
  # get evaluation statistics for training data
  eval.train <- enm@eval.train(occs.xy, bg.xy, occs.z, bg.z, mod.full, mod.full.pred, envs, other.settings)
  # make training stats table
  tune.args.col <- paste(tune.i, collapse = "_")
  train.stats.df <- data.frame(tune.args = tune.args.col, stringsAsFactors = FALSE) %>% cbind(eval.train)
  
  # define number of grp (the value of "k") for occurrences
  # k is 1 for partition "independent"
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
      # bg.test.z <- d %>% dplyr::filter(pb == 0, grp == k) %>% dplyr::select(envs.names)  
    }else{
      # assign partitions for training and validation occurrence data and for background data based on user data
      occs.val.xy <- user.val.grps %>% dplyr::filter(grp == k) %>% dplyr::select(1:2)
      occs.val.z <- user.val.grps %>% dplyr::filter(grp == k) %>% dplyr::select(all_of(envs.names))
      # bg.test.z <- d %>% dplyr::filter(pb == 0, grp == k) %>% dplyr::select(envs.names)  
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
    
    eval.validate <- enm@eval.validate(occs.val.xy, occs.train.xy, bg.xy, occs.train.z, occs.val.z, 
                               bg.z, mod.k, nk, envs, other.settings)
    
    # put into list as one-row data frame for easy binding
    cv.stats[[k]] <- data.frame(tune.args = tune.args.col, fold = k, stringsAsFactors = FALSE) %>% cbind(eval.validate)
  } 
  
  cv.stats.df <- dplyr::bind_rows(cv.stats)
  
  cv.res <- list(mod.full = mod.full, mod.full.pred = mod.full.pred, 
                 train.stats = train.stats.df, cv.stats = cv.stats.df)
  
  return(cv.res)
}
