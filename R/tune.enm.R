#' @title Functions for tuning ENMs
#' @description Internal functions to tune and summarize results for ecological niche models (ENMs) iteratively across a range of user-specified tuning settings.
#' @aliases tune.parallel tune.regular cv.enm
#' @param occs.vals matrix or data frame of environmental values corresponding
#' to occurrence localities, intended to be input when environmental rasters
#' are not used (\code{envs} is NULL) 
#' @param bg.vals matrix or data frame of environmental values corresponding
#' to background (or pseudo-absence) localities, intended to be input when 
#' environmental rasters are not used (\code{envs} is NULL) 
#' @param occs.grp numeric vector of partition group (fold) for each
#' occurrence locality, intended for user-defined partitions
#' @param bg.grp numeric vector of partition group (fold) for each background 
#' (or pseudo-absence) locality, intended for user-defined partitions
#' @param envs Raster* object of environmental variables (must be in 
#' same geographic projection as occurrence data)
#' @param enm Object of class \link{ENMdetails}.
#' @param partitions character of name of partitioning technique (see
#' \code{?partitions})
#' @param tune.tbl Data frame of tuning parameter combinations.
#' @param tune.settings Vector of tune settings from `tune.tbl`.
#' @param other.args named list of any additional model arguments not specified 
#' for tuning
#' @param categoricals character vector of names of categorical 
#' environmental variables
#' @param occs.ind matrix or data frame with two columns for longitude and latitude 
#' of occurrence localities, in that order, intended for independent evaluation;
#' when \code{partitions = "independent"}; these occurrences will be used only 
#' for evaluation, and not for model training, and thus no cross validation will 
#' be done
#' @param doClamp boolean (TRUE or FALSE); if TRUE, clamp model responses; only
#' applicable for Maxent models
#' @param skipRasters boolean (TRUE or FALSE); if TRUE, skip raster predictions
#' @param abs.auc.diff boolean (TRUE or FALSE); if TRUE, take absolute value of
#' AUCdiff; default is TRUE
#' @param numCores boolean (TRUE or FALSE); if TRUE, use specifed number of cores
#' for parallel processing

#' @name tune.enm
NULL

#' @rdname tune.enm
tune.parallel <- function(d, envs, envs.names, enm, partitions, tune.tbl, settings, numCores, parallelType) {
  # set up parallel processing functionality
  allCores <- parallel::detectCores()
  if (is.null(numCores)) {
    numCores <- allCores
  }
  cl <- parallel::makeCluster(numCores)
  n <- ifelse(nrow(tune.tbl) > 0, nrow(tune.tbl), 1)
  pb <- txtProgressBar(0, n, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  
  if(parallelType == "doParallel") {
    doParallel::registerDoParallel(cl)
    opts <- NULL
  } else if(parallelType == "doSNOW") {
    doSNOW::registerDoSNOW(cl)
    opts <- list(progress=progress)
  }
  numCoresUsed <- foreach::getDoParWorkers()
  message(paste0("\nOf ", allCores, " total cores using ", numCoresUsed, "..."))
  
  message(paste0("Running in parallel using ", parallelType, "..."))
  
  results <- foreach::foreach(i = 1:n, .packages = enm.pkgs(enm), .options.snow = opts, .export = "cv.enm") %dopar% {
    cv.enm(d, envs, envs.names, enm, tune.tbl[i,], partitions, settings)
  }
  close(pb)
  parallel::stopCluster(cl)
  return(results)
}

#' @rdname tune.enm
tune.regular <- function(d, envs, envs.names, enm, partitions, tune.tbl, settings, updateProgress) {
  results <- list()
  n <- ifelse(nrow(tune.tbl) > 0, nrow(tune.tbl), 1)
  
  # set up the console progress bar
  pb <- txtProgressBar(0, n, style = 3)
  
  for(i in 1:n) {
    # and (optionally) the shiny progress bar (updateProgress)
    if(n > 1) {
      if(is.function(updateProgress)) {
        text <- paste0('Running ', paste(as.character(tune.tbl[i,]), collapse = ""), '...')
        updateProgress(detail = text)
      }
      setTxtProgressBar(pb, i)
    }
    # set the current tune settings
    tune.i <- tune.tbl[i,]
    results[[i]] <- cv.enm(d, envs, envs.names, enm, tune.i, partitions, settings)
  }
  close(pb)
  return(results)
}

#' @rdname tune.enm
cv.enm <- function(d, envs, envs.names, enm, tune.i, partitions, settings) {
  # unpack predictor variable values for occs and bg
  occs.vals <- d %>% dplyr::filter(pb == 1) %>% dplyr::select(envs.names)
  bg.vals <- d %>% dplyr::filter(pb == 0) %>% dplyr::select(envs.names)
  # build the full model from all the data
  mod.full.args <- enm@args(occs.vals, bg.vals, tune.i, settings$other.args)
  mod.full <- do.call(enm@fun, mod.full.args)
  # calculate training auc
  e.train <- enm@eval(occs.vals, bg.vals, mod.full, settings$other.args, settings$doClamp)
  auc.train <- e.train@auc
  tune.args.col <- paste(tune.i, collapse = "_")
  train.stats.df <- data.frame(tune.args = tune.args.col, auc.train = auc.train, stringsAsFactors = FALSE)
  # if rasters selected and envs is not NULL, predict raster for the full model
  if(settings$skipRasters == FALSE & !is.null(envs)) {
    mod.full.pred <- enm@pred(mod.full, envs, settings$other.args, settings$doClamp, settings$pred.type)
    # training CBI
    d.occs.xy <- d %>% dplyr::filter(pb == 1) %>% dplyr::select(1:2)
    cbi.train <- ecospat::ecospat.boyce(mod.full.pred, d.occs.xy, PEplot = FALSE)
    train.stats.df$cbi.train <- cbi.train$Spearman.cor
  }else{
    mod.full.pred <- raster::stack()
  }
  
  # define number of grp (the value of "k") for occurrences
  nk <- length(unique(d[d$pb == 1, "grp"]))
  # k is only one for independent testing data
  if(partitions == "independent") nk <- 1
  
  # if no partitions, return results without cv.stats
  if(nk == 0) {
    cv.res <- list(mod.full = mod.full, mod.full.pred = mod.full.pred, train.stats = train.stats.df, cv.stats = NULL)
    return(cv.res)
  }
  
  # list to contain cv statistics
  cv.stats <- list()
  
  for(k in 1:nk) {
    # assign partitions for training and testing occurrence data and for background data
    occs.train.vals <- d %>% dplyr::filter(pb == 1, grp != k) %>% dplyr::select(envs.names)
    occs.test.vals <- d %>% dplyr::filter(pb == 1, grp == k) %>% dplyr::select(envs.names)
    bg.train.vals <- d %>% dplyr::filter(pb == 0, grp != k) %>% dplyr::select(envs.names)
    bg.test.vals <- d %>% dplyr::filter(pb == 0, grp == k) %>% dplyr::select(envs.names)
    
    # if bg has a grp that is not partitioned (e.g., for randomkfold, bg is given "0"),
    # use the full background for each cross validation
    # if(nrow(bg.train.vals) == 0) bg.train.vals <- d %>% dplyr::filter(pb == 0) %>% dplyr::select(envs.names)
    
    # define model arguments for current model k
    mod.args <- enm@args(occs.train.vals, bg.train.vals, tune.i, settings$other.args)
    # run the current model k
    mod <- do.call(enm@fun, mod.args)
    
    # if model is NULL for some reason, continue but report to user
    if(is.null(mod)) {
      message(paste0("\nThe model for settings ", paste(names(tune.i), tune.i, collapse = ", "), " for partition ", k, " failed (resulted in NULL). Consider changing partitions. Cross validation averages will ignore this model.\n"))
      next
    }
    
    # calculate the stats for model k
    
    # calculate auc on testing data
    # NOTE: switch to bg.test??
    e.test <- enm@eval(occs.test.vals, bg.train.vals, mod, settings$other.args, settings$doClamp)
    auc.test <- e.test@auc
    # calculate auc diff
    auc.diff <- auc.train - auc.test
    if(settings$abs.auc.diff == TRUE) auc.diff <- abs(auc.diff)
    
    # get model predictions for training and testing data
    # these predictions are used only for calculating omission rate, and
    # thus should not need any specific parameter changes for maxent/maxnet
    d.vals <- d %>% dplyr::select(envs.names)
    d.pred <- d %>% dplyr::mutate(pred = enm@pred(mod.full, d.vals, settings$other.args, settings$doClamp, settings$pred.type))
    occs.train.pred <- d.pred %>% dplyr::filter(pb == 1, grp != k) %>% dplyr::pull(pred)
    occs.test.pred <- d.pred %>% dplyr::filter(pb == 1, grp == k) %>% dplyr::pull(pred)
    # get minimum training presence threshold (expected no omission)
    min.train.thr <- min(occs.train.pred)
    or.mtp <- mean(occs.test.pred < min.train.thr)
    # get 10 percentile training presence threshold (expected 0.1 omission)
    pct10.train.thr <- calc.10p.trainThresh(occs.train.vals, occs.train.pred)
    or.10p <- mean(occs.test.pred < pct10.train.thr)
    
    # calculate continuous Boyce Index
    if(settings$cvBoyce == TRUE) {
      occs.test.xy <- d %>% dplyr::filter(pb == 1, grp == k) %>% dplyr::select(1:2)
      cbi.test <- ecospat::ecospat.boyce(mod.full.pred, occs.test.xy, PEplot = FALSE)
    }else{
      cbi.test <- NULL
    }
    
    # # calculate MESS values if bg.test values are given
    # if(!is.null(bg.test.vals) & ncol(bg.test.vals) > 1) {
    #   mess.quant <- calc.mess.kstats(occs.train.vals, bg.train.vals, occs.test.vals, bg.test.vals)
    # }else{
    #   mess.quant <- NULL
    # }
    
    # gather all evaluation statistics for k
    
    kstats <- c(fold = k, auc.test = auc.test, auc.diff = auc.diff, or.mtp = or.mtp, 
                or.10p = or.10p, cbi.test = cbi.test$Spearman.cor)
    
    # add any additional cross-validation statistics chosen by the user
    kstats <- enm@kstats(kstats, e.test, mod, occs.train.vals, occs.test.vals, bg.train.vals, 
                         bg.test.vals, occs.train.pred, occs.test.pred, settings$other.args)
    # kstats <- c(kstats, mess.quant)
    # put into list as one-row data frame for easy binding
    cv.stats[[k]] <- data.frame(tune.args = tune.args.col, rbind(kstats), row.names=NULL, stringsAsFactors = FALSE)
  } 
  
  cv.stats.df <- dplyr::bind_rows(cv.stats)
  
  cv.res <- list(mod.full = mod.full, mod.full.pred = mod.full.pred, 
                 train.stats = train.stats.df, cv.stats = cv.stats.df)
  
  return(cv.res)
}
