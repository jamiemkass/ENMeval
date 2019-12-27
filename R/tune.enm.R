#' @title Iterate tuning of ENMs
#' @description Internal functions to tune and summarize results for ecological niche models (ENMs) iteratively across a range of user-specified tuning settings. 
#' Function \code{tune.parallel()} tunes ENMs with parallelization. Function \code{cv.enm()} calculates training and testing evaluation statistics for one set of specified tuning parameters.
#' @aliases tune.parallel tune.regular cv.enm
#' @param d data frame from \code{ENMevaluate()} with occurrence and background coordinates (or coordinates plus predictor variable values) and partition group values
#' @param envs Raster* object of environmental variables (must be in same geographic projection as occurrence data)
#' @param enm Object of class \link{ENMdetails}.
#' @param partitions character of name of partitioning technique (see \code{?partitions})
#' @param tune.tbl Data frame of tuning parameter combinations.
#' @param other.settings list of settings from \code{ENMevaluate()} containing other.args, doClamp, pred.type, abs.auc.diff, cbi.cv, cbi.eval
#' @param numCores boolean (TRUE or FALSE); if TRUE, use specifed number of cores for parallel processing
#' @param parallelType character of either "doParallel" or "doSNOW" to conduct parallelization

#' @name tune.enm
NULL

#' @rdname tune.enm
tune.parallel <- function(d, envs, enm, partitions, tune.tbl, other.settings, numCores, parallelType) {
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
    cv.enm(d, envs, enm, partitions, tune.tbl[i,], other.settings)
  }
  close(pb)
  parallel::stopCluster(cl)
  return(results)
}

#' @rdname tune.enm
tune.regular <- function(d, envs, enm, partitions, tune.tbl, other.settings, updateProgress) {
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
    results[[i]] <- cv.enm(d, envs, enm, partitions, tune.i, other.settings)
  }
  close(pb)
  return(results)
}

#' @param tune.i vector of single set of tuning parameters

#' @rdname tune.enm
cv.enm <- function(d, envs, enm, partitions, tune.i, other.settings) {
  envs.names <- names(d[, 3:(ncol(d)-2)])
  # unpack predictor variable values for occs and bg
  d.vals <- d %>% dplyr::select(pb, envs.names)
  occs.vals <- d.vals %>% dplyr::filter(pb == 1) %>% dplyr::select(envs.names)
  bg.vals <- d.vals %>% dplyr::filter(pb == 0) %>% dplyr::select(envs.names)
  # build the full model from all the data
  mod.full.args <- enm@args(occs.vals, bg.vals, tune.i, other.settings$other.args)
  mod.full <- do.call(enm@fun, mod.full.args)
  if(is.null(mod.full)) stop("Training model is NULL. Consider changing the tuning parameters.")
  # calculate training auc
  e.train <- enm@eval(occs.vals, bg.vals, mod.full, other.settings$other.args, other.settings$doClamp)
  auc.train <- e.train@auc
  tune.args.col <- paste(tune.i, collapse = "_")
  train.stats.df <- data.frame(tune.args = tune.args.col, auc.train = auc.train, stringsAsFactors = FALSE)
  # if envs is not NULL, predict raster for the full model and calculate CBI.train on this raster
  if(!is.null(envs)) {
    mod.full.pred <- enm@pred(mod.full, envs, other.settings$other.args, other.settings$doClamp, other.settings$pred.type)
    # training CBI
    d.occs.xy <- d %>% dplyr::filter(pb == 1) %>% dplyr::select(1:2)
    cbi.train <- ecospat::ecospat.boyce(mod.full.pred, d.occs.xy, PEplot = FALSE)
    train.stats.df$cbi.train <- cbi.train$Spearman.cor
  }else{
    # if envs is NULL, calculate CBI.train with the occurrence + background points
    d.full.pred <- d %>% dplyr::mutate(pred = enm@pred(mod.full, d.vals %>% dplyr::select(envs.names), other.settings$other.args, other.settings$doClamp, other.settings$pred.type))
    occs.full.pred <- d.full.pred %>% dplyr::filter(pb == 1)
    cbi.train <- ecospat::ecospat.boyce(d.full.pred$pred, occs.full.pred$pred, PEplot = FALSE)
    mod.full.pred <- d.full.pred$pred
  }
    
  # define number of grp (the value of "k") for occurrences
  # k is only one for independent testing data
  nk <- ifelse(partitions == "independent", 1, length(unique(d[d$pb == 1, "grp"])))
  
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
    
    # define model arguments for current model k
    mod.k.args <- enm@args(occs.train.vals, bg.train.vals, tune.i, other.settings$other.args)
    # run the current model k
    mod.k <- do.call(enm@fun, mod.k.args)
    
    # if model is NULL for some reason, continue but report to user
    if(is.null(mod.k)) {
      message(paste0("\nThe model for settings ", paste(names(tune.i), tune.i, collapse = ", "), " for partition ", k, " failed (resulted in NULL). Consider changing partitions. Cross validation averages will ignore this model.\n"))
      next
    }
    
    # calculate the stats for model k
    
    # calculate auc on testing data: test occurrences are evaluated on full background, as in Radosavljevic & Anderson 2014
    # NOTE: switch to bg.test??
    e.test <- enm@eval(occs.test.vals, bg.vals, mod.k, other.settings$other.args, other.settings$doClamp)
    auc.test <- e.test@auc
    # calculate auc diff
    auc.diff <- auc.train - auc.test
    if(other.settings$abs.auc.diff == TRUE) auc.diff <- abs(auc.diff)
    
    # get model predictions for training and testing data
    # these predictions are used only for calculating omission rate, and
    # thus should not need any specific parameter changes for maxent/maxnet
    d.vals <- d %>% dplyr::select(envs.names)
    d.pred <- d %>% dplyr::mutate(pred = enm@pred(mod.k, d.vals, other.settings$other.args, other.settings$doClamp, other.settings$pred.type))
    occs.train.pred <- d.pred %>% dplyr::filter(pb == 1, grp != k) %>% dplyr::pull(pred) %>% as.numeric()
    occs.test.pred <- d.pred %>% dplyr::filter(pb == 1, grp == k) %>% dplyr::pull(pred) %>% as.numeric()
    # get minimum training presence threshold (expected no omission)
    min.train.thr <- min(occs.train.pred)
    or.mtp <- mean(occs.test.pred < min.train.thr)
    # get 10 percentile training presence threshold (expected 0.1 omission)
    pct10.train.thr <- calc.10p.trainThresh(occs.train.pred)
    or.10p <- mean(occs.test.pred < pct10.train.thr)
    
    # calculate continuous Boyce Index
    if(other.settings$cbi.cv == TRUE) {
      if(!is.null(envs) & other.settings$cbi.eval == "envs") {
        # use full model prediction over envs
        mod.k.pred <- enm@pred(mod.k, envs, other.settings$other.args, other.settings$doClamp, other.settings$pred.type)
        # input test occs are coordinates
        occs.test.in <- d %>% dplyr::filter(pb == 1, grp == k) %>% dplyr::select(1:2)
      }else{
        # use full background to approximate full model prediction
        mod.k.pred <- d.pred %>% dplyr::filter(pb == 0) %>% dplyr::pull(pred) %>% as.numeric()
        # input test occs are values
        occs.test.in <- occs.test.pred
      }
      cbi.test <- ecospat::ecospat.boyce(mod.k.pred, occs.test.in, PEplot = FALSE)
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
                or.10p = or.10p, cbi.test = cbi.test$Spearman.cor,
                # add any additional cross-validation statistics chosen by the user
                enm@kstats(e.test, mod.k, other.settings$other.args))
    
    # put into list as one-row data frame for easy binding
    cv.stats[[k]] <- data.frame(tune.args = tune.args.col, rbind(kstats), row.names=NULL, stringsAsFactors = FALSE)
  } 
  
  cv.stats.df <- dplyr::bind_rows(cv.stats)
  
  cv.res <- list(mod.full = mod.full, mod.full.pred = mod.full.pred, 
                 train.stats = train.stats.df, cv.stats = cv.stats.df)
  
  return(cv.res)
}
