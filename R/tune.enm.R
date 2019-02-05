#' @export

tune.enm <- function(occs.vals, bg.vals, occs.folds, bg.folds, envs, mod.fun, mod.name, 
                     tune.tbl, other.args, categoricals, doClamp, skipRasters,
                     abs.auc.diff, updateProgress) {
  results <- list()
  n <- nrow(tune.tbl)
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
    results[[i]] <- cv.enm(occs.vals, bg.vals, occs.folds, bg.folds, envs, mod.fun, mod.name, 
                           tune.tbl[i,], other.args, categoricals, doClamp, 
                           skipRasters, abs.auc.diff)
  }
  close(pb)   

  return(results)  
}

tune.enm.parallel <- function(occs.vals, bg.vals, occs.folds, bg.folds, envs, mod.fun, mod.name,
                              tune.tbl, other.args, categoricals, doClamp, skipRasters,
                              abs.auc.diff) {
  # set up parallel processing functionality
  allCores <- detectCores()
  if (is.null(numCores)) {
    numCores <- allCores
  }
  c1 <- makeCluster(numCores)
  registerDoParallel(c1)
  numCoresUsed <- getDoParWorkers()
  message(paste("Of", allCores, "total cores using", numCoresUsed))
  
  if(mod.name == "maxent") pkgs <- c("dismo", "raster", "ENMeval", "rJava")
  if(mod.name == "maxnet") pkgs <- c("dismo", "raster", "ENMeval", "maxnet")
  if(mod.name == "gbm") pkgs <- c("dismo", "raster", "ENMeval", "gbm")
  
  message("Running in parallel...")
  n <- nrow(tune.tbl)
  results <- foreach(i = 1:n, .packages = pkgs) %dopar% {
    cv.enm(occs.vals, bg.vals, occs.folds, bg.folds, envs, mod.fun, mod.name, tune.tbl[i,],
           other.args, categoricals, doClamp, skipRasters, abs.auc.diff)
  }
  stopCluster(c1)
  
  return(results)
}

cv.enm <- function(occs.vals, bg.vals, occs.folds, bg.folds, envs, mod.fun, mod.name, 
                   tune.tbl.i, other.args, categoricals, doClamp, skipRasters,
                   abs.auc.diff) {

  # build the full model from all the data
  mod.full.args <- make.args(tune.tbl.i, mod.name, occs.vals, bg.vals, other.args)
  mod.full <- do.call(mod.fun, mod.full.args)
  # calculate training auc
  auc.train <- calcAUC(occs.vals, bg.vals, mod.full, mod.name)
  
  # if rasters selected, predict for the full model
  if(skipRasters == FALSE) {
    mod.full.pred <- rasterPred(mod.full, envs, mod.name, doClamp)
  }else{
    mod.full.pred <- raster::stack()
  }
  
  # define number of folds (the value of "k")
  nk <- length(unique(occs.folds))
  
  # set up empty vectors for stats
  kstats <- data.frame(auc.test = numeric(nk), auc.diff = numeric(nk), 
                       or.mtp = numeric(nk), or.10p = numeric(nk))
  
  # cross-validation on partitions
  for(k in 1:nk) {
    # assign partitions for training and testing occurrence data and for background data
    train.k <- occs.vals[occs.folds != k,, drop = FALSE]
    test.k <- occs.vals[occs.folds == k,, drop = FALSE]
    bg.k <- bg.vals[bg.folds != k,, drop = FALSE]
    # define model arguments for current model k
    mod.k.args <- make.args(tune.tbl.i, mod.name, train.k, bg.k, other.args)
    # run the current model k
    mod.k <- do.call(mod.fun, mod.k.args)
    # calculate the stats for model k
    kstats[k,] <- evalStats(train.k, bg.k, test.k, mod.k, mod.name, abs.auc.diff)
  }
  
  cv.res <- list(mod.full = mod.full, mod.full.pred = mod.full.pred, 
                 kstats = kstats, auc.train = auc.train)
  return(cv.res)
}

evalStats <- function(occs.train, bg.train, occs.test, mod, mod.name, abs.auc.diff) {
  # calculate auc on training and testing data
  auc.train <- calcAUC(occs.train, bg.train, mod, mod.name)
  auc.test <- calcAUC(occs.test, bg.train, mod, mod.name)
  # calculate auc diff
  auc.diff <- auc.train - auc.test
  if(abs.auc.diff == TRUE) auc.diff <- abs(auc.diff)
  # get model predictions for training and testing data
  # these predictions are used only for calculating omission rate, and
  # thus should not need any specific parameter changes for maxent/maxnet
  pred.train <- vectorPred(mod, occs.train, mod.name, doClamp)
  pred.test <- vectorPred(mod, occs.test, mod.name, doClamp)
  # get 10 percentile predicted value
  occs.train.n <- nrow(occs.train)
  if(occs.train.n < 10) {
    pct10.train <- ceiling(occs.train.n * 0.1)
  } else {
    pct10.train <- floor(occs.train.n * 0.1)
  }
  pct10.train.thr <- sort(pred.train)[pct10.train]
  or.10p.test <- mean(pred.test < pct10.train.thr)
  min.train.thr <- min(pred.train)
  or.mtp.test <- mean(pred.test < min.train.thr)
  
  stats <- c(auc.test, auc.diff, or.mtp.test, or.10p.test)
  
  return(stats)
}

collateResults <- function(results, tune.tbl, mod.name, skipRasters) {
  # gather all full models into list
  mod.full.all <- lapply(results, function(x) x$mod.full)
  # gather all statistics into a data frame
  kstats.all <- lapply(results, function(x) x$kstats)
  if(skipRasters == FALSE) {
    mod.full.pred.all <- raster::stack(sapply(results, function(x) x$mod.full.pred))
  } else {
    mod.full.pred.all <- raster::stack()
  }
  
  # define a corrected variance function
  corrected.var <- function(x, nk){
    sum((x - mean(x))^2) * ((nk-1)/nk)
  }
  
  # make data frame of stats (for all folds)
  kstats.df <- dplyr::bind_rows(kstats.all)
  
  # define number of folds (the value of "k") as number of
  # rows in one of the model runs
  nk <- nrow(kstats.all[[1]])
  # define number of settings
  ns <- ncol(tune.tbl)
  # concatenate fc and rm to make settings column
  for(i in 1:ns) {
    kstats.df <- cbind(kstats.df, rep(tune.tbl[,i], each = nk))
    names(kstats.df)[ncol(kstats.df)] <- names(tune.tbl)[i]
  }
  # make new column for fold number
  kstats.df$fold <- rep(1:nk, nrow(tune.tbl))
  kstats.df <- tibble::as_tibble(kstats.df)
  
  # get number of columns in kstats.df
  nc <- ncol(kstats.df)
  # calculate number of non-zero parameters in model
  nparams <- sapply(mod.full.all, function(x) getNoParams(x, mod.name))
  
  # summarize by averaging all folds per model setting combination
  stats.df <- kstats.df %>% 
    dplyr::group_by_at(seq(nc-ns, nc-1)) %>% 
    dplyr::summarize(auc.test.mean = mean(auc.test),
                     auc.test.var = corrected.var(auc.test, nk),
                     auc.test.min = min(auc.test),
                     auc.test.max = max(auc.test),
                     auc.diff.mean = mean(auc.diff),
                     auc.diff.var = corrected.var(auc.diff, nk),
                     auc.diff.min = min(auc.diff),
                     auc.diff.max = max(auc.diff),
                     or.mtp.mean = mean(or.mtp),
                     or.mtp.var = corrected.var(or.mtp, nk),
                     or.mtp.min = min(or.mtp),
                     or.mtp.max = max(or.mtp),
                     or.10p.mean = mean(or.10p),
                     or.10p.var = corrected.var(or.10p, nk),
                     or.10p.min = min(or.10p),
                     or.10p.max = max(or.10p)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(auc.train = sapply(results, function(x) x$auc.train)) %>%
    dplyr::select(-((ns+1):(ns+16)), (ns+1):(ns+16)) %>%
    dplyr::bind_cols(calc.aicc(nparams, occs, mod.full.pred.all))
  
  # rearrange the columns for kstats
  kstats.df <- dplyr::select_at(kstats.df, c(seq(nc-ns, nc), 1:4))
  
  return(list(stats = stats.df, kstats = kstats.df, mods = mod.full.all,
              preds = mod.full.pred.all))
}


  
  # out.i <- list(mod.full = mod.full, mod.full.pred = mod.full.pred, stats = stats)
