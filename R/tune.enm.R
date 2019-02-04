#' @export

tune.enm <- function(occs.vals, bg.vals, occs.folds, bg.folds, envs, mod.fun, mod.fun.name,
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
    results[[i]] <- cv.enm(occs.vals, bg.vals, occs.folds, bg.folds, envs, mod.fun, mod.fun.name,
                           tune.tbl[i,], other.args, categoricals, doClamp, 
                           skipRasters, abs.auc.diff)
  }
  close(pb)   

  return(results)  
}

tune.enm.parallel <- function(occs.vals, bg.vals, occs.folds, bg.folds, envs, mod.fun, mod.fun.name,
                              tune.tbl, other.args, categoricals, doClamp, skipRasters,
                              abs.auc.diff) {
  # get model function's name
  mod.fun.name <- as.character(substitute(mod.fun))[3]
  # set up parallel processing functionality
  allCores <- detectCores()
  if (is.null(numCores)) {
    numCores <- allCores
  }
  c1 <- makeCluster(numCores)
  registerDoParallel(c1)
  numCoresUsed <- getDoParWorkers()
  message(paste("Of", allCores, "total cores using", numCoresUsed))
  
  if(mod.fun.name == "maxent") pkgs <- c("dismo", "raster", "ENMeval", "rJava")
  if(mod.fun.name == "maxnet") pkgs <- c("dismo", "raster", "ENMeval", "maxnet")
  
  message("Running in parallel...")
  n <- nrow(tune.tbl)
  results <- foreach(i = 1:n, .packages = pkgs) %dopar% {
    cv.enm(occs.vals, bg.vals, occs.folds, bg.folds, envs, mod.fun, mod.fun.name, tune.tbl[i,],
           other.args, categoricals, doClamp, skipRasters, abs.auc.diff)
  }
  stopCluster(c1)
  
  return(results)
}

cv.enm <- function(occs.vals, bg.vals, occs.folds, bg.folds, envs, mod.fun, mod.fun.name,
                   tune.tbl.i, other.args, categoricals, doClamp, skipRasters,
                   abs.auc.diff) {
  # get model function's name
  mod.fun.name <- as.character(substitute(mod.fun))[3]
  print(mod.fun.name)
  
  # build the full model from all the data
  mod.full.args <- make.args(tune.tbl.i, mod.fun.name, occs.vals, bg.vals, other.args)
  mod.full <- do.call(mod.fun, mod.full.args)
  auc.train.full <- dismo::evaluate(occs.vals, bg.vals, mod.full)@auc
  
  # if rasters selected, predict for the full model
  if(skipRasters == FALSE) {
    if(mod.fun.name == "maxent") {
      pred.args <- c("outputformat=raw", ifelse(doClamp==TRUE, "doclamp=true", "doclamp=false"))
      mod.full.pred <- dismo::predict(mod.full, envs, args = pred.args)
    }
    else if(mod.fun.name == "maxnet") {
      mod.full.pred <- maxnet.predictRaster(mod.full, envs, type = 'exponential', clamp = doClamp)
    }
    else{
      mod.full.pred <- predict(mod.full, envs)
    }
  } 
  else{
    mod.full.pred <- raster::stack()
  }
  
  # define number of folds (the value of "k")
  nk <- length(unique(occs.folds))
  
  # set up empty vectors for stats
  kstats <- data.frame(auc.test = numeric(nk), auc.diff = numeric(nk), 
                      or.min = numeric(nk), or.10 = numeric(nk))
  
  # cross-validation on partitions
  for(k in 1:nk) {
    # assign partitions for training and testing occurrence data and for background data
    train.k <- occs.vals[occs.folds != k,, drop = FALSE]
    test.k <- occs.vals[occs.folds == k,, drop = FALSE]
    bg.k <- bg.vals[bg.folds != k,, drop = FALSE]
    # define model arguments for current model k
    mod.k.args <- make.args(tune.tbl.i, mod.fun.name, train.k, bg.k, other.args)
    # run the current model k
    mod.k <- do.call(mod.fun, mod.k.args)
    # calculate the stats for model k
    kstats[k,] <- evalStats(train.k, bg.k, test.k, mod.k, abs.auc.diff)
  }
  
  cv.res <- list(mod.full = mod.full, mod.full.pred = mod.full.pred, 
                 kstats = kstats, auc.train.full = auc.train.full)
  return(cv.res)
}

evalStats <- function(occs.train, bg.train, occs.test, mod, abs.auc.diff) {
  # calculate auc on training and testing data
  auc.train <- dismo::evaluate(occs.train, bg.train, mod)@auc
  auc.test <- dismo::evaluate(occs.test, bg.train, mod)@auc
  # calculate auc diff
  auc.diff <- auc.train - auc.test
  if(abs.auc.diff == TRUE) auc.diff <- abs(auc.diff)
  # get model predictions for training and testing data
  # these predictions are used only for calculating omission rate, and
  # thus should not need any specific parameter changes for maxent/maxnet
  pred.train <- predict(mod, occs.train)
  pred.test <- predict(mod, occs.test)
  # get 10 percentile predicted value
  occs.train.n <- nrow(occs.train)
  if(occs.train.n < 10) {
    pct10.train <- ceiling(occs.train.n * 0.1)
  } else {
    pct10.train <- floor(occs.train.n * 0.1)
  }
  pct10.train.thr <- sort(pred.train)[pct10.train]
  or10.test <- mean(pred.test < pct10.train.thr)
  min.train.thr <- min(pred.train)
  orMin.test <- mean(pred.test < min.train.thr)
  
  stats <- c(auc.test, auc.diff, orMin.test, or10.test)
  
  return(stats)
}


  
  # out.i <- list(mod.full = mod.full, mod.full.pred = mod.full.pred, stats = stats)
