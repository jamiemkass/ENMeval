#' @export

cv.enm <- function(occs.vals, bg.vals, occs.folds, bg.folds, envs, mod.fun, mod.name, 
                   partitions, tune.tbl.i, other.args, categoricals, occs.ind, 
                   doClamp, skipRasters, abs.auc.diff) {
  
  # build the full model from all the data
  mod.full.args <- model.args(tune.tbl.i, mod.name, occs.vals, bg.vals, other.args)
  mod.full <- do.call(mod.fun, mod.full.args)
  # calculate training auc
  auc.train <- calcAUC(occs.vals, bg.vals, mod.full, mod.name)
  
  # if rasters selected and envs is not NULL, predict raster for the full model
  if(skipRasters == FALSE & !is.null(envs)) {
    mod.full.pred <- rasterPred(mod.full, envs, mod.name, other.args, doClamp)
  }else{
    mod.full.pred <- raster::stack()
  }
  
  # define number of folds (the value of "k")
  nk <- length(unique(occs.folds))
  
  # set up empty vectors for stats
  cnames <- c("fold", "auc.test", "auc.diff", "or.mtp", "or.10p", "mess.mean", "mess.min")
  kstats <- as.data.frame(matrix(nrow = nk, ncol = length(cnames), 
                                 dimnames = list(rep("", nk), cnames)), row.names = FALSE)
  
  # if there are no folds specified...
  if(nk == 0) {
    # if user selects to use independent testing data, do not do k-fold cross validation
    if(partitions == "independent") {
      occs.ind.vals <- as.data.frame(raster::extract(envs, occs.ind))
      auc.test <- calcAUC(occs.ind.vals, bg.vals, mod.full, mod.name)
      kstats[1,] <- c(1, evalStats(occs.vals, bg.vals, occs.ind.vals, bg.test = NULL, 
                              auc.train, mod.full, mod.name, other.args, doClamp, abs.auc.diff))
    }
    # # if user selects to only calculate AICc, stop here
    # if(partitions == "none") break
  }else{
    # cross-validation on partitions
    for(k in 1:nk) {
      # assign partitions for training and testing occurrence data and for background data
      occs.train.k <- occs.vals[occs.folds != k,, drop = FALSE]
      occs.test.k <- occs.vals[occs.folds == k,, drop = FALSE]
      bg.train.k <- bg.vals[bg.folds != k,, drop = FALSE]
      bg.test.k <- bg.vals[bg.folds == k,, drop = FALSE]
      # define model arguments for current model k
      mod.k.args <- model.args(tune.tbl.i, mod.name, occs.train.k, bg.train.k, other.args)
      # run the current model k
      mod.k <- do.call(mod.fun, mod.k.args)
      # calculate the stats for model k
      # kstats[k,] <- evalStats(occs.train.k, bg.train.k, occs.test.k, bg.test.k,
      #                         auc.train, mod.k, mod.name, doClamp, abs.auc.diff)
      kstats[k,] <- c(k, evalStats(occs.train.k, bg.vals, occs.test.k, bg.test.k,
                              auc.train, mod.k, mod.name, other.args, doClamp, abs.auc.diff))
    } 
  }
  
  cv.res <- list(mod.full = mod.full, mod.full.pred = mod.full.pred, 
                 kstats = kstats, auc.train = auc.train)
  
  return(cv.res)
}

evalStats <- function(occs.train, bg.train, occs.test, bg.test, auc.train, mod, mod.name, other.args, doClamp, abs.auc.diff) {
  # calculate auc on testing data
  auc.test <- calcAUC(occs.test, bg.train, mod, mod.name)
  # calculate auc diff
  auc.diff <- auc.train - auc.test
  if(abs.auc.diff == TRUE) auc.diff <- abs(auc.diff)
  # get model predictions for training and testing data
  # these predictions are used only for calculating omission rate, and
  # thus should not need any specific parameter changes for maxent/maxnet
  pred.train <- vectorPred(mod, occs.train, mod.name, other.args, doClamp)
  pred.test <- vectorPred(mod, occs.test, mod.name, other.args, doClamp)
  # get 10 percentile predicted value
  occs.train.n <- nrow(occs.train)
  if(occs.train.n < 10) {
    pct10.train <- floor(occs.train.n * 0.1)
  } else {
    pct10.train <- ceiling(occs.train.n * 0.1)
  }
  pct10.train.thr <- sort(pred.train)[pct10.train]
  or.10p.test <- mean(pred.test < pct10.train.thr)
  min.train.thr <- min(pred.train)
  or.mtp.test <- mean(pred.test < min.train.thr)
  
  # calculate MESS values if bg.test values are given
  if(!is.null(bg.test)) {
    p <- rbind(occs.train, bg.train)
    v <- rbind(occs.test, bg.test)
    cat.i <- which(names(occs.train) == categoricals)
    if(length(cat.i) > 0) {
      p <- p[,-cat.i]
      v <- v[,-cat.i]
    }
    mss <- mess.vec(p, v)
    mess.mean <- mean(mss)
    mess.min <- min(mss)  
  }else{
    mess.mean <- NA
    mess.min <- NA  
  }
  
  stats <- c(auc.test, auc.diff, or.mtp.test, or.10p.test, mess.mean, mess.min)
  
  return(stats)
}

# out.i <- list(mod.full = mod.full, mod.full.pred = mod.full.pred, stats = stats)
