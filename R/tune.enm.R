#' @export

cv.enm <- function(occs.vals, bg.vals, occs.grp, bg.grp, envs, ls, 
                   partitions, tune.tbl.i, other.args, categoricals, occs.ind, 
                   doClamp, skipRasters, abs.auc.diff) {
  
  # build the full model from all the data
  mod.full.args <- ls$args(occs.vals, bg.vals, tune.tbl.i, other.args)
  mod.full <- do.call(ls$fun, mod.full.args)
  # calculate training auc
  auc.train <- ls$calcAUC(occs.vals, bg.vals, mod.full, other.args, doClamp)
  
  # if rasters selected and envs is not NULL, predict raster for the full model
  if(skipRasters == FALSE & !is.null(envs)) {
    mod.full.pred <- ls$predict(mod.full, envs, other.args, doClamp)
    names(mod.full.pred) <- paste(tune.tbl.i, collapse = "")
  }else{
    mod.full.pred <- raster::stack()
  }
  
  # define number of grp (the value of "k")
  nk <- length(unique(occs.grp))
  
  # set up empty vectors for stats
  cnames <- c("fold", "auc.test", "auc.diff", "or.mtp", "or.10p")
  # kstats <- as.data.frame(matrix(nrow = nk, ncol = length(cnames), 
  #                                dimnames = list(rep("", nk), cnames)), row.names = FALSE)
  kstats.ls <- list()
  
  # if there are no grp specified...
  if(nk == 0) {
    # if user selects to use independent testing data, do not do k-fold cross validation
    if(partitions == "independent") {
      occs.ind.vals <- as.data.frame(raster::extract(envs, occs.ind))
      auc.test <- ls$calcAUC(occs.ind.vals, bg.vals, mod.full, other.args, doClamp)
      e <- evalStats(occs.vals, bg.vals, occs.ind.vals, bg.test = NULL, ls,
                     auc.train, mod.full, other.args, doClamp, abs.auc.diff)
      kstats.ls[[1]] <- c(fold = 1, e)
    }
    # # if user selects to only calculate AICc, stop here
    # if(partitions == "none") break
  }else{
    # cross-validation on partitions
    for(k in 1:nk) {
      # assign partitions for training and testing occurrence data and for background data
      occs.train.k <- occs.vals[occs.grp != k,, drop = FALSE]
      occs.test.k <- occs.vals[occs.grp == k,, drop = FALSE]
      bg.train.k <- bg.vals[bg.grp != k,, drop = FALSE]
      bg.test.k <- bg.vals[bg.grp == k,, drop = FALSE]
      # define model arguments for current model k
      mod.k.args <- ls$args(occs.train.k, bg.train.k, tune.tbl.i, other.args)
      # run the current model k
      mod.k <- do.call(ls$fun, mod.k.args)
      # calculate the stats for model k
      e <- evalStats(occs.train.k, bg.vals, occs.test.k, bg.test.k, ls,
                     auc.train, mod.k, categoricals, other.args, doClamp, abs.auc.diff)
      kstats.ls[[k]] <- c(fold = k, e)
    } 
  }
  
  kstats <- as.data.frame(do.call("rbind", kstats.ls))
  
  cv.res <- list(mod.full = mod.full, mod.full.pred = mod.full.pred, 
                 kstats = kstats, train.AUC = auc.train)
  
  return(cv.res)
}

evalStats <- function(occs.train, bg.train, occs.test, bg.test, ls, auc.train, mod, categoricals, other.args, doClamp, abs.auc.diff) {
  # calculate auc on testing data
  auc.test <- ls$calcAUC(occs.test, bg.train, mod, other.args, doClamp)
  # calculate auc diff
  auc.diff <- auc.train - auc.test
  if(abs.auc.diff == TRUE) auc.diff <- abs(auc.diff)
  # get model predictions for training and testing data
  # these predictions are used only for calculating omission rate, and
  # thus should not need any specific parameter changes for maxent/maxnet
  pred.train <- ls$predict(mod, occs.train, other.args, doClamp)
  pred.test <- ls$predict(mod, occs.test, other.args, doClamp)
  # get minimum training presence threshold (expected no omission)
  min.train.thr <- min(pred.train)
  or.mtp <- mean(pred.test < min.train.thr)
  # get 10 percentile training presence threshold (expected 0.1 omission)
  occs.train.n <- nrow(occs.train)
  if(occs.train.n < 10) {
    pct10.train <- floor(occs.train.n * 0.1)
  } else {
    pct10.train <- ceiling(occs.train.n * 0.1)
  }
  pct10.train.thr <- sort(pred.train)[pct10.train]
  or.10p <- mean(pred.test < pct10.train.thr)
  
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
    mess.quant <- quantile(mss)
    names(mess.quant) <- paste0("mess_", names(mess.quant))
  }else{
    mess.quant <- NULL
  }
  
  stats <- c(auc.test = auc.test, auc.diff = auc.diff, or.mtp = or.mtp, 
             or.10p = or.10p, mess.quant)
  
  return(stats)
}
