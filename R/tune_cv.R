#' @export

modelTune.maxentJar <- function(occs.vals, bg.vals, envs, nk, occs.folds, bg.folds, args.i, userArgs, 
                                rasterPreds, clamp, categoricals) {
  # set up data: x is coordinates of occs and bg.vals, 
  # p is vector of 0's and 1's designating occs and bg.vals
  x <- rbind(occs.vals, bg.vals)
  p <- c(rep(1, nrow(occs.vals)), rep(0, nrow(bg.vals)))
  
  # build the full model from all the data
  mod.full <- dismo::maxent(x, p, args = c(args.i, userArgs),
                            factors = categoricals)  
  pred.args <- c("outputformat=raw", ifelse(clamp==TRUE, "doclamp=true", "doclamp=false"))
  
  # if rasters selected, predict for the full model
  if (rasterPreds == TRUE) {
    mod.full.pred <- predict(mod.full, envs, args = pred.args)  
  } else {
    mod.full.pred <- stack()
  }
  
  # set up empty vectors for stats
  AUC.TEST <- double()
  AUC.DIFF <- double()
  OR10 <- double()
  ORmin <- double()
  
  # define number of folds (the value of "k")
  nk <- length(unique(occs.folds))
  
  # cross-validation on partitions
  for (k in 1:nk) {
    # set up training and testing data groups
    train.k <- occs.vals[occs.folds != k,, drop = FALSE]
    test.k <- occs.vals[occs.folds == k,, drop = FALSE]
    bg.k <- bg.vals[bg.folds != k,, drop = FALSE]
    # redefine x and p for partition groups
    x <- rbind(train.k, bg.k)
    p <- c(rep(1, nrow(train.k)), rep(0, nrow(bg.k)))
    
    # run the current test model
    mod.k <- dismo::maxent(x, p, args = c(args.i, userArgs), factors = categoricals)  
    
    AUC.TEST[k] <- dismo::evaluate(test.k, bg.vals, mod.k)@auc
    AUC.DIFF[k] <- max(0, dismo::evaluate(train.k, bg.vals, mod.k)@auc - AUC.TEST[k])
    
    # predict values for training and testing data
    train.k.pred <- predict(mod.k, train.k, args = pred.args)
    test.k.pred <- predict(mod.k, test.k, args = pred.args)  
    
    # figure out 90% of total no. of training records
    if (nrow(train.k) < 10) {
      pct.90 <- floor(nrow(train.k) * 0.9)
    } else {
      pct.90 <- ceiling(nrow(train.k) * 0.9)
    }
    train.thr.10 <- rev(sort(train.k.pred))[pct.90]
    OR10[k] <- mean(test.k.pred < train.thr.10)
    train.thr.min <- min(train.k.pred)
    ORmin[k] <- mean(test.k.pred < train.thr.min)
  }
  stats <- c(AUC.DIFF, AUC.TEST, OR10, ORmin)
  out.i <- list(mod.full, stats, mod.full.pred)
  return(out.i)
}

#' @export

modelTune.maxnet <- function(occs.vals, bg.vals, envs, occs.folds, bg.folds, args.i,  
                             rasterPreds, clamp) {
  # set up data: x is coordinates of occs and bg.vals, 
  # p is vector of 0's and 1's designating occs and bg.vals
  x <- rbind(occs.vals, bg.vals)
  p <- c(rep(1, nrow(occs.vals)), rep(0, nrow(bg.vals)))
  
  # build the full model from all the data
  mod.full <- maxnet::maxnet(p, x, f=maxnet::maxnet.formula(p=p, data=x, classes=args.i[1]), 
                             regmult = as.numeric(args.i[2]))
  
  # if not skipping rasters, predict for the full model
  if (skipRasters == FALSE) {
    mod.full.pred <- maxnet.predictRaster(mod.full, envs, type = 'exponential', clamp = clamp)
  } else {
    mod.full.pred <- stack()
  }
  
  # define number of folds (the value of "k")
  nk <- length(unique(occs.folds))
  
  # set up empty vectors for stats
  stats <- data.frame(auc.test = numeric(nk), auc.diff = numeric(nk), 
                or.min = numeric(nk), or.10 = numeric(nk))
  
  # cross-validation on partitions
  for(k in 1:nk) {
    # assign partitions for training and testing occurrence data and for background data
    train.k <- occs.vals[occs.folds != k,, drop = FALSE]
    test.k <- occs.vals[occs.folds == k,, drop = FALSE]
    bg.k <- bg.vals[bg.folds != k,, drop = FALSE]
    # define x and p for this fold
    x.k <- rbind(train.k, bg.k)
    p.k <- c(rep(1, nrow(train.k)), rep(0, nrow(bg.k)))
    
    # run the current test model
    mod.k <- maxnet::maxnet(p.k, x.k, f=maxnet::maxnet.formula(p=p.k, data=x.k, classes=args.i[1]), 
                          regmult = as.numeric(args.i[2]))
    
    stats$auc.test[k] <- dismo::evaluate(test.k, bg.vals, mod.k)@auc
    stats$auc.diff[k] <- max(0, dismo::evaluate(train.k, bg.vals, mod.k)@auc - stats$auc.test[k])
    
    # predict values for training and testing data
    train.k.pred <- dismo::predict(mod.k, train.k, type = 'exponential')
    test.k.pred <- dismo::predict(mod.k, test.k, type = 'exponential')  
    
    # figure out 90th percentile of total no. of training records
    pct.10 <- floor(nrow(train.k) * 0.1)
    sort(train.k.pred)[pct.10]
    
    
    if(nrow(train.k) < 10) {
      pct.10 <- ceiling(nrow(train.k) * 0.1)
    } else {
      pct.10 <- floor(nrow(train.k) * 0.1)
    }
    pct.10.thr <- sort(train.k.pred)[pct.10]
    stats$or.10[k] <- mean(test.k.pred < pct.10.thr)
    min.thr <- min(train.k.pred)
    stats$or.min[k] <- mean(test.k.pred < min.thr)
  }
  
  out.i <- list(mod.full = mod.full, mod.full.pred = mod.full.pred, stats = stats)
  return(out.i)
}
