#################################################
#########	MODEL TUNE for maxnet #############
#################################################

modelTune.maxnet <- function(pres, bg, env, nk, group.data, args, userArgs, 
                             rasterPreds, clamp, progbar, updateProgress) {
  
  if (progbar==TRUE) {
    pb <- txtProgressBar(0, length(args), style = 3)
  }
  out <- list()
  
  for (i in 1:length(args)) {
    # set up the console progress bar (progbar), 
    # or the shiny progress bar (updateProgress)
    if (length(args) > 1) {
      if (is.function(updateProgress)) {
        text <- paste0('Running ', args.lab[[1]][i], args.lab[[2]][i], '...')
        updateProgress(detail = text)
      } else if (progbar==TRUE) {
        setTxtProgressBar(pb, i)
      }
    }
    
    # set up data: x is coordinates of occs and bg, 
    # p is vector of 0's and 1's designating occs and bg
    x <- rbind(pres, bg)
    p <- c(rep(1, nrow(pres)), rep(0, nrow(bg)))
    
    # build the full model from all the data
      full.mod <- maxnet::maxnet(p, x, f=maxnet::maxnet.formula(p=p, data=x, classes=args[[i]][1]), 
                         regmult = as.numeric(args[[i]][2]))
    
    # if rasters selected, predict for the full model
    if (rasterPreds == TRUE) {
      predictive.map <- predict_maxnetRas(full.mod, env, type = 'exponential', clamp = clamp)
    } else {
      predictive.map <- stack()
    }
    
    # set up empty vectors for stats
    AUC.TEST <- double()
    AUC.DIFF <- double()
    OR10 <- double()
    ORmin <- double()
    
    # cross-validation on partitions
    for (k in 1:nk) {
      # set up training and testing data groups
      train.val <- pres[group.data$occ.grp != k, ]
      test.val <- pres[group.data$occ.grp == k, ]
      bg.val <- bg[group.data$bg.grp != k, ]
      # redefine x and p for partition groups
      x <- rbind(train.val, bg.val)
      p <- c(rep(1, nrow(train.val)), rep(0, nrow(bg.val)))
      
      # run the current test model
      mod <- maxnet::maxnet(p, x, f=maxnet::maxnet.formula(p=p, data=x, classes=args[[i]][1]), 
                            regmult = as.numeric(args[[i]][2]))
      
      AUC.TEST[k] <- dismo::evaluate(test.val, bg, mod)@auc
      AUC.DIFF[k] <- max(0, dismo::evaluate(train.val, bg, mod)@auc - AUC.TEST[k])
      
      # predict values for training and testing data
      p.train <- predict(mod, train.val, type = 'exponential')
      p.test <- predict(mod, test.val, type = 'exponential')  
      
      # figure out 90% of total no. of training records
      if (nrow(train.val) < 10) {
        n90 <- floor(nrow(train.val) * 0.9)
      } else {
        n90 <- ceiling(nrow(train.val) * 0.9)
      }
      train.thr.10 <- rev(sort(p.train))[n90]
      OR10[k] <- mean(p.test < train.thr.10)
      train.thr.min <- min(p.train)
      ORmin[k] <- mean(p.test < train.thr.min)
    }
    stats <- c(AUC.DIFF, AUC.TEST, OR10, ORmin)
    out[[i]] <- list(full.mod, stats, predictive.map)
  }
  if (progbar == TRUE) close(pb)
  return(out)
}

