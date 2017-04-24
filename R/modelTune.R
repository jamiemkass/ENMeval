#################################################
#########	MODEL TUNE #############
#################################################

modelTune <- function(pres, bg, env, nk, group.data, progbar, maxent.args, 
                      userArgs, rasterPreds, clamp, java, updateProgress) {
  
  if (progbar==TRUE) {
    pb <- txtProgressBar(0, length(maxent.args), style = 3)
  }
  out <- list()
  
  for (i in 1:length(maxent.args)) {
    # set up the console progress bar (progbar), 
    # or the shiny progress bar (updateProgress)
    if (length(maxent.args) > 1) {
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
    if (java == FALSE) {
      full.mod <- maxnet::maxnet(p, x, f=maxnet::maxnet.formula(p=p, data=x, classes=maxent.args[[i]][1]), 
                         regmult = as.numeric(maxent.args[[i]][2]))
    } else {
      # set up temp folder to delete later
      tmpfolder <- tempfile()
      full.mod <- dismo::maxent(x, p, args = c(maxent.args[[i]], userArgs),
                         factors = categoricals, path = tmpfolder)  
    }
    
    # if rasters selected, predict for the full model
    if (rasterPreds == TRUE) {
      if (java == FALSE) {
        predictive.map <- predict_maxnetRas(full.mod, env, type = 'exponential', clamp = clamp)
      } else {
        pred.args <- c("outputformat=raw", ifelse(clamp==TRUE, "doclamp=true", "doclamp=false"))
        predictive.map <- predict(full.mod, env, args = pred.args)  
      }
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
      if (java == FALSE) {
        mod <- maxnet::maxnet(p, x, f=maxnet::maxnet.formula(p=p, data=x, classes=maxent.args[[i]][1]), 
                           regmult = as.numeric(maxent.args[[i]][2]))
      } else {
        mod <- dismo::maxent(x, p, args = c(maxent.args[[i]], userArgs), factors = categoricals,
                      path = tmpfolder)  
      }
      
      AUC.TEST[k] <- dismo::evaluate(test.val, bg, mod)@auc
      AUC.DIFF[k] <- max(0, dismo::evaluate(train.val, bg, mod)@auc - AUC.TEST[k])
      
      # predict values for training and testing data
      if (java == FALSE) {
        p.train <- predict(mod, train.val, type = 'exponential')
        p.test <- predict(mod, test.val, type = 'exponential')  
      } else {
        p.train <- predict(mod, train.val, args = pred.args)
        p.test <- predict(mod, test.val, args = pred.args)  
      }
      
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
    if (java == TRUE) unlink(tmpfolder, recursive = TRUE)
    stats <- c(AUC.DIFF, AUC.TEST, OR10, ORmin)
    out[[i]] <- list(full.mod, stats, predictive.map)
  }
  if (progbar == TRUE) close(pb)
  return(out)
}

