#' @title Functions for tuning ENMs
#' @description Internal functions to tune and summarize results for ecological niche models (ENMs) iteratively across a range of user-specified tuning settings.
#' @aliases tune.parallel tune.regular cv.enm evalStats
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
tune.parallel <- function(occs.vals, bg.vals, occs.grp, bg.grp, envs, enm, 
                          partitions, tune.tbl, other.args, categoricals, 
                          occs.ind, doClamp, skipRasters, abs.auc.diff, numCores) {
  # set up parallel processing functionality
  allCores <- parallel::detectCores()
  if (is.null(numCores)) {
    numCores <- allCores
  }
  cl <- parallel::makeCluster(numCores)
  doSNOW::registerDoSNOW(cl)
  numCoresUsed <- foreach::getDoParWorkers()
  message(paste0("Of ", allCores, " total cores using ", numCoresUsed, "..."))
  
  message("Running in parallel...")
  n <- ifelse(nrow(tune.tbl) > 0, nrow(tune.tbl), 1)
  pb <- txtProgressBar(0, n, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  results <- foreach::foreach(i = 1:n, .packages = enm.pkgs(enm), .options.snow = opts) %dopar% {
    cv.enm(occs.vals, bg.vals, occs.grp, bg.grp, envs, enm,
           partitions, tune.tbl[i,], other.args, categoricals, 
           occs.ind, doClamp, skipRasters, abs.auc.diff)
  }
  close(pb)
  parallel::stopCluster(cl)
  return(results)
}

#' @rdname tune.enm
tune.regular <- function(occs.vals, bg.vals, occs.grp, bg.grp, envs, enm, 
                         partitions, tune.tbl, other.args, categoricals, 
                         occs.ind, doClamp, skipRasters, abs.auc.diff, updateProgress) {
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
    results[[i]] <- cv.enm(occs.vals, bg.vals, occs.grp, bg.grp, envs, enm,
                           partitions, tune.settings = tune.tbl[i,], other.args, 
                           categoricals, occs.ind, doClamp, skipRasters, abs.auc.diff)
  }
  close(pb)
  return(results)
}

#' @rdname tune.enm
cv.enm <- function(occs.vals, bg.vals, occs.grp, bg.grp, envs, enm, 
                   partitions, tune.settings, other.args, categoricals, 
                   occs.ind, doClamp, skipRasters, abs.auc.diff) {
  
  # build the full model from all the data
  mod.full.args <- enm@args(occs.vals, bg.vals, tune.settings, other.args)
  mod.full <- do.call(enm@fun, mod.full.args)
  # calculate training auc
  auc.train <- enm@auc(occs.vals, bg.vals, mod.full, other.args, doClamp)
  
  # if rasters selected and envs is not NULL, predict raster for the full model
  if(skipRasters == FALSE & !is.null(envs)) {
    mod.full.pred <- enm@pred(mod.full, envs, other.args, doClamp)
  }else{
    mod.full.pred <- raster::stack()
  }
  
  # define number of grp (the value of "k")
  nk <- length(unique(occs.grp))
  
  # set up empty vectors for stats
  cnames <- c("fold", "auc.test", "auc.diff", "or.mtp", "or.10p")
  # kstats <- as.data.frame(matrix(nrow = nk, ncol = length(cnames), 
  #                                dimnames = list(rep("", nk), cnames)), row.names = FALSE)
  kstats.enm <- list()
  
  # if there are no grp specified...
  if(nk == 0) {
    # if user selects to use independent testing data, do not do k-fold cross validation
    if(partitions == "independent") {
      occs.ind.vals <- as.data.frame(raster::extract(envs, occs.ind))
      auc.test <- enm@auc(occs.ind.vals, bg.vals, mod.full, other.args, doClamp)
      e <- evalStats(occs.vals, bg.vals, occs.ind.vals, bg.test = NULL, enm,
                     auc.train, mod.full, categoricals, other.args, doClamp, abs.auc.diff)
      kstats.enm[[1]] <- c(fold = 1, e)
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
      mod.k.args <- enm@args(occs.train.k, bg.train.k, tune.settings, other.args)
      # run the current model k
      mod.k <- do.call(enm@fun, mod.k.args)
      # calculate the stats for model k
      e <- evalStats(occs.train.k, bg.vals, occs.test.k, bg.test.k, enm,
                     auc.train, mod.k, categoricals, other.args, doClamp, abs.auc.diff)
      kstats.enm[[k]] <- c(fold = k, e)
    } 
  }
  
  kstats <- as.data.frame(do.call("rbind", kstats.enm))
  
  cv.res <- list(mod.full = mod.full, mod.full.pred = mod.full.pred, 
                 kstats = kstats, train.AUC = auc.train)
  
  return(cv.res)
}

#' @rdname tune.enm
evalStats <- function(occs.train, bg.train, occs.test, bg.test, enm, auc.train, 
                      mod, categoricals, other.args, doClamp, abs.auc.diff) {
  # calculate auc on testing data
  auc.test <- enm@auc(occs.test, bg.train, mod, other.args, doClamp)
  # calculate auc diff
  auc.diff <- auc.train - auc.test
  if(abs.auc.diff == TRUE) auc.diff <- abs(auc.diff)
  # get model predictions for training and testing data
  # these predictions are used only for calculating omission rate, and
  # thus should not need any specific parameter changes for maxent/maxnet
  pred.train <- enm@pred(mod, occs.train, other.args, doClamp)
  pred.test <- enm@pred(mod, occs.test, other.args, doClamp)
  # get minimum training presence threshold (expected no omission)
  min.train.thr <- min(pred.train)
  or.mtp <- mean(pred.test < min.train.thr)
  # get 10 percentile training presence threshold (expected 0.1 omission)
  occs.train.n <- nrow(occs.train)
  if (occs.train.n < 10) {
    pct90.train <- floor(occs.train.n * 0.9)
  } else {
    pct90.train <- ceiling(occs.train.n * 0.9)
  }
  pct10.train.thr <- rev(sort(pred.train))[pct90.train]
  or.10p <- mean(pred.test < pct10.train.thr)
  
  # calculate MESS values if bg.test values are given
  if(!is.null(bg.test) & ncol(bg.test) > 1) {
    p <- rbind(occs.train, bg.train)
    v <- rbind(occs.test, bg.test)
    cat.i <- which(names(occs.train) == categoricals)
    if(length(cat.i) > 0) {
      p <- p[,-cat.i]
      v <- v[,-cat.i]
    }
    mss <- mess.vec(p, v)
    mess.quant <- quantile(mss)
    names(mess.quant) <- paste0("mess.", gsub("%", "p", names(mess.quant)))
  }else{
    mess.quant <- NULL
  }
  
  stats <- c(auc.test = auc.test, auc.diff = auc.diff, or.mtp = or.mtp, 
             or.10p = or.10p, other = 3, mess.quant)
  
  return(stats)
}
