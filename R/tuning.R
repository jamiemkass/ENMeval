#############################################
#########	TUNING FUNCTION	#############
#############################################
# THIS FUNCTION DOES SPATIALLY-INDEPENDENT EVALUATIONS
# INPUT ARGUMENTS COME FROM WRAPPER FUNCTION

tuning <- function (occs, envs, bg, folds, algorithm, args,
                    args.lab, categoricals, aggregation.factor, bin.output,
                    clamp, rasterPreds, parallel, numCores, progbar, updateProgress,
                    userArgs) {

  # unpack occs.folds and bg.folds
  occs.folds <- folds$occs.folds
  bg.folds <- folds$bg.folds
  
  # extract predictor variable values at coordinates for occs and bg
  occs.vals <- as.data.frame(extract(envs, occs))
  bg.vals <- as.data.frame(extract(envs, bg))

  # remove rows from occs, occs.vals, occs.folds with NA for any predictor variable
  occs.vals.na <- sum(rowSums(is.na(occs.vals)) > 0)
  bg.vals.na <- sum(rowSums(is.na(bg)) > 0)
  if(occs.vals.na > 0) {
    i <- !apply(occs.vals, 1, anyNA)
    occs <- occs[i,]
    occs.folds <- occs.folds[i]
    occs.vals <- occs.vals[i,]
    message(paste("There were", occs.vals.na, "occurrence records with NA for at least
                  one predictor variable. Removed these from analysis,
                  resulting in", nrow(occs.vals), "occurrence records."))
  }
  # do the same for associated bg variables
  if(bg.vals.na > 0) {
    i <- !apply(bg.vals, 1, anyNA)
    occs <- occs[i,]
    bg.folds <- bg.folds[i]
    bg.vals <- bg.vals[i,]
    message(paste("There were", bg.vals.na, "background records with NA for at least
                  one predictor variable. Removed these from analysis,
                  resulting in", nrow(bg.vals), "background records."))
  }

  # convert fields for categorical data to factor class
  if(!is.null(categoricals)) {
    for (i in 1:length(categoricals)) {
      occs.vals[, categoricals[i]] <- as.factor(occs.vals[, categoricals[i]])
      bg.vals[, categoricals[i]] <- as.factor(bg.vals[, categoricals[i]])
    }
  }

  # define number of folds (the value of "k")
  nk <- length(unique(occs.folds))

  # differential behavior for parallel and default
  if(parallel == TRUE) {
    # set up parallel computing
    allCores <- detectCores()
    if (is.null(numCores)) {
      numCores <- allCores
    }
    c1 <- makeCluster(numCores)
    registerDoParallel(c1)
    numCoresUsed <- getDoParWorkers()
    message(paste("Of", allCores, "total cores using", numCoresUsed))
    #cat(paste("Of", allCores, "total cores using", numCoresUsed, "\n"))

    # log file to record status of parallel loops
    message("Running in parallel...")
    if (algorithm == 'maxnet') {
      # load maxnet package because it is not loaded automatically when embedded in foreach
      # if not loaded here, dismo::evaluate fails later
      # library(maxnet)
      out <- foreach(i = seq_len(length(args)),
                     .packages = c("dismo", "raster", "ENMeval", "maxnet")) %dopar% {
                       modelTune.maxnet(occs.vals, bg.vals, envs, nk, folds, args[[i]],
                                        rasterPreds, clamp)
                     }
    } else if (algorithm == 'maxent.jar') {
      out <- foreach(i = seq_len(length(args)),
                     .packages = c("dismo", "raster", "ENMeval", "rJava")) %dopar% {
                       modelTune.maxentJar(occs.vals, bg.vals, envs, nk, folds, args[[i]],
                                           userArgs, rasterPreds, clamp)
                     }
    }
    stopCluster(c1)
  } else {
    out <- list()
    if (progbar == TRUE & !is.function(updateProgress)) {
      pb <- txtProgressBar(0, length(args), style = 3)
    }
    for (i in 1:length(args)) {
      # set up the console progress bar (progbar),
      # or the shiny progress bar (updateProgress)
      if (length(args) > 1) {
        if (is.function(updateProgress)) {
          text <- paste0('Running ', args.lab[[1]][i], args.lab[[2]][i], '...')
          updateProgress(detail = text)
        } else if (progbar == TRUE) {
          setTxtProgressBar(pb, i)
        }
      }
      if (algorithm == 'maxnet') {
        out[[i]] <- modelTune.maxnet(occs.vals, bg.vals, envs, nk, folds, args[[i]],
                                     rasterPreds, clamp)
      } else if (algorithm == 'maxent.jar') {
        out[[i]] <- modelTune.maxentJar(occs.vals, bg.vals, envs, nk, folds, args[[i]],
                                        userArgs, rasterPreds, clamp)
      }
    }
    if (progbar==TRUE) close(pb)
  }

  # gather all full models into list
  full.mods <- lapply(out, function(x) x[[1]])
  # gather all statistics into a data frame
  statsTbl <- as.data.frame(t(sapply(out, function(x) x[[2]])))
  if (rasterPreds) {
    predictive.maps <- stack(sapply(out, function(x) x[[3]]))
  } else {
    predictive.maps <- stack()
  }

  AUC.DIFF <- statsTbl[,1:nk]
  AUC.TEST <- statsTbl[,(nk+1):(2*nk)]
  OR10 <- statsTbl[,((2*nk)+1):(3*nk)]
  ORmin <- statsTbl[,((3*nk)+1):(4*nk)]
  # rename column fields
  names(AUC.DIFF) <- paste("diff.AUC_bin", 1:nk, sep = ".")
  Mean.AUC.DIFF <- rowMeans(AUC.DIFF)
  Var.AUC.DIFF <- corrected.var(AUC.DIFF, nk)
  names(AUC.TEST) <- paste("AUC_bin", 1:nk, sep = ".")
  Mean.AUC <- rowMeans(AUC.TEST)
  Var.AUC <- corrected.var(AUC.TEST, nk)
  names(OR10) <- paste("test.or10pct_bin", 1:nk, sep = ".")
  Mean.OR10 <- rowMeans(OR10)
  Var.OR10 <- apply(OR10, 1, var)
  names(ORmin) <- paste("test.orMTP_bin", 1:nk, sep = ".")
  Mean.ORmin <- rowMeans(ORmin)
  Var.ORmin <- apply(ORmin, 1, var)

  # get training AUCs for each model
  full.AUC <- double()

  for (i in 1:length(full.mods)) {
    if (algorithm == 'maxnet') {
      full.AUC[i] <- dismo::evaluate(occs.vals, bg.vals, full.mods[[i]])@auc
    } else if (algorithm == 'maxent.jar') {
      full.AUC[i] <- full.mods[[i]]@results[5]
    }
  }

  # get total number of parameters
  nparam <- numeric()
  for (i in 1:length(full.mods)) {
    if (algorithm == 'maxnet') {
      nparam[i] <- length(full.mods[[i]]$betas)
    } else if (algorithm == 'maxent.jar') {
      nparam[i] <- get.params(full.mods[[i]])
    }
  }

#  if (rasterPreds==TRUE) { # this should now work even if rasterPreds==F
    aicc <- calc.aicc(nparam, occs, predictive.maps)
#  } else {
#    aicc <- rep(NaN, length(full.AUC))
#  }

  features <- args.lab[[1]]
  rm <- args.lab[[2]]
  settings <- paste(args.lab[[1]], args.lab[[2]], sep = "_")

  res <- data.frame(settings, features, rm, train.AUC = full.AUC,
                    avg.test.AUC = Mean.AUC, var.test.AUC = Var.AUC,
                    avg.diff.AUC = Mean.AUC.DIFF, var.diff.AUC = Var.AUC.DIFF,
                    avg.test.orMTP = Mean.ORmin, var.test.orMTP = Var.ORmin,
                    avg.test.or10pct = Mean.OR10, var.test.or10pct = Var.OR10, aicc)
  if (bin.output == TRUE) {
    res <- as.data.frame(cbind(res, AUC.TEST, AUC.DIFF, OR10, ORmin))
  }

  if (rasterPreds==TRUE) {
    names(predictive.maps) <- settings
  }

  results <- ENMevaluation(algorithm = algorithm, results = res, predictions = predictive.maps,
                           models = full.mods, partition.method = method,
                           occs.pts = occs, occs.folds = folds[[1]],
                           bg.pts = bg, bg.folds = folds[[2]])
  return(results)
}
