#############################################
#########	TUNING FUNCTION	#############
#############################################
# THIS FUNCTION DOES SPATIALLY-INDEPENDENT EVALUATIONS
# INPUT ARGUMENTS COME FROM WRAPPER FUNCTION

tuning <- function (occs, envs, bg, occs.grp, bg.grp, method, algorithm, args,
                    args.lab, categoricals, aggregation.factor, kfolds, bin.output,
                    clamp, alg, rasterPreds, parallel, numCores, progbar, updateProgress,
                    userArgs) {

  # extract predictor variable values at coordinates for occs and bg
  occs.vals <- as.data.frame(extract(envs, occs))
  bg.vals <- as.data.frame(extract(envs, bg))

  # remove rows with NA for any predictor variable
  # also redefine occs and bg without NA rows
  occs.vals.na <- sum(rowSums(is.na(occs.vals)) > 0)
  bg.vals.na <- sum(rowSums(is.na(bg)) > 0)
  if(occs.vals.na > 0) {
    i <- !apply(occs.vals, 1, anyNA)
    occs <- occs[i,]
    occs.grp <- occs.grp[i]
    occs.vals <- occs.vals[i,]
    message(paste("There were", occs.vals.na, "occurrence records with NA for at least
                  one predictor variable. Removed these from analysis,
                  resulting in", nrow(occs.vals), "occurrence records."))
  }
  if(bg.vals.na > 0) {
    i <- !apply(bg.vals, 1, anyNA)
    occs <- occs[i,]
    bg.grp <- bg.grp[i]
    bg.vals <- bg.vals[i,]
    message(paste("There were", bg.vals.na, "background records with NA for at least
                  one predictor variable. Removed these from analysis,
                  resulting in", nrow(bg.vals), "background records."))
  }

  if (!is.null(categoricals)) {
    for (i in 1:length(categoricals)) {
      pres[, categoricals[i]] <- as.factor(pres[, categoricals[i]])
      bg[, categoricals[i]] <- as.factor(bg[, categoricals[i]])
    }
  }

  # assign partition groups based on choice of method
  if ("checkerboard1" %in% method) {
    method <- c(method = "checkerboard1", aggregation.factor = aggregation.factor)
    group.data <- get.checkerboard1(occs, envs, bg, aggregation.factor)
  }
  if ("checkerboard2" %in% method) {
    method <- c(method = "checkerboard2", aggregation.factor = aggregation.factor)
    group.data <- get.checkerboard2(occs, envs, bg, aggregation.factor)
  }
  if ("block" %in% method)
    group.data <- get.block(occs, bg)
  if ("jackknife" %in% method)
    group.data <- get.jackknife(occs, bg)
  if ("randomkfold" %in% method) {
    method <- c(method = "randomkfold", number.folds = kfolds)
    group.data <- get.randomkfold(occs, bg, kfolds)
  }
  if ("user" %in% method) {
    method <- c(method= "user", number.folds = length(unique(occs.grp)))
    group.data <- get.user(occs.grp, bg.grp)
  }


  # define number of groups (the value of "k")
  nk <- length(unique(group.data$occs.grp))

  # differential behavior for parallel and default
  if (parallel == TRUE) {
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
                       modelTune.maxnet(pres, bg, envs, nk, group.data, args[[i]],
                                        rasterPreds, clamp)
                     }
    } else if (algorithm == 'maxent.jar') {
      out <- foreach(i = seq_len(length(args)),
                     .packages = c("dismo", "raster", "ENMeval", "rJava")) %dopar% {
                       modelTune.maxentJar(pres, bg, envs, nk, group.data, args[[i]],
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
        out[[i]] <- modelTune.maxnet(pres, bg, envs, nk, group.data, args[[i]],
                                     rasterPreds, clamp)
      } else if (algorithm == 'maxent.jar') {
        out[[i]] <- modelTune.maxentJar(pres, bg, envs, nk, group.data, args[[i]],
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
      full.AUC[i] <- dismo::evaluate(pres, bg, full.mods[[i]])@auc
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

  results <- ENMevaluation(algorithm = alg, results = res, predictions = predictive.maps,
                           models = full.mods, partition.method = method,
                           occs.pts = occs, occs.grp = group.data[[1]],
                           bg.pts = bg, bg.grp = group.data[[2]])
  return(results)
}
