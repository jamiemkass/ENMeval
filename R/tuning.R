#############################################
#########	TUNING FUNCTION	#############
#############################################
# THIS FUNCTION DOES SPATIALLY-INDEPENDENT EVALUATIONS
# INPUT ARGUMENTS COME FROM WRAPPER FUNCTION

tuning <- function (occ, env, bg.coords, occ.grp, bg.grp, method, maxent.args,
                    args.lab, categoricals, aggregation.factor, kfolds, bin.output,
                    clamp, alg, rasterPreds, java, parallel, numCores, progbar, updateProgress,
                    userArgs) {

  noccs <- nrow(occ)

  # extract predictor variable values at coordinates for occs and bg
  pres <- as.data.frame(extract(env, occ))
  bg <- as.data.frame(extract(env, bg.coords))
  
  # if rows have NA for any predictor, toss them out
  # also redefine occ and bg without NA rows
  numNA.pres <- sum(rowSums(is.na(pres)))
  numNA.bg <- sum(rowSums(is.na(bg)))
  if (numNA.pres > 0) {
    message(paste("There are", numNA.pres, "occurrence records with NA for at least 
                  one predictor variable. Removing these records from analysis,
                  resulting in", nrow(pres) - numNA.pres, "records..."))
    occ <- occ[-which(is.na(pres)),]
    pres <- na.omit(pres)
  }
  if (numNA.bg > 0) {
    message(paste("There are", numNA.bg, "background records with NA for at least one predictor variable. 
                  Removing these records from analysis, resulting in", nrow(bg) - numNA.bg, "records..."))
    bg.coords <- bg.coords[-which(is.na(bg)),]
    bg <- na.omit(bg)
  }
  
  if (!is.null(categoricals)) {
    for (i in 1:length(categoricals)) {
      pres[, categoricals[i]] <- as.factor(pres[, categoricals[i]])
      bg[, categoricals[i]] <- as.factor(bg[, categoricals[i]])
    }
  }
  
  # assign partition groups based on choice of method
  if (method == "checkerboard1")
    group.data <- get.checkerboard1(occ, env, bg.coords, aggregation.factor)
  if (method == "checkerboard2")
    group.data <- get.checkerboard2(occ, env, bg.coords, aggregation.factor)
  if (method == "block")
    group.data <- get.block(occ, bg.coords)
  if (method == "jackknife")
    group.data <- get.jackknife(occ, bg.coords)
  if (method == "randomkfold")
    group.data <- get.randomkfold(occ, bg.coords, kfolds)
  if (method == "user")
    group.data <- get.user(occ.grp, bg.grp)
  # define number of groups (the value of "k")
  nk <- length(unique(group.data$occ.grp))

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
    #cat("Running in parallel...\n")
    out <- foreach(i = seq_len(length(maxent.args)), .packages = c("dismo", "raster", "ENMeval")) %dopar% {
      modelTune(pres, bg, env, nk, group.data, progbar, maxent.args, 
                userArgs, rasterPreds, clamp, java, updateProgress)
    }
    stopCluster(c1)
  } else {
      out <- modelTune(pres, bg, env, nk, group.data, progbar, maxent.args, 
                       userArgs, rasterPreds, clamp, java, updateProgress)
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
  names(AUC.DIFF) <- paste("AUC.DIFF_bin", 1:nk, sep = ".")
  Mean.AUC.DIFF <- rowMeans(AUC.DIFF)
  Var.AUC.DIFF <- corrected.var(AUC.DIFF, noccs)
  names(AUC.TEST) <- paste("AUC_bin", 1:nk, sep = ".")
  Mean.AUC <- rowMeans(AUC.TEST)
  Var.AUC <- corrected.var(AUC.TEST, noccs)
  names(OR10) <- paste("OR10_bin", 1:nk, sep = ".")
  Mean.OR10 <- rowMeans(OR10)
  Var.OR10 <- apply(OR10, 1, var)
  names(ORmin) <- paste("ORmin_bin", 1:nk, sep = ".")
  Mean.ORmin <- rowMeans(ORmin)
  Var.ORmin <- apply(ORmin, 1, var)

  # get training AUCs for each model
  full.AUC <- double()
  
  for (i in 1:length(full.mods)) {
    if (java == TRUE) {
      full.AUC[i] <- full.mods[[i]]@results[5]  
    } else {
      full.AUC[i] <- evaluate(pres, bg, full.mods[[i]])@auc
    }
  }
  
  # get total number of parameters
  nparam <- numeric()
  for (i in 1:length(full.mods)) {
    if (java ==TRUE) {
      nparam[i] <- get.params(full.mods[[i]])  
    } else {
      nparam[i] <- length(full.mods[[i]]$betas)
    }
  }
    
#  if (rasterPreds==TRUE) { # this should now work even if rasterPreds==F
    aicc <- calc.aicc(nparam, occ, predictive.maps)
#  } else {
#    aicc <- rep(NaN, length(full.AUC))
#  }

  features <- args.lab[[1]]
  rm <- args.lab[[2]]
  settings <- paste(args.lab[[1]], args.lab[[2]], sep = "_")

  res <- data.frame(settings, features, rm, full.AUC, Mean.AUC,
                    Var.AUC, Mean.AUC.DIFF, Var.AUC.DIFF, Mean.OR10, Var.OR10,
                    Mean.ORmin, Var.ORmin, aicc)
  if (bin.output == TRUE) {
    res <- as.data.frame(cbind(res, AUC.TEST, AUC.DIFF, OR10, ORmin))
  }

  if (rasterPreds==TRUE) {
    names(predictive.maps) <- settings
  }

  results <- ENMevaluation(algorithm = alg, results = res, predictions = predictive.maps,
                           models = full.mods, partition.method = method, 
                           occ.pts = occ, occ.grp = group.data[[1]],
                           bg.pts = bg.coords, bg.grp = group.data[[2]])
  return(results)
}
