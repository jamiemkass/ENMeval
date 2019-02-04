#############################################
#########	TUNING FUNCTION	#############
#############################################
# THIS FUNCTION DOES SPATIALLY-INDEPENDENT EVALUATIONS
# INPUT ARGUMENTS COME FROM WRAPPER FUNCTION

tuning <- function (occs, envs, bg, folds, algorithm, args,
                    args.lab, categoricals, aggregation.factor, bin.output,
                    clamp, rasterPreds, parallel, numCores, progbar, updateProgress,
                    userArgs) {

  

  
    
  

  

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
