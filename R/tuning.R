#############################################
#########	TUNING FUNCTION	#############
#############################################
# THIS FUNCTION DOES SPATIALLY-INDEPENDENT EVALUATIONS
# INPUT ARGUMENTS COME FROM WRAPPER FUNCTION

tuning <- function (occ, env, bg.coords, occ.grp, bg.grp, method, maxent.args, 
                    args.lab, categoricals, aggregation.factor, kfolds, bin.output, 
                    rasterPreds, updateProgress) {
  
  noccs <- nrow(occ)
  if (method == "checkerboard1") 
    group.data <- get.checkerboard1(occ, env, bg.coords, 
                                    aggregation.factor)
  if (method == "checkerboard2") 
    group.data <- get.checkerboard2(occ, env, bg.coords, 
                                    aggregation.factor)
  if (method == "block") 
    group.data <- get.block(occ, bg.coords)
  if (method == "jackknife") 
    group.data <- get.jackknife(occ, bg.coords)
  if (method == "randomkfold") 
    group.data <- get.randomkfold(occ, bg.coords, kfolds)
  if (method == "user") 
    group.data <- get.user(occ.grp, bg.grp)
  nk <- length(unique(group.data$occ.grp))
  pres <- as.data.frame(extract(env, occ))
  bg <- as.data.frame(extract(env, bg.coords))
  if (!is.null(categoricals)) {
    pres[, categoricals] <- as.factor(pres[, categoricals])
    bg[, categoricals] <- as.factor(bg[, categoricals])
  }
  if (length(maxent.args) > 1 & !is.function(updateProgress)) {
    pb <- txtProgressBar(0, length(maxent.args), style = 3) 
  }
  full.mod <- list()
  AUC.TEST <- data.frame()
  AUC.DIFF <- data.frame()
  OR10 <- data.frame()
  ORmin <- data.frame()
  predictive.maps <- stack()
  nparm <- vector()
  full.AUC <- vector()
  for (a in 1:length(maxent.args)) {
    if (length(maxent.args) > 1) {
      if (is.function(updateProgress)) {
        text <- paste0('Running ', args.lab[[1]][a], args.lab[[2]][a], '...')
        updateProgress(detail = text)
      } else {
        setTxtProgressBar(pb, a) 
      }
    }
      
    x <- rbind(pres, bg)
    p <- c(rep(1, nrow(pres)), rep(0, nrow(bg)))
    tmpfolder <- tempfile()
    full.mod[a] <- maxent(x, p, args = maxent.args[[a]], 
                          factors = categoricals, path = tmpfolder)
    pred.args <- c("outputformat=raw", "doclamp=true")
    if (rasterPreds==TRUE) {
      predictive.maps <- stack(predictive.maps, predict(full.mod[[a]], env, args = pred.args))  
    }
    full.AUC[a] <- full.mod[[a]]@results[5]
    for (k in 1:nk) {
      train.val <- pres[group.data$occ.grp != k, ]
      test.val <- pres[group.data$occ.grp == k, ]
      bg.val <- bg[group.data$bg.grp != k, ]
      x <- rbind(train.val, bg.val)
      p <- c(rep(1, nrow(train.val)), rep(0, nrow(bg.val)))
      mod <- maxent(x, p, args = maxent.args[[a]], factors = categoricals, 
                    path = tmpfolder)
      AUC.TEST[a, k] <- evaluate(test.val, bg, mod)@auc
      AUC.DIFF[a, k] <- max(0, evaluate(train.val, bg, 
                                        mod)@auc - AUC.TEST[a, k])
      p.train <- predict(mod, train.val, args = pred.args)
      p.test <- predict(mod, test.val, args = pred.args)
      if (nrow(train.val) < 10) {
        n90 <- floor(nrow(train.val) * 0.9)
      }
      else {
        n90 <- ceiling(nrow(train.val) * 0.9)
      }
      train.thr.10 <- rev(sort(p.train))[n90]
      OR10[a, k] <- mean(p.test < train.thr.10)
      train.thr.min <- min(p.train)
      ORmin[a, k] <- mean(p.test < train.thr.min)
    }
    unlink(tmpfolder, recursive = TRUE)
  }
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
  for (i in 1:length(full.mod)) nparm[i] <- get.params(full.mod[[i]])
  if (rasterPreds==TRUE) {
    aicc <- calc.aicc(nparm, occ, predictive.maps)
  } else {
    aicc <- rep(NaN, length(full.AUC))
  }
  features <- args.lab[[1]]
  rm <- args.lab[[2]]
  settings <- paste(args.lab[[1]], args.lab[[2]], sep = "_")
  res <- data.frame(settings, features, rm, full.AUC, Mean.AUC, 
                    Var.AUC, Mean.AUC.DIFF, Var.AUC.DIFF, Mean.OR10, Var.OR10, 
                    Mean.ORmin, Var.ORmin, aicc)
  if (bin.output == TRUE) {
    res <- as.data.frame(cbind(res, AUC.TEST, AUC.DIFF, OR10, 
                               ORmin))
  }
  if (rasterPreds==TRUE) {
    names(predictive.maps) <- settings
  }
  results <- ENMevaluation(results = res, predictions = predictive.maps, models = full.mod,
                           partition.method = method, occ.pts = occ, occ.grp = group.data[[1]], 
                           bg.pts = bg.coords, bg.grp = group.data[[2]])
  return(results)
}