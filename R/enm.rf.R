################################# #
# random forest ENMdetails object ####
################################# #

name <- "random forest"

fun <- randomForest::randomForest

pkgs <- c("randomForest", "dismo", "raster")

msgs <- function(tune.args) {
  if(!all("ntree" %in% names(tune.args), "mtry" %in% names(tune.args))) {
    stop("RF settings must include 'ntree' and 'mtry'.")
  }
  # construct user message with version info
  msg <- paste0("random forest using the randomForest() function from randomForest package v", 
                packageVersion('randomForest')) 
  return(msg)
}

args <- function(occs.z, bg.z, tune.i, other.settings) {
  out <- list()
  out$x <- rbind(occs.z, bg.z)
  out$y <- factor(c(rep(1, nrow(occs.z)), rep(0, nrow(bg.z))))
  out$ntree <- tune.i$ntree
  out$mtry <- tune.i$mtry
  # out$sampsize <- tune.i$sampsize
  # out$nodesize <- tune.i$nodesize
  out$importance <- TRUE
  out <- c(out, other.settings$other.args)
  return(out)
}

evaluate <- function(occs.z, bg.z, mod, other.settings) {
  dismo::evaluate(occs.z, bg.z, mod, type = "prob", index = 2)
}

train <- function(occs.xy, bg.xy, occs.z, bg.z, mod.full, mod.full.pred, envs, other.settings) {
  # training AUC
  e <- enm.rf@evaluate(occs.z, bg.z, mod.full, other.settings)
  auc.train <- e@auc
  # training CBI
  if(!is.null(envs)) {
    # if raster envs exists, use it to calculate CBI
    cbi.train.out <- ecospat::ecospat.boyce(mod.full.pred, occs.xy, PEplot = FALSE)
  }else{
    # if no raster envs exists, calculate CBI with occs and bg predictions
    d.full.pred <- dplyr::bind_rows(occs.z, bg.z)
    d.full.pred$pred <- mod.full.pred
    d.full.pred$pb <- c(rep(1, nrow(occs.z)), rep(0, nrow(bg.z)))
    cbi.train.out <- ecospat::ecospat.boyce(d.full.pred$pred, d.full.pred[d.full.pred$pb == 1, "pred"], PEplot = FALSE)
  }
  cbi.train <- cbi.train.out$Spearman.cor
  out.df <- data.frame(auc.train = auc.train, cbi.train = cbi.train)
  return(out.df)
}

validate <- function(occs.val.xy, occs.train.xy, bg.xy, occs.train.z, occs.val.z, bg.z, bg.val.z, mod.k, nk, other.settings) {
  ## validation AUC
  # calculate auc on validation data: validation occurrences are evaluated on full background, as in Radosavljevic & Anderson 2014
  # for auc.diff calculation, to perform the subtraction, it is essential that both stats are calculated over the same background
  e.train <- enm.rf@evaluate(occs.train.z, bg.z, mod.k, other.settings)
  auc.train <- e.train@auc
  # AUC validation
  # calculate AUC on validation data: if training and validation occurrences are evaluated on same background (full), auc.diff can
  # also be calculated (Radosavljevic & Anderson 2014); if validation occurrences are evaluated on partitioned background, auc.diff is NULL
  if(other.settings$validation.bg == "full") {
    e.val <- enm.rf@evaluate(occs.val.z, bg.z, mod.k, other.settings)
    auc.val <- e.val@auc
    # calculate auc diff as auc train (partition not k) minus auc validation (partition k)
    auc.diff <- auc.train - auc.val
  } else if(other.settings$validation.bg == "partition") {
    e.val <- enm.rf@evaluate(occs.val.z, bg.val.z, mod.k, other.settings)
    auc.val <- e.val@auc
    auc.diff <- NA
  }
  
  if(other.settings$abs.auc.diff == TRUE & !is.na(auc.diff)) auc.diff <- abs(auc.diff)
  
  ## omission rates
  # get model predictions for training and validation data
  occs.train.pred <- enm.rf@predict(mod.k, occs.train.z, other.settings)
  occs.val.pred <- enm.rf@predict(mod.k, occs.val.z, other.settings)
  # get minimum training presence threshold (expected no omission)
  min.train.thr <- min(occs.train.pred)
  or.mtp <- mean(occs.val.pred < min.train.thr)
  # get 10 percentile training presence threshold (expected 0.1 omission)
  pct10.train.thr <- calc.10p.trainThresh(occs.train.pred)
  or.10p <- mean(occs.val.pred < pct10.train.thr)
  
  ## validation CBI
  if(other.settings$cbi.cv == TRUE) {
    if(other.settings$validation.bg == "full") {
      mod.k.pred <- enm.rf@predict(mod.k, bg.z, other.settings)  
    }else if(other.settings$validation.bg == "partition") {
      mod.k.pred <- enm.rf@predict(mod.k, bg.val.z, other.settings)  
    }
    cbi.val <- ecospat::ecospat.boyce(mod.k.pred, occs.val.pred, PEplot = FALSE)$Spearman.cor
  }else{
    cbi.val <- NA
  }
  
  # gather all evaluation statistics for k
  out.df <- data.frame(auc.val = auc.val, auc.diff = auc.diff, cbi.val = cbi.val, or.mtp = or.mtp, or.10p = or.10p)
  
  return(out.df)
}

predict <- function(mod, envs, other.settings) {
  if(inherits(envs, "BasicRaster") == TRUE) {
    pred <- raster::predict(envs, mod, type = "prob", index = 2, na.rm = TRUE)
  }else{
    pred <- dismo::predict(mod, envs, type = "prob", index = 2, na.rm = TRUE)  
  }
  return(pred)
}

nparams <- function(mod) {
  # as no L1 regularization occurs, no parameters are dropped
  nrow(mod$importance)
}

#' @export
enm.rf <- ENMdetails(name = name, fun = fun, pkgs = pkgs, msgs = msgs, args = args, 
                      evaluate = evaluate, train = train, validate = validate,
                      predict = predict, nparams = nparams)
