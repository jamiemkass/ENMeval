################################# #
# brt ENMdetails object ####
################################# #

name <- "brt"

fun <- dismo::gbm.step

pkgs <- c("dismo", "raster", "gbm")

msgs <- function(tune.args) {
  if(!all("tree.complexity" %in% names(tune.args), "learning.rate" %in% names(tune.args), "bag.fraction" %in% names(tune.args))) {
    stop("BRT settings must include 'tree.complexity', 'learning.rate', and 'bag.fraction'.")
  }
  # construct user message with version info
  msg <- paste0("boosted regression trees (BRTs) using the gbm.step() function from gbm package v", 
                packageVersion('gbm'), " and dismo package v", packageVersion('dismo')) 
  return(msg)
}

args <- function(occs.z, bg.z, tune.i, other.settings) {
  out <- list()
  d <- rbind(occs.z, bg.z)
  p <- c(rep(1, nrow(occs.z)), rep(0, nrow(bg.z)))
  out$data <- cbind(p, d)
  out$tree.complexity <- tune.i$tree.complexity
  out$learning.rate <- tune.i$learning.rate
  out$bag.fraction <- tune.i$bag.fraction
  out$gbm.x <- 2:ncol(out$data)
  out$gbm.y <- 1
  out$silent <- TRUE
  out <- c(out, other.settings$other.args)
  return(out)
}

eval.train <- function(occs.xy, bg.xy, occs.z, bg.z, mod.full, mod.full.pred, envs, other.settings) {
  # training AUC
  e <- dismo::evaluate(occs.z, bg.z, mod.full, n.trees = length(mod.full$trees))
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

eval.validate <- function(occs.val.xy, occs.train.xy, bg.xy, occs.train.z, occs.val.z, bg.z, bg.val.z, mod.k, nk, other.settings) {
  ## validation AUC
  # calculate auc on validation data: validation occurrences are evaluated on full background, as in Radosavljevic & Anderson 2014
  # for auc.diff calculation, to perform the subtraction, it is essential that both stats are calculated over the same background
  e.train <- dismo::evaluate(occs.train.z, bg.z, mod.k, n.trees = length(mod.k$trees))
  auc.train <- e.train@auc
  # AUC validation
  # calculate AUC on validation data: if training and validation occurrences are evaluated on same background (full), auc.diff can
  # also be calculated (Radosavljevic & Anderson 2014); if validation occurrences are evaluated on partitioned background, auc.diff is NULL
  if(other.settings$validation.bg == "full") {
    e.val <- dismo::evaluate(occs.val.z, bg.z, mod.k, n.trees = length(mod.k$trees))
    auc.val <- e.val@auc
    # calculate auc diff as auc train (partition not k) minus auc validation (partition k)
    auc.diff <- auc.train - auc.val
  } else if(other.settings$validation.bg == "partition") {
    e.val <- dismo::evaluate(occs.val.z, bg.val.z, mod.k, n.trees = length(mod.k$trees))
    auc.val <- e.val@auc
    auc.diff <- NA
  }
  
  if(other.settings$abs.auc.diff == TRUE & !is.null(auc.diff)) auc.diff <- abs(auc.diff)
  
  ## omission rates
  # get model predictions for training and validation data
  occs.train.pred <- enm.brt@pred(mod.k, occs.train.z, other.settings)
  occs.val.pred <- enm.brt@pred(mod.k, occs.val.z, other.settings)
  # get minimum training presence threshold (expected no omission)
  min.train.thr <- min(occs.train.pred)
  or.mtp <- mean(occs.val.pred < min.train.thr)
  # get 10 percentile training presence threshold (expected 0.1 omission)
  pct10.train.thr <- calc.10p.trainThresh(occs.train.pred)
  or.10p <- mean(occs.val.pred < pct10.train.thr)
  
  ## validation CBI
  if(other.settings$cbi.cv == TRUE) {
    if(other.settings$validation.bg == "full") {
      mod.k.pred <- enm.brt@pred(mod.k, bg.z, other.settings)  
    }else if(other.settings$validation.bg == "partition") {
      mod.k.pred <- enm.brt@pred(mod.k, bg.val.z, other.settings)  
    }
    cbi.val <- ecospat::ecospat.boyce(mod.k.pred, occs.val.pred, PEplot = FALSE)$Spearman.cor
  }else{
    cbi.val <- NA
  }
  
  # gather all evaluation statistics for k
  out.df <- data.frame(auc.val = auc.val, auc.diff = auc.diff, cbi.val = cbi.val, or.mtp = or.mtp, or.10p = or.10p)
  
  return(out.df)
}

pred <- function(mod, envs, other.settings) {
  if(inherits(envs, "BasicRaster") == TRUE) {
    pred <- raster::predict(envs, mod, type = "response", n.trees = mod$gbm.call$best.trees, na.rm = TRUE)
  }else{
    pred <- dismo::predict(mod, envs, type = "response", n.trees = mod$gbm.call$best.trees, na.rm = TRUE)  
  }
  return(pred)
}

nparams <- function(mod) {
  # as no L1 regularization occurs, no parameters are dropped
  length(mod$var.names)
}

#' @export
enm.brt <- ENMdetails(name = name, fun = fun, pkgs = pkgs, msgs = msgs, args = args, 
                      eval.train = eval.train, eval.validate = eval.validate,
                      pred = pred, nparams = nparams)
