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

args <- function(occs.vals, bg.vals, tune.i, other.settings) {
  out <- list()
  d <- rbind(occs.vals, bg.vals)
  p <- c(rep(1, nrow(occs.vals)), rep(0, nrow(bg.vals)))
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

eval.train <- function(occs.xy, bg.xy, occs.vals, bg.vals, mod.full, mod.full.pred, envs, other.settings) {
  # training AUC
  e <- dismo::evaluate(occs.vals, bg.vals, mod.full, n.trees = length(mod.full$trees))
  auc.train <- e@auc
  # training CBI
  if(!is.null(envs)) {
    # if raster envs exists, use it to calculate CBI
    cbi.train.out <- ecospat::ecospat.boyce(mod.full.pred, occs.xy, PEplot = FALSE)
  }else{
    # if no raster envs exists, calculate CBI with occs and bg predictions
    d.full.pred <- dplyr::bind_rows(occs.vals, bg.vals)
    d.full.pred$pred <- mod.full.pred
    d.full.pred$pb <- c(rep(1, nrow(occs.vals)), rep(0, nrow(bg.vals)))
    cbi.train.out <- ecospat::ecospat.boyce(d.full.pred$pred, d.full.pred[d.full.pred$pb == 1, "pred"], PEplot = FALSE)
  }
  cbi.train <- cbi.train.out$Spearman.cor
  out.df <- data.frame(auc.train = auc.train, cbi.train = cbi.train)
  return(out.df)
}

eval.test <- function(occs.test.xy, occs.train.xy, bg.xy, occs.train.vals, occs.test.vals, bg.vals, mod.k, nk, envs, other.settings) {
  ## testing AUC
  # calculate auc on testing data: test occurrences are evaluated on full background, as in Radosavljevic & Anderson 2014
  # for auc.diff calculation, do perform the subtraction, it is essential that both stats are calculated over the same background
  e.train <- dismo::evaluate(occs.train.vals, bg.vals, mod.k, n.trees = length(mod.k$trees))
  e.test <- dismo::evaluate(occs.test.vals, bg.vals, mod.k, n.trees = length(mod.k$trees))
  auc.train <- e.train@auc
  auc.test <- e.test@auc
  # calculate auc diff as auc train (partition not k) minus auc test (partition k)
  auc.diff <- auc.train - auc.test
  if(other.settings$abs.auc.diff == TRUE) auc.diff <- abs(auc.diff)
  
  ## omission rates
  # get model predictions for training and testing data
  occs.train.pred <- enm.brt@pred(mod.k, occs.train.vals, other.settings)
  occs.test.pred <- enm.brt@pred(mod.k, occs.test.vals, other.settings)
  # get minimum training presence threshold (expected no omission)
  min.train.thr <- min(occs.train.pred)
  or.mtp <- mean(occs.test.pred < min.train.thr)
  # get 10 percentile training presence threshold (expected 0.1 omission)
  pct10.train.thr <- calc.10p.trainThresh(occs.train.pred)
  or.10p <- mean(occs.test.pred < pct10.train.thr)
  
  ## testing CBI
  if(other.settings$cbi.cv == TRUE) {
    if(!is.null(envs)) {
      # predict to raster
      mod.k.pred <- enm.brt@pred(mod.k, envs, other.settings)
    }else{
      # use full background to approximate full model prediction
      mod.k.pred <- enm.brt@pred(mod.k, bg.vals, other.settings)
    }
    cbi.test <- ecospat::ecospat.boyce(mod.k.pred, occs.test.pred, PEplot = FALSE)
  }else{
    cbi.test <- NULL
  }
  
  # gather all evaluation statistics for k
  out.df <- data.frame(auc.test = auc.test, auc.diff = auc.diff, or.mtp = or.mtp, or.10p = or.10p)
  if(!is.null(cbi.test)) out.df <- cbind(out.df, cbi.test = cbi.test$Spearman.cor)
  
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
                      eval.train = eval.train, eval.test = eval.test,
                      pred = pred, nparams = nparams)