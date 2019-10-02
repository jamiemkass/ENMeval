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
  msg <- paste("Boosted regression trees (BRTs) using the gbm.step() function from gbm package v.", 
               packageVersion('gbm'), "and dismo package v.", packageVersion('dismo')) 
  return(msg)
}

args <- function(occs.vals, bg.vals, tune.tbl.i, other.args) {
  out <- list()
  d <- rbind(occs.vals, bg.vals)
  p <- c(rep(1, nrow(occs.vals)), rep(0, nrow(bg.vals)))
  out$data <- cbind(p, d)
  out$tree.complexity <- tune.tbl.i$tree.complexity
  out$learning.rate <- tune.tbl.i$learning.rate
  out$bag.fraction <- tune.tbl.i$bag.fraction
  out$gbm.x <- 2:ncol(out$data)
  out$gbm.y <- 1
  out$silent <- TRUE
  out <- c(out, other.args)
  return(out)
}

aic <- function(occs, nparam, mod.full.pred.all) NULL

auc <- function(occs.vals, bg.vals, mod, other.args, doClamp) {
  e <- dismo::evaluate(occs.vals, bg.vals, mod, n.trees = length(mod$trees))@auc
  return(e)
}

kstats <- function(occs.train, bg.train, occs.test, bg.test, categoricals,
                   auc.train, mod, other.args, doClamp, abs.auc.diff) {
  # calculate auc on testing data
  auc.test <- auc(occs.test, bg.train, mod, other.args, doClamp)
  # calculate auc diff
  auc.diff <- auc.train - auc.test
  if(abs.auc.diff == TRUE) auc.diff <- abs(auc.diff)
  # get model predictions for training and testing data
  # these predictions are used only for calculating omission rate, and
  # thus should not need any specific parameter changes for maxent/maxnet
  pred.train <- pred(mod, occs.train, other.args, doClamp)
  pred.test <- pred(mod, occs.test, other.args, doClamp)
  # get minimum training presence threshold (expected no omission)
  min.train.thr <- min(pred.train)
  or.mtp <- mean(pred.test < min.train.thr)
  # get 10 percentile training presence threshold (expected 0.1 omission)
  pct10.train.thr <- calc.10p.trainThresh(occs.train, pred.train)
  or.10p <- mean(pred.test < pct10.train.thr)
  
  # calculate MESS values if bg.test values are given
  if(!is.null(bg.test) & ncol(bg.test) > 1) {
    mess.quant <- calc.mess.kstats(occs.train, bg.train, occs.test, bg.test, categoricals)
  }else{
    mess.quant <- NULL
  }
  
  stats <- c(auc.test = auc.test, auc.diff = auc.diff, or.mtp = or.mtp, 
             or.10p = or.10p, mess.quant)
  
  return(stats)
}

pred <- function(mod, envs, other.args, doClamp) {
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
enm.brt <- ENMdetails(name = name, fun = fun, pkgs = pkgs, msgs = msgs, 
                      args = args, aic = aic, auc = auc, kstats = kstats,
                      pred = pred, nparams = nparams)