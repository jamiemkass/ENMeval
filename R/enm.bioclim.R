################################# #
# bioclim ENMdetails object ####
################################# #

fun <- dismo::bioclim

pkgs <- c("dismo", "raster")

msgs <- function(tune.args) {
  msg <- paste("BIOCLIM from dismo v.", packageVersion('dismo'))
  return(msg)
}

args <- function(occs.vals, bg.vals, tune.tbl.i, other.args) {
  out <- list()
  out$x <- occs.vals 
  out <- c(out, other.args)
  return(out)
}

auc <- function(occs.vals, bg.vals, mod, other.args, doClamp) {
  e <- dismo::evaluate(occs.vals, bg.vals, mod, tails = other.args$tails)@auc
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
  # if no tails in other.args, defaults to NULL
  pred <- dismo::predict(mod, envs, tails = other.args$tails, na.rm = TRUE)
  return(pred)
}

nparams <- function(mod) {
  # as no L1 regularization occurs, no parameters are dropped
  length(mod@min)
}

enm.bioclim <- ENMdetails(fun = fun, pkgs = pkgs, msgs = msgs, 
                     args = args, auc = auc, kstats = kstats, 
                     pred = pred, nparams = nparams)