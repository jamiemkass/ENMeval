################################# #
# maxnet ENMdetails object ####
################################# #

fun <- maxnet::maxnet

pkgs <- c("dismo", "raster", "maxnet")

msgs <- function(tune.args) {
  if(!("rm" %in% names(tune.args)) | !("fc" %in% names(tune.args))) {
    stop("For Maxent, please specify both 'rm' and 'fc' settings. See ?tune.args for help.")
  }else{
    if(!is.numeric(tune.args[["rm"]])) {
      stop("Please input numeric values for 'rm' settings for Maxent.")
    }
    all.fc <- unlist(sapply(1:5, function(x) apply(combn(c("L","Q","H","P","T"), x), 2, function(y) paste(y, collapse = ""))))
    if(any(!tune.args[["fc"]] %in% all.fc)) {
      stop("Please input accepted values for 'fc' settings for Maxent.")
    }
    msg <- paste("maxnet from maxnet package v.", packageVersion('maxnet'))
    return(msg)
  }
}

args <- function(occs.vals, bg.vals, tune.tbl.i, other.args) {
  out <- list()
  out$data <- rbind(occs.vals, bg.vals)
  out$p <- c(rep(1, nrow(occs.vals)), rep(0, nrow(bg.vals)))
  out$f <- maxnet::maxnet.formula(out$p, out$data, classes = tolower(tune.tbl.i$fc))
  out$regmult <- tune.tbl.i$rm
  out <- c(out, other.args)
  return(out)
}

auc <- function(occs.vals, bg.vals, mod, other.args, doClamp) {
  e <- dismo::evaluate(occs.vals, bg.vals, mod, type = 'exponential', clamp = doClamp)@auc
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

pred <- function(mod, envs, other.args, doClamp, type) {
  pred <- maxnet.predictRaster(mod, envs, doClamp, type = 'exponential', other.args)
  return(pred)
}

nparams <- function(mod) {
  length(mod$betas)
}

enm.maxnet <- ENMdetails(fun = fun, pkgs = pkgs, msgs = msgs, 
                        args = args, auc = auc, kstats = kstats, 
                        pred = pred, nparams = nparams)