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
                     args = args, auc = auc, pred = pred, 
                     nparams = nparams)