################################# #
# bioclim ENMdetails object ####
################################# #

name <- "bioclim"

fun <- dismo::bioclim

pkgs <- c("dismo", "raster")

msgs <- function(tune.args) {
  msg <- paste0("BIOCLIM from dismo v", packageVersion('dismo'))
  return(msg)
}

args <- function(occs.vals, bg.vals, tune.tbl.i, other.args) {
  out <- list()
  out$x <- occs.vals 
  out <- c(out, other.args)
  return(out)
}

aic <- function(occs, nparam, mod.full.pred.all) NULL

eval <- function(occs.vals, bg.vals, mod, other.args, doClamp) {
  e <- dismo::evaluate(occs.vals, bg.vals, mod, tails = other.args$tails)
  return(e)
}

kstats <- function(kstats, e.test, mod, occs.train.vals, occs.test.vals, 
                   bg.train.vals, bg.test.vals, occs.train.pred, occs.test.pred, other.args) {
  
  return(kstats)
}

pred <- function(mod, envs, other.args, doClamp, pred.type) {
  # if no tails in other.args, defaults to NULL
  pred <- dismo::predict(mod, envs, tails = other.args$tails, na.rm = TRUE)
  return(pred)
}

nparams <- function(mod) {
  # as no L1 regularization occurs, no parameters are dropped
  length(mod@min)
}

#' @export
enm.bioclim <- ENMdetails(name = name, fun = fun, pkgs = pkgs, msgs = msgs, 
                     args = args, aic = aic, eval = eval, kstats = kstats, 
                     pred = pred, nparams = nparams)