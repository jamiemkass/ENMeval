################################# #
# bioclim ENMdetails object ####
################################# #

name <- "bioclim"

fun <- dismo::bioclim

msgs <- function(tune.args) {
  msg <- paste0("BIOCLIM from dismo v", packageVersion('dismo'))
  return(msg)
}

args <- function(occs.z, bg.z, tune.tbl.i, other.settings) {
  out <- list()
  out$x <- occs.z 
  out <- c(out, other.settings$other.args)
  return(out)
}

predict <- function(mod, envs, other.settings) {
  # if no tails in other.args, defaults to NULL
  pred <- dismo::predict(mod, envs, tails = other.settings$tails, na.rm = TRUE)
  return(pred)
}

nparams <- function(mod) {
  # as no L1 regularization occurs, no parameters are dropped
  length(mod@min)
}

#' @export
enm.bioclim <- ENMdetails(name = name, fun = fun, msgs = msgs, args = args, 
                          predict = predict, nparams = nparams)
