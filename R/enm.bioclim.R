################################# #
# bioclim ENMdetails object ####
################################# #

bioclim.name <- "bioclim"

bioclim.fun <- dismo::bioclim

bioclim.errors <- function(occs, envs, bg, tune.args, partitions, algorithm, 
                           partition.settings, other.settings, 
                           categoricals, doClamp, clamp.directions) {
  if(!is.null(categoricals)) {
    stop("* BIOCLIM cannot run with categorical variables. Please remove these variables before running.")
  }
}

bioclim.msgs <- function(tune.args, other.settings) {
  msg <- paste0("BIOCLIM from dismo v", packageVersion('dismo'))
  return(msg)
}

bioclim.args <- function(occs.z, bg.z, tune.tbl.i, other.settings) {
  out <- list()
  out$x <- occs.z 
  out <- c(out, other.settings$other.args)
  return(out)
}

bioclim.predict <- function(mod, envs, tune.tbl.i, other.settings) {
  # if no tails in other.args, defaults to NULL
  pred <- dismo::predict(mod, envs, tails = tune.tbl.i, na.rm = TRUE)
  return(pred)
}

bioclim.ncoefs <- function(mod) {
  # as no L1 regularization occurs, no parameters are dropped
  length(mod@min)
}

# no existing method in model object for variable importance
bioclim.varimp <- function(mod) {
  NULL
}

#' @export
enm.bioclim <- ENMdetails(name = bioclim.name, fun = bioclim.fun, errors = bioclim.errors,
                          msgs = bioclim.msgs, args = bioclim.args,
                          predict = bioclim.predict, ncoefs = bioclim.ncoefs, varimp = bioclim.varimp)
