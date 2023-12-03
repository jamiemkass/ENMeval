################################# #
# bioclim ENMdetails object ####
################################# #

bioclim.name <- "bioclim"

bioclim.fun <- predicts::envelope

bioclim.errors <- function(occs, envs, bg, tune.args, partitions, algorithm, 
                           partition.settings, other.settings, 
                           categoricals, doClamp, clamp.directions) {
  if(!is.null(categoricals)) {
    stop("* BIOCLIM cannot run with categorical variables. Please remove these variables before running.")
  }
}

bioclim.msgs <- function(tune.args, other.settings) {
  msg <- paste0("BIOCLIM from predicts v", packageVersion('predicts'))
  return(msg)
}

bioclim.args <- function(occs.z, bg.z, tune.tbl.i, other.settings) {
  out <- list()
  out$x <- occs.z 
  out <- c(out, other.settings$other.args)
  return(out)
}

# NOTE: clamping is not needed for BIOCLIM predictions because predictions
# are always clamped for this algorithm
bioclim.predict <- function(mod, envs, other.settings) {
  require(predicts)
  # if no tails in other.args, defaults to NULL
  pred <- predict(mod, envs, tails = other.settings$tails)
  return(pred)
}

bioclim.ncoefs <- function(mod) {
  # as no L1 regularization occurs, no parameters are dropped
  length(mod@min)
}

# no existing method in model object for variable importance
bioclim.variable.importance <- function(mod) {
  NULL
}

#' @title ENMdetails bioclim
#' @description This is the ENMdetails implementation for the BIOCLIM climate envelope model, implemented by predicts.
#' @export
enm.bioclim <- ENMdetails(name = bioclim.name, fun = bioclim.fun, errors = bioclim.errors,
                          msgs = bioclim.msgs, args = bioclim.args,
                          predict = bioclim.predict, ncoefs = bioclim.ncoefs, variable.importance = bioclim.variable.importance)
