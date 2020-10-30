################################# #
# bioclim ENMdetails object ####
################################# #

bioclim.name <- "bioclim"

bioclim.fun.train <- dismo::bioclim

bioclim.fun.val <- bioclim.fun.train

bioclim.msgs <- function(tune.args, other.settings) {
  if(!is.null(other.settings$categoricals)) {
    stop("* BIOCLIM cannot run with categorical variables. Please remove these variables before running.")
  }
  msg <- paste0("BIOCLIM from dismo v", packageVersion('dismo'))
  return(msg)
}

bioclim.args.train <- function(occs.z, bg.z, tune.tbl.i, other.settings) {
  out <- list()
  out$x <- occs.z 
  out <- c(out, other.settings$other.args)
  return(out)
}

bioclim.args.val <- function(occs.z, bg.z, tune.tbl.i, other.settings, mod = NULL) {
  bioclim.args.train(occs.z, bg.z, tune.tbl.i, other.settings)
}

bioclim.predict <- function(mod, envs, doClamp, other.settings) {
  # if no tails in other.args, defaults to NULL
  pred <- dismo::predict(mod, envs, tails = other.settings$tails, na.rm = TRUE)
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
enm.bioclim <- ENMdetails(name = bioclim.name, fun.train = bioclim.fun.train, fun.val = bioclim.fun.val,
                          msgs = bioclim.msgs, args.train = bioclim.args.train, args.val = bioclim.args.val,
                          predict = bioclim.predict, ncoefs = bioclim.ncoefs, varimp = bioclim.varimp)
