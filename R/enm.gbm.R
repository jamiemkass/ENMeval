################################# #
# GBM ENMdetails object ####
################################# #

gbm.name <- "GBM"

gbm.fun <- gbm::gbm

gbm.errors <- function(occs, envs, bg, tune.args, partitions, algorithm, 
                          partition.settings, other.settings, 
                          categoricals, doClamp, clamp.directions) {
  if(!("interaction.depth" %in% names(tune.args))) {
    stop("GBM settings must include 'interaction.depth' (highest level of variable interactions allowed for each tree). See ?tune.args for details.")
  }
  if(any(tune.args$interaction.depth <= 0)) {
    stop("Please input positive integer values for 'interaction.depth' settings for gbm.")
  }
}

gbm.msgs <- function(tune.args, other.settings) {
  msg <- paste0("randomForest from randomForest package v", packageVersion('randomForest'))
  return(msg)
}

gbm.args <- function(occs.z, bg.z, tune.tbl.i, other.settings) {
  out <- list()
  out$formula <- formula(p ~.)
  out$data <- rbind(occs.z, bg.z)
  p <- c(rep(1, nrow(occs.z)), rep(0, nrow(bg.z)))
  out$data <- cbind(p, out$data)
  out$distribution <- "bernoulli"
  out$verbose <- FALSE
  # set to algorithm default if none specified
  out$n.trees <- ifelse(is.null(tune.tbl.i$n.trees), 100, tune.tbl.i$n.trees)
  # set to algorithm default if none specified
  out$bag.fraction <- ifelse(is.null(tune.tbl.i$bag.fraction), 0.5, tune.tbl.i$bag.fraction)
  out$interaction.depth <- tune.tbl.i$interaction.depth
  out$shrinkage <- tune.tbl.i$shrinkage
  out <- c(out, other.settings$other.args)
  return(out)
}

gbm.predict <- function(mod, envs, other.settings) {
  requireNamespace("gbm", quietly = TRUE)
  if(inherits(envs, "SpatRaster") == TRUE) {
    pred <- terra::predict(envs, mod, type = "response", fun = predict)
  }else{
    # otherwise, envs is data frame, so return data frame of predicted values
    pred <- predict(mod, envs, type = "response")
  }
  return(pred)
}

gbm.ncoefs <- function(mod) {
  length(mod$var.names)
}

# no existing method in model object for variable importance
gbm.variable.importance <- function(mod) {
  NULL
}

#' @title ENMdetails gbm
#' @description This is the ENMdetails implementation for gradient boosting 
#' trees, implemented with the gbm package.
#' @export
enm.gbm <- ENMdetails(name = gbm.name, fun = gbm.fun, errors = gbm.errors,
                         msgs = gbm.msgs, args = gbm.args,
                         predict = gbm.predict, ncoefs = gbm.ncoefs, variable.importance = gbm.variable.importance)
