################################# #
# random forest ENMdetails object ####
################################# #

rf.name <- "randomForest"

rf.fun <- ranger::ranger

rf.msgs <- function(tune.args, other.settings) {
  if(!all("num.trees" %in% names(tune.args), "mtry" %in% names(tune.args))) {
    stop("Random forest settings must include 'num.trees' and 'mtry'. See ?tune.args for details")
  }
  # construct user message with version info
  msg <- paste0("Random forest using the ranger() function from ranger package v", 
                packageVersion('ranger')) 
  return(msg)
}

rf.args <- function(occs.z, bg.z, tune.i, other.settings) {
  out <- list()
  d <- rbind(occs.z, bg.z)
  p <- c(rep(1, nrow(occs.z)), rep(0, nrow(bg.z)))
  out$data <- cbind(p, d)
  out$formula <- formula(out$data)
  out$num.trees <- tune.i$num.trees
  out$mtry <- tune.i$mtry
  out$sample.fraction <- nrow(occs.z) / nrow(bg.z)
  out$case.weights <- ifelse(out$data$p == 1, 1, nrow(occs.z) / nrow(bg.z))
  out$probability <- TRUE
  out$importance <- "permutation"
  out <- c(out, other.settings$other.args)
  return(out)
}

rf.predict <- function(mod, envs, other.settings) {
  if(inherits(envs, "BasicRaster") == TRUE) {
    pred <- raster::predict(envs, mod, type = "response", fun = function(model, ...) predict(model, ...)$predictions, na.rm = TRUE)
  }else{
    pred <- dismo::predict(mod, envs, type = "response", na.rm = TRUE)$predictions[,1]
  }
  return(pred)
}

rf.ncoefs <- function(mod) {
  # as no L1 regularization occurs, no parameters are dropped
  return(nrow(mod$importance))
}

# see ?randomForest for detailed description of how to interpret the table
rf.varimp <- function(mod) {
  return(mod$importance)
}

#' @export
enm.randomForest <- ENMdetails(name = rf.name, fun = rf.fun,
                               msgs = rf.msgs, args = rf.args,
                               predict = rf.predict, ncoefs = rf.ncoefs, varimp = rf.varimp)
