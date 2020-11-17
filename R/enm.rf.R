################################# #
# random forest ENMdetails object ####
################################# #

rf.name <- "randomForest"

rf.fun.train <- ranger::ranger

rf.fun.val <- rf.fun.train

rf.msgs <- function(tune.args, other.settings) {
  if(!all("ntree" %in% names(tune.args), "mtry" %in% names(tune.args))) {
    stop("Random forest settings must include 'ntree' and 'mtry'. See ?tune.args for details")
  }
  # construct user message with version info
  msg <- paste0("Random forest using the randomForest() function from randomForest package v", 
                packageVersion('randomForest')) 
  return(msg)
}

rf.args.train <- function(occs.z, bg.z, tune.i, other.settings) {
  out <- list()
  d <- rbind(occs.z, bg.z)
  p <- c(rep(1, nrow(occs.z)), rep(0, nrow(bg.z)))
  out$data <- cbind(p, d)
  out$formula <- formula(out$data)
  out$num.trees <- tune.i$num.trees
  out$mtry <- tune.i$mtry
  out$sample.fraction <- tune.i$sample.fraction
  out$importance <- "permutation"
  out$case.weights <- c(rep(1, nrow(occs.z)), rep(nrow(occs.z)/nrow(bg.z), nrow(bg.z)))
  out <- c(out, other.settings$other.args)
  return(out)
}

rf.args.val <- function(occs.z, bg.z, tune.i, other.settings, mod = NULL) {
  rf.args.train(occs.z, bg.z, tune.i, other.settings)
}

rf.predict <- function(mod, envs, other.settings) {
  if(inherits(envs, "BasicRaster") == TRUE) {
    pred <- raster::predict(envs, mod, type = "response", fun = function(model, ...) predict(model, ...)$predictions, na.rm = TRUE)
  }else{
    pred <- dismo::predict(mod, envs, type = "response", na.rm = TRUE)$predictions
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
enm.randomForest <- ENMdetails(name = rf.name, fun.train = rf.fun.train, fun.val = rf.fun.val,
                               msgs = rf.msgs, args.train = rf.args.train, args.val = rf.args.val,
                               predict = rf.predict, ncoefs = rf.ncoefs, varimp = rf.varimp)
