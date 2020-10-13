################################# #
# random forest ENMdetails object ####
################################# #

rf.name <- "randomForest"

rf.fun.train <- randomForest::randomForest

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
  out$x <- rbind(occs.z, bg.z)
  out$y <- factor(c(rep(1, nrow(occs.z)), rep(0, nrow(bg.z))))
  out$ntree <- tune.i$ntree
  out$mtry <- tune.i$mtry
  # out$sampsize <- tune.i$sampsize
  # out$nodesize <- tune.i$nodesize
  out$importance <- TRUE
  out <- c(out, other.settings$other.args)
  return(out)
}

rf.args.val <- function(occs.z, bg.z, tune.i, other.settings, mod = NULL) {
  rf.args.train(occs.z, bg.z, tune.i, other.settings)
}

rf.predict <- function(mod, envs, clamp, other.settings) {
  if(inherits(envs, "BasicRaster") == TRUE) {
    pred <- raster::predict(envs, mod, type = "prob", na.rm = TRUE)
  }else{
    pred <- dismo::predict(mod, envs, type = "prob", na.rm = TRUE)[,2]
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
