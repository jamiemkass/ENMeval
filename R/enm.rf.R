################################# #
# random forest ENMdetails object ####
################################# #

name <- "randomForest"

fun <- randomForest::randomForest

msgs <- function(tune.args) {
  if(!all("ntree" %in% names(tune.args), "mtry" %in% names(tune.args))) {
    stop("Random forest settings must include 'ntree' and 'mtry'. See ?tune.args for details")
  }
  # construct user message with version info
  msg <- paste0("Random forest using the randomForest() function from randomForest package v", 
                packageVersion('randomForest')) 
  return(msg)
}

args <- function(occs.z, bg.z, tune.i, other.settings) {
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

predict <- function(mod, envs, other.settings) {
  if(inherits(envs, "BasicRaster") == TRUE) {
    pred <- raster::predict(envs, mod, type = "prob", na.rm = TRUE)
  }else{
    pred <- dismo::predict(mod, envs, type = "prob", na.rm = TRUE)[,2]
  }
  return(pred)
}

ncoefs <- function(mod) {
  # as no L1 regularization occurs, no parameters are dropped
  return(nrow(mod$importance))
}

# see ?randomForest for detailed description of how to interpret the table
varimp <- function(mod) {
  return(mod$importance)
}

#' @export
enm.randomForest <- ENMdetails(name = name, fun = fun, msgs = msgs, args = args, 
                      predict = predict, ncoefs = ncoefs, varimp = varimp)
