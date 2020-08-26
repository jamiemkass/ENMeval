################################# #
# random forest ENMdetails object ####
################################# #

name <- "random forest"

fun <- randomForest::randomForest

pkgs <- c("randomForest", "dismo", "raster")

msgs <- function(tune.args) {
  if(!all("ntree" %in% names(tune.args), "mtry" %in% names(tune.args))) {
    stop("RF settings must include 'ntree' and 'mtry'. See ?tune.args for details")
  }
  # construct user message with version info
  msg <- paste0("random forest using the randomForest() function from randomForest package v", 
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

evaluate <- function(occs.z, bg.z, mod, other.settings) {
  dismo::evaluate(occs.z, bg.z, mod, type = "prob", index = 2)
}

predict <- function(mod, envs, other.settings) {
  if(inherits(envs, "BasicRaster") == TRUE) {
    pred <- raster::predict(envs, mod, type = "prob", index = 2, na.rm = TRUE)
  }else{
    pred <- dismo::predict(mod, envs, type = "prob", index = 2, na.rm = TRUE)  
  }
  return(pred)
}

nparams <- function(mod) {
  # as no L1 regularization occurs, no parameters are dropped
  nrow(mod$importance)
}

varimp <- function(mod) {
  mod@importance
}

#' @export
enm.rf <- ENMdetails(name = name, fun = fun, pkgs = pkgs, msgs = msgs, args = args, 
                      evaluate = evaluate, predict = predict, nparams = nparams)
