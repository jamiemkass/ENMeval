################################# #
# boostedRegressionTrees ENMdetails object ####
################################# #

name <- "boostedRegressionTrees"

fun <- gbm::gbm

pkgs <- c("dismo", "raster", "gbm")

msgs <- function(tune.args) {
  if(!all("tc" %in% names(tune.args), "lr" %in% names(tune.args))) {
    stop('Boosted regression trees settings must include "ntree", "tc" (tree complexity, or "interaction depth"), "lr" (learning rate, or "shrinkage"). See ?tune.args for details')
  }
  # construct user message with version info
  msg <- paste0("Boosted regression trees using the gbm() function from gbm package v", packageVersion('gbm'))
  return(msg)
}

args <- function(occs.z, bg.z, tune.i, other.settings) {
  out <- list()
  d <- rbind(occs.z, bg.z)
  p <- c(rep(1, nrow(occs.z)), rep(0, nrow(bg.z)))
  out$data <- cbind(p, d)
  out$formula <- formula(out$data)
  out$n.trees <- tune.i$ntree
  out$interaction.depth <- tune.i$tc
  out$shrinkage <- tune.i$lr
  out$distribution <- "bernoulli"
  out <- c(out, other.settings$other.args)
  return(out)
}

predict <- function(mod, envs, other.settings) {
  if(inherits(envs, "BasicRaster") == TRUE) {
    pred <- raster::predict(envs, mod, type = "response", n.trees = length(mod$trees), na.rm = TRUE)
  }else{
    pred <- dismo::predict(mod, envs, type = "response", n.trees = length(mod$trees), na.rm = TRUE)  
  }
  return(pred)
}

nparams <- function(mod) {
  # as no L1 regularization occurs, no parameters are dropped
  length(mod$var.names)
}

#' @export
enm.boostedRegressionTrees <- ENMdetails(name = name, fun = fun, pkgs = pkgs, msgs = msgs, args = args, 
                      predict = predict, nparams = nparams)
