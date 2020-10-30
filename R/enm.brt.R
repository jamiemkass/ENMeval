################################# #
# boostedRegressionTrees ENMdetails object ####
################################# #

brt.name <- "boostedRegressionTrees"

brt.fun.train <- dismo::gbm.step 

brt.fun.val <- gbm::gbm

brt.msgs <- function(tune.args, other.settings) {
  if(!all("tc" %in% names(tune.args), "lr" %in% names(tune.args))) {
    stop('Boosted regression trees settings must include "tc" (tree complexity, or "interaction depth") and "lr" (learning rate, or "shrinkage"). See ?tune.args for details')
  }
  # construct user message with version info
  msg <- paste0("Boosted regression trees tuned for optimal tree number with gbm.step() from dismo package v", packageVersion('dismo'), " and validated with gbm() function from gbm package v", packageVersion('gbm'))
  return(msg)
}

brt.args.train <- function(occs.z, bg.z, tune.i, other.settings) {
  out <- list()
  d <- rbind(occs.z, bg.z)
  p <- c(rep(1, nrow(occs.z)), rep(0, nrow(bg.z)))
  out$data <- cbind(p, d)
  out$tree.complexity <- tune.i$tc
  out$learning.rate <- tune.i$lr
  # out$bag.fraction <- tune.i$bag.fraction
  out$gbm.x <- 2:ncol(out$data)
  out$gbm.y <- 1
  out$family <- "bernoulli"
  out$site.weights <- c(rep(1, nrow(occs.z)), rep(nrow(occs.z)/nrow(bg.z), nrow(bg.z)))
  out$silent <- TRUE
  out <- c(out, other.settings$other.args)
  return(out)
}

brt.args.val <- function(occs.z, bg.z, tune.i, other.settings, mod.full = NULL) {
  out <- list()
  d <- rbind(occs.z, bg.z)
  p <- c(rep(1, nrow(occs.z)), rep(0, nrow(bg.z)))
  out$data <- cbind(p, d)
  out$formula <- formula(out$data)
  # get number of trees from existing mod tuned with gbm.step
  out$n.trees <- mod.full$n.trees
  out$interaction.depth <- tune.i$tc
  out$shrinkage <- tune.i$lr
  out$distribution <- "bernoulli" 
  out$weights <- c(rep(1, nrow(occs.z)), rep(nrow(occs.z)/nrow(bg.z), nrow(bg.z)))
  out <- c(out, other.settings$other.args)
  return(out)
}

brt.predict <- function(mod, envs, doClamp, other.settings) {
  if(inherits(envs, "BasicRaster") == TRUE) {
    pred <- raster::predict(envs, mod, type = "response", n.trees = length(mod$trees), na.rm = TRUE)
  }else{
    pred <- dismo::predict(mod, envs, type = "response", n.trees = length(mod$trees), na.rm = TRUE)  
  }
  return(pred)
}

brt.ncoefs <- function(mod) {
  # as no L1 regularization occurs, no coefficients are dropped
  return(length(mod$var.names))
}

# see ?summary.gbm for detailed description of how to interpret the table
brt.varimp <- function(mod) {
  return(summary(mod, plotit=FALSE))
}

#' @export
enm.boostedRegressionTrees <- ENMdetails(name = brt.name, fun.train = brt.fun.train, fun.val = brt.fun.val,
                                         msgs = brt.msgs, args.train = brt.args.train, args.val = brt.args.val,
                                         predict = brt.predict, ncoefs = brt.ncoefs, varimp = brt.varimp)
