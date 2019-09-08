################################# #
# brt ENMdetails object ####
################################# #

fun <- dismo::gbm.step

pkgs <- c("dismo", "raster", "gbm")

msgs <- function(tune.args) {
  if(!all("tree.complexity" %in% names(tune.args), "learning.rate" %in% names(tune.args), "bag.fraction" %in% names(tune.args))) {
    stop("BRT settings must include 'tree.complexity', 'learning.rate', and 'bag.fraction'.")
  }
  # construct user message with version info
  msg <- paste("Boosted regression trees (BRTs) using the gbm.step() function from gbm package v.", 
               packageVersion('gbm'), "and dismo package v.", packageVersion('dismo')) 
  return(msg)
}

args <- function(occs.vals, bg.vals, tune.tbl.i, other.args) {
  out <- list()
  d <- rbind(occs.vals, bg.vals)
  p <- c(rep(1, nrow(occs.vals)), rep(0, nrow(bg.vals)))
  out$data <- cbind(p, d)
  out$tree.complexity <- tune.tbl.i$tree.complexity
  out$learning.rate <- tune.tbl.i$learning.rate
  out$bag.fraction <- tune.tbl.i$bag.fraction
  out$gbm.x <- 2:ncol(out$data)
  out$gbm.y <- 1
  out$silent <- TRUE
  out <- c(out, other.args)
  return(out)
}

auc <- function(occs.vals, bg.vals, mod, other.args, doClamp) {
  e <- dismo::evaluate(occs.vals, bg.vals, mod, n.trees = length(mod$trees))@auc
  return(e)
}

pred <- function(mod, envs, other.args, doClamp) {
  if(inherits(envs, "BasicRaster") == TRUE) {
    pred <- raster::predict(envs, mod, type = "response", n.trees = mod$gbm.call$best.trees, na.rm = TRUE)
  }else{
    pred <- dismo::predict(mod, envs, type = "response", n.trees = mod$gbm.call$best.trees, na.rm = TRUE)  
  }
  return(pred)
}

nparams <- function(mod) {
  # as no L1 regularization occurs, no parameters are dropped
  length(mod$var.names)
}

enm.brt <- ENMdetails(fun = fun, pkgs = pkgs, msgs = msgs, 
                      args = args, auc = auc, pred = pred, 
                      nparams = nparams)