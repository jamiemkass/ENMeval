################################# #
# boostedRegressionTrees ENMdetails object ####
################################# #

brt.name <- "boostedRegressionTrees"

brt.fun <- dismo::gbm.step

brt.msgs <- function(tune.args, other.settings) {
  if(!all("tc" %in% names(tune.args), "lr" %in% names(tune.args))) {
    stop('Boosted regression trees settings must include "tc" (tree complexity) and "lr" (learning rate). See ?tune.args for details')
  }
  # construct user message with version info
  msg <- paste0("Boosted regression trees tuned with gbm.step() from dismo package v", packageVersion('dismo'))
  return(msg)
}

brt.args <- function(occs.z, bg.z, tune.i, other.settings) {
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
  # grab partition groups from other.settings for gbm.step to do CV internally with custom folds
  out$fold.vector <- c(as.numeric(as.character(other.settings$occs.grp)), as.numeric(as.character(other.settings$bg.grp)))
  out$n.folds <- length(unique(other.settings$occs.grp))
  out$keep.fold.vector <- TRUE
  out$keep.fold.fit <- TRUE
  out$keep.fold.models <- TRUE
  out$silent <- TRUE
  out <- c(out, other.settings$other.args)
  return(out)
}

brt.predict <- function(mod, envs, other.settings) {
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

#' @title ENMdetails: Boosted Regression Trees
#' @description Boosted regression trees are implemented with the function \code{gbm.step} from the \code{dismo}
#' package. Details on this algorithm sand its uses in ecology and evolution can be found in Elith et al. (2008).
#' The function \code{gbm.step} runs a cross-validation procedure to determine the optimal number of trees by minimizing
#' holdout deviance (see ?dismo::gbm.step for details). As the number of trees (\code{n.trees}) is decided by \code{gbm.step} 
#' for each combination of the tunable parameters, \code{n.trees} is not itself tunable in ENMeval.\cr
#' \cr
#' The tunable parameters in this implementation are:\cr
#' 1. \code{tree.complexity}: the number of variable interactions to allow (a value of 1 means no interactions, 
#' or an additive model); analogous to \code{interaction.depth} in the \code{gbm} package\cr
#' 2. \code{learning.rate}: the contribution of each tree to the model (a smaller value increases the number of 
#' trees required, and smaller is generally preferable); analogous to \code{shrinkage} in the \code{gbm} package\cr
#' 
#' In order to run \code{gbm.step} with custom partitions, the partitions chosen by the user in ENMeval are specified 
#' with argument \code{fold.vector}. See Details for internal changes to partitions that are made for random k-fold 
#' so the function can run. In addition, \code{site.weights} are included for background records to balance the occurrence 
#' and background data during model fitting, equal to the number of occurrence records divided by the number of background records.
#' Besides these specifications, the function \code{gbm.step} is run with all defaults
#' 
#' Jackknife (leave-one-out) partitioning is not available for boosted regression trees mainly because the models often fail
#' with small datasets, and even if they do succeed to completion these models in general are not advisable for data-poor species
#' (pers. comm. Jane Elith). 
#' 
#' @details As \code{gbm.step} requires that each class is partitioned for cross validation, both occurrence and background
#' records must have associated folds for the optimization procedure. Thus, for random k-fold partitioning, background records
#' are internally partitioned for \code{gbm.step}. In addition, validation models are built internally by \code{gbm.step}, and 
#' therefore are not built sequentially within \code{tune.enm.R} as for other algorithms -- this results in a implementation that
#' avoids running the optimization procedure multiple times.
#' @references Elith, J., Leathwick, J. R., & Hastie, T. (2008). A working guide to boosted regression trees. \emph{Journal of Animal Ecology}, 77(4), 802-813.
#' @export
enm.boostedRegressionTrees <- ENMdetails(name = brt.name, fun = brt.fun,
                                         msgs = brt.msgs, args = brt.args,
                                         predict = brt.predict, ncoefs = brt.ncoefs, varimp = brt.varimp)
