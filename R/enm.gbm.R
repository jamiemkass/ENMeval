################################# #
# rf ENMdetails object ####
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
  #implementation of down-sampled RF from Valavi et al. 2021
  out <- list()
  out$formula <- formula(p ~.)
  out$data <- rbind(occs.z, bg.z)
  p <- as.factor(c(rep(1, nrow(occs.z)), rep(0, nrow(bg.z))))
  out$data <- cbind(p, out$data)
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
    pred <- terra::predict(envs, mod, type = "prob", fun = predict,
                           other.settings$other.args)[[2]]
  }else{
    # otherwise, envs is data frame, so return data frame of predicted values
    pred <- predict(mod, envs, type = "prob", 
                    other.settings$other.args)[,2] |> as.numeric()
  }
  return(pred)
}

gbm.ncoefs <- function(mod) {
  nrow(mod$importance)
}

# no existing method in model object for variable importance
gbm.variable.importance <- function(mod) {
  # remove mean decrease in accuracy for absences (1st column)
  # and mean decrease in accuracy over all classes (3rd column)
  # this leaves mean decrease in accuracy for presences and mean decrease
  # in Gini index
  imp <- data.frame(mod$importance[,c(-1,-3)])
  names(imp)[1] <- "MeanDecreaseAccuracyPres"
  imp <- dplyr::arrange(imp, dplyr::desc(MeanDecreaseAccuracyPres))
  return(imp)
}

#' @title ENMdetails rf
#' @description This is the ENMdetails implementation for random forest, the R version of
#' the Maxent algorithm. The configuration for running the model now includes addsamplestobackground = TRUE,
#' which explicitly adds presences to the background for model training, though as the current 
#' version of rf has this set to TRUE as default, behavior between ENMeval versions should not differ.
#' @export
enm.rf <- ENMdetails(name = gbm.name, fun = gbm.fun, errors = gbm.errors,
                         msgs = gbm.msgs, args = gbm.args,
                         predict = gbm.predict, ncoefs = gbm.ncoefs, variable.importance = gbm.variable.importance)
