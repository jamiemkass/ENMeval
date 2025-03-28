################################# #
# rf ENMdetails object ####
################################# #

rf.name <- "rf"

rf.fun <- randomForest::randomForest

rf.errors <- function(occs, envs, bg, tune.args, partitions, algorithm, 
                          partition.settings, other.settings, 
                          categoricals, doClamp, clamp.directions) {
  if(!("mtry" %in% names(tune.args))) {
    stop("RF settings must include 'mtry' (number of variables randomly sampled at each tree split). See ?tune.args for details.")
  }
  if(any(tune.args$mtry <= 0)) {
    stop("Please input positive integer values for 'mtry' settings for RF.")
  }
}

rf.msgs <- function(tune.args, other.settings) {
  msg <- paste0("randomForest from randomForest package v", packageVersion('randomForest'))
  return(msg)
}

rf.args <- function(occs.z, bg.z, tune.tbl.i, other.settings) {
  #implementation of down-sampled RF from Valavi et al. 2021
  out <- list()
  out$formula <- formula(p ~.)
  out$data <- rbind(occs.z, bg.z)
  p <- as.factor(c(rep(1, nrow(occs.z)), rep(0, nrow(bg.z))))
  out$data <- cbind(p, out$data)
  prNum <- as.numeric(table(p)["1"])
  bgNum <- as.numeric(table(p)["0"])
  out$sampsize <- c("0" = prNum, "1" = prNum)
  out$replace <- TRUE
  out$ntree <- 1000
  out$importance <- TRUE
  out <- c(out, other.settings$other.args)
  return(out)
}

rf.predict <- function(mod, envs, other.settings) {
  requireNamespace("randomForest", quietly = TRUE)
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

rf.ncoefs <- function(mod) {
  nrow(mod$importance)
}

# no existing method in model object for variable importance
rf.variable.importance <- function(mod) {
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
enm.rf <- ENMdetails(name = rf.name, fun = rf.fun, errors = rf.errors,
                         msgs = rf.msgs, args = rf.args,
                         predict = rf.predict, ncoefs = rf.ncoefs, variable.importance = rf.variable.importance)
