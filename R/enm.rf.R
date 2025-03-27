################################# #
# rf ENMdetails object ####
################################# #

rf.name <- "random_forest"

rf.fun <- randomForest::randomForest

rf.errors <- function(occs, envs, bg, tune.args, partitions, algorithm, 
                          partition.settings, other.settings, 
                          categoricals, doClamp, clamp.directions) {
  if(!("rm" %in% names(tune.args)) | !("fc" %in% names(tune.args))) {
    stop("Maxent settings must include 'rm' (regularization multiplier) and 'fc' (feature class) settings. See ?tune.args for details.")
  }else{
    if(!is.numeric(tune.args[["rm"]])) {
      stop("Please input numeric values for 'rm' settings for maxnet.")
    }
    all.fc <- unlist(sapply(1:5, function(x) apply(combn(c("L","Q","H","P","T"), x), 2, function(y) paste(y, collapse = ""))))
    if(any(!tune.args[["fc"]] %in% all.fc)) {
      stop("Please input accepted values for 'fc' settings for maxnet.")
    }
  }
  if(any(tune.args$rm <= 0)) {
    stop("Please input a positive value for 'rm' settings for maxnet.")
  }
}

rf.msgs <- function(tune.args, other.settings) {
  msg <- paste0("maxnet from maxnet package v", packageVersion('maxnet'))
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
  out <- c(out, other.settings$other.args)
  return(out)
}

rf.predict <- function(mod, envs, other.settings) {
  requireNamespace("randomForest", quietly = TRUE)
  # function to generate a prediction Raster* when raster data is specified as envs,
  # and a prediction data frame when a data frame is specified as envs
  if(inherits(envs, "SpatRaster") == TRUE) {
    pred <- maxnet.predictRaster(mod, envs, other.settings$pred.type,
                                 other.settings$doClamp, 
                                 other.settings$other.args)
  }else{
    # otherwise, envs is data frame, so return data frame of predicted values
    pred <- predict(mod, envs, type = other.settings$pred.type, na.rm = TRUE, 
                    clamp = other.settings$doClamp, 
                    other.settings$other.args) |> as.numeric()
  }
  return(pred)
}

rf.ncoefs <- function(mod) {
  length(mod$betas)
}

# no existing method in model object for variable importance
rf.variable.importance <- function(mod) {
  NULL
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
