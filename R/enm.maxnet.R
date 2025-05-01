################################# #
# maxnet ENMdetails object ####
################################# #

maxnet.name <- "maxnet"

maxnet.fun <- maxnet::maxnet

maxnet.errors <- function(occs, envs, bg, tune.args, partitions, algorithm, 
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

maxnet.msgs <- function(tune.args, other.settings) {
  msg <- paste0("maxnet from maxnet package v", packageVersion('maxnet'))
  return(msg)
}

maxnet.args <- function(occs.z, bg.z, tune.tbl.i, other.settings) {
  out <- list()
  out$data <- rbind(occs.z, bg.z)
  out$p <- c(rep(1, nrow(occs.z)), rep(0, nrow(bg.z)))
  out$f <- maxnet::maxnet.formula(out$p, out$data, classes = tolower(tune.tbl.i$fc))
  out$regmult <- tune.tbl.i$rm
  # some models fail to converge if this parameter is not set to TRUE
  # usually the case with sparse datasets
  if(is.null(other.settings$addsamplestobackground)) {
    out$addsamplestobackground <- TRUE
  }else{
    out$addsamplestobackground <- other.settings$addsamplestobackground
  }
  out <- c(out, other.settings$other.args)
  return(out)
}

maxnet.predict <- function(mod, envs, other.settings) {
  requireNamespace("maxnet", quietly = TRUE)
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

maxnet.ncoefs <- function(mod) {
  length(mod$betas)
}

# no existing method in model object for variable importance
maxnet.variable.importance <- function(mod) {
  NULL
}

#' @title ENMdetails maxnet
#' @description This is the ENMdetails implementation for maxnet, the R version of
#' the Maxent algorithm. The configuration for running the model now includes addsamplestobackground = TRUE,
#' which explicitly adds presences to the background for model training, though as the current 
#' version of maxnet has this set to TRUE as default, behavior between ENMeval versions should not differ.
#' @export
enm.maxnet <- ENMdetails(name = maxnet.name, fun = maxnet.fun, errors = maxnet.errors,
                         msgs = maxnet.msgs, args = maxnet.args,
                         predict = maxnet.predict, ncoefs = maxnet.ncoefs, variable.importance = maxnet.variable.importance)
