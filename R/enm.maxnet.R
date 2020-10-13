################################# #
# maxnet ENMdetails object ####
################################# #

maxnet.name <- "maxnet"

maxnet.fun.train <- maxnet::maxnet

maxnet.fun.val <- maxnet.fun.train

maxnet.msgs <- function(tune.args, other.settings) {
  if(!("rm" %in% names(tune.args)) | !("fc" %in% names(tune.args))) {
    stop("Maxent settings must include 'rm' (regularization multiplier) and 'fc' (feature class) settings. See ?tune.args for details.")
  }else{
    if(!is.numeric(tune.args[["rm"]])) {
      stop("Please input numeric values for 'rm' settings for Maxent.")
    }
    all.fc <- unlist(sapply(1:5, function(x) apply(combn(c("L","Q","H","P","T"), x), 2, function(y) paste(y, collapse = ""))))
    if(any(!tune.args[["fc"]] %in% all.fc)) {
      stop("Please input accepted values for 'fc' settings for Maxent.")
    }
    msg <- paste0("maxnet from maxnet package v", packageVersion('maxnet'))
    return(msg)
  }
}

maxnet.args.train <- function(occs.z, bg.z, tune.i, other.settings) {
  out <- list()
  out$data <- rbind(occs.z, bg.z)
  out$p <- c(rep(1, nrow(occs.z)), rep(0, nrow(bg.z)))
  out$f <- maxnet::maxnet.formula(out$p, out$data, classes = tolower(tune.i$fc))
  out$regmult <- tune.i$rm
  # some models fail to converge if this parameter is not set to TRUE
  # usually the case with sparse datasets
  out$addsamplestobackground <- TRUE
  out <- c(out, other.settings$other.args)
  return(out)
}

maxnet.args.val <- function(occs.z, bg.z, tune.i, other.settings, mod = NULL) {
  maxnet.args.train(occs.z, bg.z, tune.i, other.settings)
}

maxnet.predict <- function(mod, envs, doClamp, other.settings) {
  # function to generate a prediction Raster* when raster data is specified as envs,
  # and a prediction data frame when a data frame is specified as envs
  if(inherits(envs, "BasicRaster") == TRUE) {
    envs.n <- raster::nlayers(envs)
    envs.pts <- raster::getValues(envs) %>% as.data.frame()
    mxnet.p <- dismo::predict(mod, envs.pts, type = other.settings$pred.type, 
                              clamp = doClamp,  other.settings$other.args)
    envs.pts[as.numeric(row.names(mxnet.p)), "pred"] <- mxnet.p
    pred <- raster::rasterFromXYZ(cbind(raster::coordinates(envs), envs.pts$pred), res=raster::res(envs), crs = raster::crs(envs)) 
  }else{
    # otherwise, envs is data frame, so return data frame of predicted values
    pred <- dismo::predict(mod, envs, type = other.settings$pred.type, na.rm = TRUE, clamp = doClamp, other.settings$other.args) %>% as.numeric()
  }
  return(pred)
}

maxnet.ncoefs <- function(mod) {
  length(mod$betas)
}

# no existing method in model object for variable importance
maxnet.varimp <- function(mod) {
  NULL
}

#' @export
enm.maxnet <- ENMdetails(name = maxnet.name, fun.train = maxnet.fun.train, fun.val = maxnet.fun.val,
                         msgs = maxnet.msgs, args.train = maxnet.args.train, args.val = maxnet.args.val,
                         predict = maxnet.predict, ncoefs = maxnet.ncoefs, varimp = maxnet.varimp)
