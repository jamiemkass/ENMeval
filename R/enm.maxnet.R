################################# #
# maxnet ENMdetails object ####
################################# #

name <- "maxnet"

fun <- maxnet::maxnet

pkgs <- c("dismo", "raster", "maxnet")

msgs <- function(tune.args) {
  if(!("rm" %in% names(tune.args)) | !("fc" %in% names(tune.args))) {
    stop("For Maxent, please specify both 'rm' and 'fc' settings. See ?tune.args for help.")
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

args <- function(occs.z, bg.z, tune.i, other.settings) {
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

eval.train <- function(occs.xy, bg.xy, occs.z, bg.z, mod.full, mod.full.pred, envs, other.settings) {
  # training AUC
  e <- dismo::evaluate(occs.z, bg.z, mod.full, clamp = other.settings$clamp, type = other.settings$pred.type)
  auc.train <- e@auc
  # training CBI
  if(!is.null(envs)) {
    # if raster envs exists, use it to calculate CBI
    cbi.train.out <- ecospat::ecospat.boyce(mod.full.pred, occs.xy, PEplot = FALSE)
  }else{
    # if no raster envs exists, calculate CBI with occs and bg predictions
    d.full.pred <- dplyr::bind_rows(occs.z, bg.z)
    d.full.pred$pred <- mod.full.pred
    d.full.pred$pb <- c(rep(1, nrow(occs.z)), rep(0, nrow(bg.z)))
    cbi.train.out <- ecospat::ecospat.boyce(d.full.pred$pred, d.full.pred[d.full.pred$pb == 1, "pred"], PEplot = FALSE)
  }
  cbi.train <- cbi.train.out$Spearman.cor
  out.df <- data.frame(auc.train = auc.train, cbi.train = cbi.train)
  return(out.df)
}

eval.validate <- function(occs.test.xy, occs.train.xy, bg.xy, occs.train.z, occs.val.z, bg.z, mod.k, nk, envs, other.settings) {
  ## validation AUC
  # calculate auc on validation data: test occurrences are evaluated on full background, as in Radosavljevic & Anderson 2014
  # for auc.diff calculation, do perform the subtraction, it is essential that both stats are calculated over the same background
  e.train <- dismo::evaluate(occs.train.z, bg.z, mod.k, clamp = other.settings$clamp, type = other.settings$pred.type)
  e.test <- dismo::evaluate(occs.val.z, bg.z, mod.k, clamp = other.settings$clamp, type = other.settings$pred.type)
  auc.train <- e.train@auc
  auc.val <- e.test@auc
  # calculate auc diff as auc train (partition not k) minus auc test (partition k)
  auc.diff <- auc.train - auc.val
  if(other.settings$abs.auc.diff == TRUE) auc.diff <- abs(auc.diff)
  
  ## omission rates
  # get model predictions for training and validation data
  occs.train.pred <- enm.maxnet@pred(mod.k, occs.train.z, other.settings)
  occs.val.pred <- enm.maxnet@pred(mod.k, occs.val.z, other.settings)
  # get minimum training presence threshold (expected no omission)
  min.train.thr <- min(occs.train.pred)
  or.mtp <- mean(occs.val.pred < min.train.thr)
  # get 10 percentile training presence threshold (expected 0.1 omission)
  pct10.train.thr <- calc.10p.trainThresh(occs.train.pred)
  or.10p <- mean(occs.val.pred < pct10.train.thr)
  
  # custom stats
  # maxTSS.test <- max(e.test@TPR + e.test@TNR) - 1 
  # maxKappa.test <- max(e.test@kappa)
  
  ## validation CBI
  if(other.settings$cbi.cv == TRUE) {
    if(!is.null(envs)) {
      # predict to raster
      mod.k.pred <- enm.maxnet@pred(mod.k, envs, other.settings)
    }else{
      # use full background to approximate full model prediction
      mod.k.pred <- enm.maxnet@pred(mod.k, bg.z, other.settings)  
    }
    cbi.test <- ecospat::ecospat.boyce(mod.k.pred, occs.val.pred, PEplot = FALSE)
  }else{
    cbi.test <- NULL
  }
  
  # gather all evaluation statistics for k
  out.df <- data.frame(auc.val = auc.val, auc.diff = auc.diff, or.mtp = or.mtp, or.10p = or.10p)
  # maxTSS.test = maxTSS.test, maxKappa.test = maxKappa.test)
  if(!is.null(cbi.test)) out.df <- cbind(out.df, cbi.test = cbi.test$Spearman.cor)
  
  return(out.df)
}

pred <- function(mod, envs, other.settings) {
  # function to generate a prediction Raster* when raster data is specified as envs,
  # and a prediction data frame when a data frame is specified as envs
  if(inherits(envs, "BasicRaster") == TRUE) {
    envs.n <- raster::nlayers(envs)
    envs.pts <- raster::getValues(envs) %>% as.data.frame()
    mxnet.p <- predict(mod, envs.pts, type = other.settings$pred.type, 
                       clamp = other.settings$clamp,  other.settings$other.args)
    envs.pts[as.numeric(row.names(mxnet.p)), "pred"] <- mxnet.p
    pred <- raster::rasterFromXYZ(cbind(raster::coordinates(envs), envs.pts$pred), res=raster::res(envs), crs = raster::crs(envs)) 
  }else{
    # otherwise, envs is data frame, so return data frame of predicted values
    pred <- predict(mod, envs, type = other.settings$pred.type, 
                    clamp = other.settings$clamp, na.rm = TRUE, other.settings$other.args)
  }
  return(pred)
}

nparams <- function(mod) {
  length(mod$betas)
}

#' @export
enm.maxnet <- ENMdetails(name = name, fun = fun, pkgs = pkgs, msgs = msgs, args = args, 
                         eval.train = eval.train, eval.validate = eval.validate,
                         pred = pred, nparams = nparams)
