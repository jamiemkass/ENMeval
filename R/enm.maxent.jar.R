################################# #
# maxent.jar ENMdetails object ####
################################# #

name <- "maxent.jar"

fun <- dismo::maxent

pkgs <- c("dismo", "raster", "rJava")

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
  }
  if(is.null(getOption('dismo_rJavaLoaded'))) {
    # to avoid trouble on macs
    Sys.setenv(NOAWT=TRUE)
    if ( requireNamespace('rJava') ) {
      rJava::.jpackage('dismo')
      options(dismo_rJavaLoaded=TRUE)
    } else {
      stop('rJava cannot be loaded')
    }
  }
  mxe <- rJava::.jnew("meversion") 
  v <- try(rJava::.jcall(mxe, "S", "meversion"))
  msg <- paste0("maxent.jar v", v, " from dismo package v", packageVersion('dismo'))
  return(msg)
}

args <- function(occs.z, bg.z, tune.i, other.settings) {
  out <- list()
  out$x <- rbind(occs.z, bg.z)
  out$p <- c(rep(1, nrow(occs.z)), rep(0, nrow(bg.z)))
  out$args <- c("noaddsamplestobackground", "noremoveDuplicates", "noautofeature")
  if(!grepl("L", tune.i$fc)) out$args <- c(out$args, "nolinear")
  if(!grepl("Q", tune.i$fc)) out$args <- c(out$args, "noquadratic")
  if(!grepl("H", tune.i$fc)) out$args <- c(out$args, "nohinge")
  if(!grepl("P", tune.i$fc)) out$args <- c(out$args, "noproduct")
  if(!grepl("T", tune.i$fc)) out$args <- c(out$args, "nothreshold") else out$args <- c(out$args, "threshold")
  out$args <- c(out$args, paste0("betamultiplier=", tune.i$rm, sep=""))
  out <- c(out, other.settings$other.args)
  return(out)
}

eval.train <- function(occs.xy, bg.xy, occs.z, bg.z, mod.full, mod.full.pred, envs, other.settings) {
  # training AUC
  clamp <- ifelse(other.settings$clamp == TRUE, "doclamp=true", "doclamp=false")
  output.format <- paste0("outputformat=", other.settings$pred.type)
  e <- dismo::evaluate(occs.z, bg.z, mod.full, args = c(output.format, clamp))
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

eval.validate <- function(occs.val.xy, occs.train.xy, bg.xy, occs.train.z, occs.val.z, bg.z, bg.val.z, mod.k, nk, other.settings) {
  ## validation AUC
  # calculate auc on validation data: validation occurrences are evaluated on full background, as in Radosavljevic & Anderson 2014
  # for auc.diff calculation, to perform the subtraction, it is essential that both stats are calculated over the same background
  clamp <- ifelse(other.settings$clamp == TRUE, "doclamp=true", "doclamp=false")
  output.format <- paste0("outputformat=", other.settings$pred.type)
  e.train <- dismo::evaluate(occs.train.z, bg.z, mod.k, args = c(output.format, clamp))
  auc.train <- e.train@auc
  # AUC validation
  # calculate AUC on validation data: if training and validation occurrences are evaluated on same background (full), auc.diff can
  # also be calculated (Radosavljevic & Anderson 2014); if validation occurrences are evaluated on partitioned background, auc.diff is NULL
  if(other.settings$validation.bg == "full") {
    e.val <- dismo::evaluate(occs.val.z, bg.z, mod.k, args = c(output.format, clamp))
    auc.val <- e.val@auc
    # calculate auc diff as auc train (partition not k) minus auc validation (partition k)
    auc.diff <- auc.train - auc.val
  } else if(other.settings$validation.bg == "partition") {
    e.val <- dismo::evaluate(occs.val.z, bg.val.z, mod.k, args = c(output.format, clamp))
    auc.val <- e.val@auc
    auc.diff <- NA
  }
  
  if(other.settings$abs.auc.diff == TRUE & !is.null(auc.diff)) auc.diff <- abs(auc.diff)
  
  ## omission rates
  # get model predictions for training and validation data
  occs.train.pred <- enm.maxent.jar@pred(mod.k, occs.train.z, other.settings)
  occs.val.pred <- enm.maxent.jar@pred(mod.k, occs.val.z, other.settings)
  # get minimum training presence threshold (expected no omission)
  min.train.thr <- min(occs.train.pred)
  or.mtp <- mean(occs.val.pred < min.train.thr)
  # get 10 percentile training presence threshold (expected 0.1 omission)
  pct10.train.thr <- calc.10p.trainThresh(occs.train.pred)
  or.10p <- mean(occs.val.pred < pct10.train.thr)
  
  ## validation CBI
  if(other.settings$cbi.cv == TRUE) {
    if(other.settings$validation.bg == "full") {
      mod.k.pred <- enm.maxent.jar@pred(mod.k, bg.z, other.settings)  
    }else if(other.settings$validation.bg == "partition") {
      mod.k.pred <- enm.maxent.jar@pred(mod.k, bg.val.z, other.settings)  
    }
    cbi.val <- ecospat::ecospat.boyce(mod.k.pred, occs.val.pred, PEplot = FALSE)$Spearman.cor
  }else{
    cbi.val <- NA
  }  
  # gather all evaluation statistics for k
  out.df <- data.frame(auc.val = auc.val, auc.diff = auc.diff, cbi.val = cbi.val, or.mtp = or.mtp, or.10p = or.10p)
  
  return(out.df)
}

pred <- function(mod, envs, other.settings) {
  type.arg <- paste("outputformat", other.settings$pred.type, sep = "=")
  clamp.arg <- ifelse(other.settings$clamp == TRUE, "doclamp=true", "doclamp=false")
  pred <- dismo::predict(mod, envs, args = c(type.arg, clamp.arg), na.rm = TRUE)
  return(pred)
}

nparams <- function(mod) {
  lambdas <- mod@lambdas[1:(length(mod@lambdas)-4)]
  countNonZeroParams <- function(x) if(strsplit(x, split=", ")[[1]][2] != '0.0') 1
  np <- sum(unlist(sapply(lambdas, countNonZeroParams)))
  return(np)
}

#' @export
enm.maxent.jar <- ENMdetails(name = name, fun = fun, pkgs = pkgs, msgs = msgs, args = args, 
                             eval.train = eval.train, eval.validate = eval.validate, 
                             pred = pred, nparams = nparams)
