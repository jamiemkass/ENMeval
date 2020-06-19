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

args <- function(occs.vals, bg.vals, tune.tbl.i, other.settings) {
  out <- list()
  out$x <- rbind(occs.vals, bg.vals)
  out$p <- c(rep(1, nrow(occs.vals)), rep(0, nrow(bg.vals)))
  out$args <- c("noaddsamplestobackground", "noremoveDuplicates", "noautofeature")
  if(!grepl("L", tune.tbl.i$fc)) out$args <- c(out$args, "nolinear")
  if(!grepl("Q", tune.tbl.i$fc)) out$args <- c(out$args, "noquadratic")
  if(!grepl("H", tune.tbl.i$fc)) out$args <- c(out$args, "nohinge")
  if(!grepl("P", tune.tbl.i$fc)) out$args <- c(out$args, "noproduct")
  if(!grepl("T", tune.tbl.i$fc)) out$args <- c(out$args, "nothreshold") else out$args <- c(out$args, "threshold")
  out$args <- c(out$args, paste0("betamultiplier=", tune.tbl.i$rm, sep=""))
  out <- c(out, other.settings$other.args)
  return(out)
}

aic <- function(occs, nparam, mod.full.pred.all) {
  calc.aicc(occs, nparam, mod.full.pred.all)
}

eval.train <- function(occs.xy, bg.xy, occs.vals, bg.vals, mod.full, mod.full.pred, envs, other.settings) {
  # training AUC
  clamp <- ifelse(other.settings$doClamp == TRUE, "doclamp=true", "doclamp=false")
  output.format <- paste0("outputformat=", other.settings$pred.type)
  e <- dismo::evaluate(occs.vals, bg.vals, mod.full, args = c(output.format, clamp))
  auc.train <- e@auc
  # training CBI
  if(!is.null(envs)) {
    # if raster envs exists, use it to calculate CBI
    cbi.train.out <- ecospat::ecospat.boyce(mod.full.pred, occs.xy, PEplot = FALSE)
  }else{
    # if no raster envs exists, calculate CBI with occs and bg predictions
    d.full.pred <- dplyr::bind_rows(occs.vals, bg.vals)
    d.full.pred$pred <- mod.full.pred
    d.full.pred$pb <- c(rep(1, nrow(occs.vals)), rep(0, nrow(bg.vals)))
    cbi.train.out <- ecospat::ecospat.boyce(d.full.pred$pred, d.full.pred[d.full.pred$pb == 1, "pred"], PEplot = FALSE)
  }
  cbi.train <- cbi.train.out$Spearman.cor
  out.df <- data.frame(auc.train = auc.train, cbi.train = cbi.train)
  return(out.df)
}

eval.test <- function(occs.test.xy, occs.train.xy, bg.xy, occs.train.vals, occs.test.vals, bg.vals, mod.k, nk, envs, other.settings) {
  ## testing AUC
  # calculate auc on testing data: test occurrences are evaluated on full background, as in Radosavljevic & Anderson 2014
  # for auc.diff calculation, do perform the subtraction, it is essential that both stats are calculated over the same background
  clamp <- ifelse(other.settings$doClamp == TRUE, "doclamp=true", "doclamp=false")
  output.format <- paste0("outputformat=", other.settings$pred.type)
  e.train <- dismo::evaluate(occs.train.vals, bg.vals, mod.k, args = c(output.format, clamp))
  e.test <- dismo::evaluate(occs.test.vals, bg.vals, mod.k, args = c(output.format, clamp))
  auc.train <- e.train@auc
  auc.test <- e.test@auc
  # calculate auc diff as auc train (partition not k) minus auc test (partition k)
  auc.diff <- auc.train - auc.test
  if(other.settings$abs.auc.diff == TRUE) auc.diff <- abs(auc.diff)
  
  ## omission rates
  # get model predictions for training and testing data
  occs.train.pred <- enm.maxent.jar@pred(mod.k, occs.train.vals, other.settings)
  occs.test.pred <- enm.maxent.jar@pred(mod.k, occs.test.vals, other.settings)
  # get minimum training presence threshold (expected no omission)
  min.train.thr <- min(occs.train.pred)
  or.mtp <- mean(occs.test.pred < min.train.thr)
  # get 10 percentile training presence threshold (expected 0.1 omission)
  pct10.train.thr <- calc.10p.trainThresh(occs.train.pred)
  or.10p <- mean(occs.test.pred < pct10.train.thr)
  
  ## testing CBI
  if(other.settings$cbi.cv == TRUE) {
    if(other.settings$cbi.eval == "envs") {
      # use full model prediction over envs
      mod.k.pred <- enm.maxent.jar@pred(mod.k, envs, other.settings)
      cbi.test <- ecospat::ecospat.boyce(mod.k.pred, occs.test.xy, PEplot = FALSE)
    }else{
      # use full background to approximate full model prediction
      mod.k.pred <- enm.maxent.jar@pred(mod.k, bg.vals, other.settings)
      cbi.test <- ecospat::ecospat.boyce(mod.k.pred, occs.test.pred, PEplot = FALSE)
    }
  }else{
    cbi.test <- NULL
  }
  
  # gather all evaluation statistics for k
  out.df <- data.frame(auc.test = auc.test, auc.diff = auc.diff, or.mtp = or.mtp, or.10p = or.10p)
  if(!is.null(cbi.test)) out.df <- cbind(out.df, cbi.test = cbi.test$Spearman.cor)
  
  return(out.df)
}

pred <- function(mod, envs, other.settings) {
  type.arg <- paste("outputformat", other.settings$pred.type, sep = "=")
  clamp.arg <- ifelse(other.settings$doClamp == TRUE, "doclamp=true", "doclamp=false")
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
enm.maxent.jar <- ENMdetails(name = name, fun = fun, pkgs = pkgs, msgs = msgs, 
                        args = args, aic = aic, 
                        eval.train = eval.train, eval.test = eval.test, 
                        pred = pred, nparams = nparams)