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

args <- function(occs.vals, bg.vals, tune.tbl.i, other.args) {
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
  out <- c(out, other.args)
  return(out)
}

aic <- function(occs, nparam, mod.full.pred.all) {
  calc.aicc(occs, nparam, mod.full.pred.all)
}

eval <- function(occs.vals, bg.vals, mod, other.args, doClamp) {
  clamp <- ifelse(doClamp == TRUE, "doclamp=true", "doclamp=false")
  e <- dismo::evaluate(occs.vals, bg.vals, mod, args = c("outputformat=cloglog", clamp))
  return(e)
}

kstats <- function(e.test, mod, other.args) {
  user.kstats <- c(maxTSS.test = max(e.test@TPR + e.test@TNR) - 1, maxKappa.test = max(e.test@kappa))
  return(user.kstats)
}

pred <- function(mod, envs, other.args, doClamp, pred.type) {
  type.arg <- paste("outputformat", pred.type, sep = "=")
  clamp.arg <- ifelse(doClamp == TRUE, "doclamp=true", "doclamp=false")
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
                        args = args, aic = aic, eval = eval, kstats = kstats, 
                        pred = pred, nparams = nparams)