################################# #
# maxent.jar ENMdetails object ####
################################# #

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
  msg <- paste("maxent.jar v.", v, "from dismo package v.", packageVersion('dismo'))
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
  if(!grepl("T", tune.tbl.i$fc)) out$args <- c(out$args, "nothreshold")
  out$args <- c(out$args, paste0("betamultiplier=", tune.tbl.i$rm, sep=""))
  out <- c(out, other.args)
  return(out)
}

auc <- function(occs.vals, bg.vals, mod, other.args, doClamp) {
  e <- dismo::evaluate(occs.vals, bg.vals, mod, args = c("outputformat=raw", ifelse(doClamp == TRUE, "doclamp=true", "doclamp=false")))@auc
  return(e)
}

pred <- function(mod, envs, other.args, doClamp) {
  pred <- dismo::predict(mod, envs, args = c("outputformat=raw", ifelse(doClamp == TRUE, "doclamp=true", "doclamp=false")), na.rm = TRUE)
  return(pred)
}

nparams <- function(mod) {
  lambdas <- mod@lambdas[1:(length(mod@lambdas)-4)]
  countNonZeroParams <- function(x) if(strsplit(x, split=", ")[[1]][2] != '0.0') 1
  np <- sum(unlist(sapply(lambdas, countNonZeroParams)))
  return(np)
}

enm.maxent.jar <- ENMdetails(fun = fun, pkgs = pkgs, msgs = msgs, 
                        args = args, auc = auc, pred = pred, 
                        nparams = nparams)