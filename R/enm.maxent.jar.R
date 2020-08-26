################################# #
# maxent.jar ENMdetails object ####
################################# #

name <- "maxent.jar"

fun <- dismo::maxent

pkgs <- c("dismo", "raster", "rJava")

msgs <- function(tune.args) {
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

evaluate <- function(occs.z, bg.z, mod, other.settings) {
  output.format <- paste0("outputformat=", other.settings$pred.type)
  clamp <- ifelse(other.settings$clamp == TRUE, "doclamp=true", "doclamp=false")
  dismo::evaluate(occs.z, bg.z, mod, args = c(output.format, clamp))
}

predict <- function(mod, envs, other.settings) {
  output.format <- paste0("outputformat=", other.settings$pred.type)
  clamp <- ifelse(other.settings$clamp == TRUE, "doclamp=true", "doclamp=false")
  pred <- dismo::predict(mod, envs, args = c(output.format, clamp), na.rm = TRUE)
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
                             evaluate = evaluate, predict = predict, nparams = nparams)
