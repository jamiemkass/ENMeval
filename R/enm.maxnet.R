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

args <- function(occs.vals, bg.vals, tune.tbl.i, other.args) {
  out <- list()
  out$data <- rbind(occs.vals, bg.vals)
  out$p <- c(rep(1, nrow(occs.vals)), rep(0, nrow(bg.vals)))
  out$f <- maxnet::maxnet.formula(out$p, out$data, classes = tolower(tune.tbl.i$fc))
  out$regmult <- tune.tbl.i$rm
  out <- c(out, other.args)
  return(out)
}

aic <- function(occs, nparams, mod.full.pred.all) {
  calc.aicc(occs, nparams, mod.full.pred.all)
}

eval <- function(occs.vals, bg.vals, mod, other.args, doClamp) {
  e <- dismo::evaluate(occs.vals, bg.vals, mod, clamp = doClamp, type = "cloglog")
  return(e)
}

kstats <- function(e.test, mod, other.args) {
  user.kstats <- c(maxTSS.test = max(e.test@TPR + e.test@TNR) - 1, maxKappa.test = max(e.test@kappa))
  return(user.kstats)
}

pred <- function(mod, envs, other.args, doClamp, pred.type) {
  pred <- maxnet.predictRaster(mod, envs, doClamp, type = pred.type, other.args)
  return(pred)
}

nparams <- function(mod) {
  length(mod$betas)
}

#' @export
enm.maxnet <- ENMdetails(name = name, fun = fun, pkgs = pkgs, msgs = msgs, 
                        args = args, aic = aic, eval = eval, kstats = kstats, 
                        pred = pred, nparams = nparams)