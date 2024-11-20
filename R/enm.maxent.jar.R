################################# #
# maxent.jar ENMdetails object ####
################################# #

maxent.jar.name <- "maxent.jar"

maxent.jar.fun <- predicts::MaxEnt

maxent.jar.errors <- function(occs, envs, bg, tune.args, partitions, algorithm, 
                          partition.settings, other.settings, 
                          categoricals, doClamp, clamp.directions) {
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
  # if(is.null(getOption('dismo_rJavaLoaded'))) {
  #   # to avoid trouble on macs
  #   Sys.setenv(NOAWT=TRUE)
  #   if ( requireNamespace('rJava') ) {
  #     rJava::.jpackage('predicts')
  #     options(dismo_rJavaLoaded=TRUE)
  #   } else {
  #     stop('rJava cannot be loaded.')
  #   }
  # }
}

maxent.jar.msgs <- function(tune.args, other.settings) {
  # rJava::.jpackage("predicts")
  mxe <- rJava::.jnew("meversion") 
  v <- try(rJava::.jcall(mxe, "S", "meversion"))
  msg <- paste0("maxent.jar v", v, " using the predicts package v", packageVersion('predicts'))
  return(msg)
}

maxent.jar.args <- function(occs.z, bg.z, tune.tbl.i, other.settings) {
  out <- list()
  out$x <- rbind(occs.z, bg.z)
  out$p <- c(rep(1, nrow(occs.z)), rep(0, nrow(bg.z)))
  out$args <- c("noremoveDuplicates", "noautofeature", other.settings$other.args)
  if(!grepl("L", tune.tbl.i$fc)) out$args <- c(out$args, "nolinear")
  if(!grepl("Q", tune.tbl.i$fc)) out$args <- c(out$args, "noquadratic")
  if(!grepl("H", tune.tbl.i$fc)) out$args <- c(out$args, "nohinge")
  if(!grepl("P", tune.tbl.i$fc)) out$args <- c(out$args, "noproduct")
  if(!grepl("T", tune.tbl.i$fc)) out$args <- c(out$args, "nothreshold") else out$args <- c(out$args, "threshold")
  out$args <- c(out$args, paste0("betamultiplier=", tune.tbl.i$rm, sep=""))
  out$path <- other.settings$path
  return(out)
}

#' @importFrom predicts predict
maxent.jar.predict <- function(mod, envs, other.settings) {
  output.format <- paste0("outputformat=", other.settings$pred.type)
  model.clamp <- ifelse(other.settings$doClamp == TRUE, "doclamp=true", "doclamp=false")
  pred <- predicts::predict(mod, envs, args = c(output.format, model.clamp))
  return(pred)
}

maxent.jar.ncoefs <- function(mod) {
  lambdas <- mod@lambdas[1:(length(mod@lambdas)-4)]
  countNonZeroParams <- function(x) if(strsplit(x, split=", ")[[1]][2] != '0.0') 1
  np <- sum(unlist(sapply(lambdas, countNonZeroParams)))
  return(np)
}

maxent.jar.variable.importance <- function(mod) {
  res <- mod@results
  # percent contribution is a heuristic measure of variable importance
  # quoting A Brief Tutorial on Maxent (Phillips 2017):
  # "(It) depend(s) on the particular path that the Maxent code uses to get to the optimal solution, 
  # and a different algorithm could get to the same solution via a different path, resulting in 
  # different percent contribution values. In addition, when there are highly correlated environmental 
  # variables, the percent contributions should be interpreted with caution."
  # ref: https://biodiversityinformatics.amnh.org/open_source/maxent/Maxent_tutorial2017.pdf
  pc <- res[grepl('contribution', rownames(res)),]
  # permutation importance is a non-heuristic measure of variable importance
  # From the older 2006 version of A Brief Tutorial:
  # "The permutation importance measure depends only on the final Maxent model, not the path used to 
  # obtain it. The contribution for each variable is determined by randomly permuting the values of 
  # that variable among the training points (both presence and background) and measuring the resulting 
  # decrease in training AUC. A large decrease indicates that the model depends heavily on that variable. 
  # Values are normalized to give percentages."
  pi <- res[grepl('permutation', rownames(res)),]
  varnames <- sapply(strsplit(names(pc), '.contribution'), function(x) x[1])
  df <- data.frame(variable=varnames, percent.contribution=pc, permutation.importance=pi, row.names=NULL)
  return(df)
}

#' @title ENMdetails maxent.jar
#' @description This is the ENMdetails implementation for maxent.jar, the Java version of
#' the Maxent algorithm. The configuration for running the model differs slightly from that
#' in previous versions of ENMeval (0.3.0 and before) in that this version (>=2.0.0) uses the
#' default of adding presences to the background for model training, while previous versions
#' had turned this off. Specifically, previous versions ran maxent() with "noaddsamplestobackground"
#' in the "args" vector argument, while this version does not.
#' @export
enm.maxent.jar <- ENMdetails(name = maxent.jar.name, fun = maxent.jar.fun, errors = maxent.jar.errors,
                             msgs = maxent.jar.msgs, args = maxent.jar.args,
                             predict = maxent.jar.predict, ncoefs = maxent.jar.ncoefs, variable.importance = maxent.jar.variable.importance)
