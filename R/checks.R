checks.partitions <- function(partitions, partition.settings, occs.testing) {
  if(!(partitions %in% c("jackknife", "randomkfold", "block", "checkerboard", 
                         "user", "testing", "none"))) {
    stop("Please enter an accepted partition method. See ?partitions for details.")
  }
  
  if(partitions == "testing" & is.null(occs.testing)) {
    stop("If doing testing evaluations, please provide testing data (occs.testing).")
  }
  
  if(partitions == "block") {
    if(is.null(partition.settings$orientation)) {
      stop("For block partitioning, please provide an orientation in partition.settings. See ?partitions for details.")
    }
  }
  
  if(partitions == "checkerboard") {
    if(is.null(envs)) {
      stop('For checkerboard partitioning, predictor variable rasters (the envs argument) are required.')
    }
    if(is.null(partition.settings$aggregation.factor)) {
      stop("For checkerboard partitioning, please provide an aggregation factor in partition.settings. See ?partitions for details.")
    }
  }
  
  if(partitions == "randomkfold") {
    if(is.null(partition.settings$kfolds)) {
      stop('For random k-fold partitioning, please provide a value for kfolds in partition.settings. See ?partitions for details.')  
    }else{
      if(partition.settings$kfolds == 0) {
        stop('For random k-fold partitioning, please provide a numeric, non-zero value of kfolds.')  
      }
    }
  }
  
  # if occs.testing input, coerce partitions to 'testing'
  if(partitions == "testing") {
    if(is.null(occs.testing)) {
      stop('If performing fully withheld testing, enter occs.testing dataset and assign partitions to "testing".')
    }
    if(nrow(occs.testing) == 0) {
      stop('If performing fully withheld testing, enter occs.testing dataset and assign partitions to "testing".')
    }
  }
}

checks.other.settings <- function(other.settings) {
  if(is.null(other.settings$removeduplicates)) other.settings$removeduplicates = TRUE
  if(is.null(other.settings$abs.auc.diff)) other.settings$abs.auc.diff <- TRUE
  if(is.null(other.settings$pred.type)) other.settings$pred.type <- "cloglog"
  if(is.null(other.settings$validation.bg)) other.settings$validation.bg <- "full"
  # add whether to use ecospat to other.settings to avoid multiple requires
  other.settings <- c(other.settings, ecospat.use = ecospat.use)
  return(other.settings)
}

checks.envs <- function(envs, algorithm) {
  if(!is.null(envs)) {
    # environmental raster data checks
    if(inherits(envs, "SpatRaster") == FALSE) {
      stop('From this version of ENMeval, the package will only use "terra" raster data types. Please convert from "raster" to "terra" with terra::rast(r), where r is a RasterStack.')
    }else{
      if(terra::nlyr(envs) < 2 & algorithm %in% c("maxent.jar", "maxnet")) {
        stop('Maxent is generally not designed to be run with a single predictor variable. Please rerun with multiple predictors.')
      }
      if(terra::nlyr(envs) < 2 & algorithm == "bioclim") {
        stop('BIOCLIM is not designed to be run with a single predictor variable. Please rerun with multiple predictors.')
      }
    }
  }
}

checks.tuning <- function(tune.args) {
  # if a vector of tuning arguments is numeric, make sure it is sorted 
  # (for results table and plotting)
  tune.args.num <- which((sapply(tune.args, class) %in% c("numeric", "integer")) 
                         & sapply(tune.args, length) > 1)
  for(i in tune.args.num) {
    tune.args[[i]] <- sort(tune.args[[i]])
  }
  return(tune.args)
}

checks.cats <- function(occs, envs, categoricals) {
  # find factor rasters or columns and identify them as categoricals
  if(!is.null(envs)) {
    cat.extra <- unique(c(categoricals, names(envs)[which(terra::is.factor(envs))]))
  }else{
    cat.extra <- unique(c(categoricals, names(occs)[which(sapply(occs, is.factor))]))
  }
  if(length(cat.extra) > categoricals) {
    stop(paste0("There are ", length(cat.extra), " variables with class factor, but only ", length(categoricals), " categorical variables specified. Please enter all categorical variable names in argument 'categoricals'."))
  }
}