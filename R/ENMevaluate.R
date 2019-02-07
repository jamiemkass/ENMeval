#' @importFrom magrittr %>%
#' @export 

ENMevaluate <- function(occs, envs = NULL, bg = NULL, occs.vals = NULL, bg.vals = NULL, 
                        mod.fun, tune.args, other.args = NULL, categoricals = NULL, 
                        partitions = NULL, occs.folds = NULL, bg.folds = NULL, n.bg = 10000, 
                        overlap = FALSE, aggregation.factor = c(2, 2), kfolds = NA, bin.output = FALSE,
                        doClamp = TRUE, skipRasters = FALSE, abs.auc.diff = TRUE, parallel = FALSE, 
                        numCores = NULL, updateProgress = FALSE) {

  # record start time
  start.time <- proc.time()
  # get model function's name
  mod.name <- as.character(substitute(mod.fun))[3]
  print(mod.name)
  
  ########### #
  # CHECKS ####
  ########### #
  
  ## general parameter checks
  if(is.null(partitions)) {
    stop("Please specify a partition method for cross validation.")
  }
  if(!(partitions %in% c("jackknife", "randomkfold", "block", "checkerboard1", "checkerboard2", "user"))) {
    stop("Please make sure partition method is one of the available options.")
  }
  
  # print model-specific message
  model.msgs(tune.args, mod.name)
  
  # make table for all tuning parameter combinations
  tune.tbl <- expand.grid(tune.args, stringsAsFactors = FALSE)

  ## data checks and formatting
  # if occs is combined occurrence and background with environmental
  # predictor values (SWD format)
  if(is.null(envs)) {
    if(is.null(occs.vals)) {
      stop("If inputting data without rasters (SWD), please specify both occs.vals and bg.vals.")
    }
    occs.vals <- as.data.frame(occs.vals)
    if(is.null(bg.vals)) {
      stop("If inputting data without rasters (SWD), please specify both occs.vals and bg.vals.")
    }
    bg.vals <- as.data.frame(bg.vals)
    if(is.null(bg)) {
      stop("If inputting data without rasters (SWD), please specify bg in addition to occs.")
    }
    message("Data without rasters were input (SWD format), so no raster predictions will be generated and AICc cannot be calculated.")
  }else{
    # if no background points specified, generate random ones
    if(is.null(bg)) {
      bg <- dismo::randomPoints(envs, n = n.bg)
    }  
  }
  
  # make sure occs and bg are data frames with identical column names
  occs <- data.frame(longitude = occs[,1], latitude = occs[,2])
  bg <- data.frame(longitude = bg[,1], latitude = bg[,2])
  
  # print message for selected cross-validation method
  if(partitions == "jackknife") {
    ptn.msg <- "k-1 jackknife"
    folds <- get.jackknife(occs, bg)
  }
  if(partitions == "randomkfold") {
    ptn.msg <- paste0("random ", kfolds, "-fold")
    folds <- get.randomkfold(occs, bg, kfolds)
  }
  if(partitions == "block") {
    ptn.msg <- "spatial block (4-fold)"
    folds <- get.block(occs, bg)
  }
  if(partitions == "checkerboard1") {
    if(is.null(envs)) stop("Cannot use checkerboard partitioning is envs in NULL.")
    ptn.msg <- "checkerboard (2-fold)"
    folds <- get.checkerboard1(occs, envs, bg, aggregation.factor)
  }
  if(partitions == "checkerboard2") {
    if(is.null(envs)) stop("Cannot use checkerboard partitioning is envs in NULL.")
    ptn.msg <- "hierarchical checkerboard (4-fold)"
    folds <- get.checkerboard2(occs, envs, bg, aggregation.factor)
  }
  if(partitions == "user") {
    ptn.msg <- "user-defined k-fold"
    folds <- list(occs.folds = occs.folds, bg.folds = bg.folds)
  }
  
  # unpack occs.folds and bg.folds
  occs.folds <- folds$occs.folds
  bg.folds <- folds$bg.folds
  
  message(paste("Doing evaluations with", ptn.msg, "cross validation..."))
  
  ############# #
  # ANALYSIS ####
  ############# #
  
  # extract predictor variable values at coordinates for occs and bg
  if(!is.null(envs)) {
    occs.vals <- as.data.frame(raster::extract(envs, occs))
    bg.vals <- as.data.frame(raster::extract(envs, bg))  
  }
  
  # remove rows from occs, occs.vals, occs.folds with NA for any predictor variable
  occs.vals.na <- sum(rowSums(is.na(occs.vals)) > 0)
  bg.vals.na <- sum(rowSums(is.na(bg.vals)) > 0)
  if(occs.vals.na > 0) {
    i <- !apply(occs.vals, 1, anyNA)
    occs <- occs[i,]
    occs.folds <- occs.folds[i]
    occs.vals <- occs.vals[i,]
    message(paste("There were", occs.vals.na, "occurrence records with NA for at least one predictor variable. Removed these from analysis, resulting in", nrow(occs.vals), "occurrence records."))
  }
  # do the same for associated bg variables
  if(bg.vals.na > 0) {
    i <- !apply(bg.vals, 1, anyNA)
    occs <- occs[i,]
    bg.folds <- bg.folds[i]
    bg.vals <- bg.vals[i,]
    message(paste("There were", bg.vals.na, "background records with NA for at least one predictor variable. Removed these from analysis, resulting in", nrow(bg.vals), "background records."))
  }
  
  # convert fields for categorical data to factor class
  if(!is.null(categoricals)) {
    for(i in 1:length(categoricals)) {
      occs.vals[, categoricals[i]] <- as.factor(occs.vals[, categoricals[i]])
      bg.vals[, categoricals[i]] <- as.factor(bg.vals[, categoricals[i]])
    }
  }
  
  if(parallel == TRUE) {
    results <- tune.enm.parallel(occs.vals, bg.vals, occs.folds, bg.folds, envs, 
                                 mod.fun, mod.name, tune.tbl, other.args, categoricals, 
                                 doClamp, skipRasters, abs.auc.diff)
  } 
  else {
    results <- tune.enm(occs.vals, bg.vals, occs.folds, bg.folds, envs, 
                        mod.fun, mod.name, tune.tbl, other.args, categoricals, 
                        doClamp, skipRasters, abs.auc.diff, updateProgress)
  }
  
  res <- collateResults(results, tune.tbl, envs, mod.name, skipRasters)
  e <- ENMevaluation(algorithm = mod.name, 
                     stats = res$stats, kstats = res$kstats,
                     predictions = res$preds, models = res$mods, 
                     partition.method = partitions,
                     occs = occs, occs.folds = occs.folds,
                     bg = bg, bg.folds = bg.folds)
  
  
  # if niche overlap selected, calculate and add the resulting matrix to results
  if (overlap == TRUE) {
    if (length(args) > 1) {
      if(nlayers(results@predictions) > 1) {
        message("Calculating niche overlap")
        overlap.mat <- calc.niche.overlap(results@predictions, "D")
        results@overlap <- overlap.mat
      }
      else {
        message("Cannot calculate niche overlap without rasterPreds")
      }
    }
    else {
      message("Need >1 settings to do niche overlap")
    }
  }
  
  # calculate time expended and print message
  timed <- proc.time() - start.time
  t.min <- floor(timed[3] / 60)
  t.sec <- timed[3] - (t.min * 60)
  message(paste("ENMeval completed in", t.min, "minutes", round(t.sec, 1), "seconds."))
  
  return(e)
}
