#' @title Tuning and evaluation of ecological niche models
#' @description \code{ENMevaluate()} builds ecological niche models iteratively across a range of user-specified tuning settings. Users can choose to tune with cross validation or an independent occurrence dataset. \code{ENMevaluate()} returns an \code{ENMevaluation} object with slots containing evaluation statistics (for each combination of settings and for each cross validation fold therein) and raster predictions for each model. The evaluation statistics are meant to aid users in identifying settings that balance model fit and predictive ability.
#' @param occs matrix or data frame with two columns for longitude and latitude 
#' of occurrence localities, in that order
#' @param envs Raster* object of environmental variables (must be in 
#' same geographic projection as occurrence data)
#' @param bg matrix or data frame with two columns for longitude and latitude of 
#' background (or pseudo-absence) localities, in that order; if NULL, points 
#' will be randomly sampled across \code{envs} with the number specified by 
#' parameter \code{n.bg}
#' @param occs.vals matrix or data frame of environmental values corresponding
#' to occurrence localities, intended to be input when environmental rasters
#' are not used (\code{envs} is NULL) 
#' @param bg.vals matrix or data frame of environmental values corresponding
#' to background (or pseudo-absence) localities, intended to be input when 
#' environmental rasters are not used (\code{envs} is NULL) 
#' @param mod.name character of the name of the chosen model
#' @param tune.args named list of model settings to be tuned
#' @param other.args named list of any additional model arguments not specified 
#' for tuning
#' @param categoricals character vector of names of categorical 
#' environmental variables
#' @param partitions character of name of partitioning technique (see
#' \code{?partitions})
#' @param occs.folds numeric vector of partition group (fold) for each
#' occurrence locality, intended for user-defined partitions
#' @param bg.folds numeric vector of partition group (fold) for each background 
#' (or pseudo-absence) locality, intended for user-defined partitions
#' @param occs.ind matrix or data frame with two columns for longitude and latitude 
#' of occurrence localities, in that order, intended for independent evaluation;
#' when \code{partitions = "independent"}; these occurrences will be used only 
#' for evaluation, and not for model training, and thus no cross validation will 
#' be done
#' @param kfolds numeric for number of partition groups (folds), only for random
#' k-fold partitioning
#' @param aggregation.factor numeric vector with length 2 that specifies the
#' factors for aggregating \code{envs} in order to perform checkerboard
#' partitioning
#' @param n.bg numeric for number of random background (or pseudo-absence) points
#' to sample; necessary if \code{bg} is NULL
#' @param overlap boolean (TRUE or FALSE); if TRUE, calculate niche overlap 
#' statistics
#' @param doClamp boolean (TRUE or FALSE); if TRUE, clamp model responses; only
#' applicable for Maxent models
#' @param skipRasters boolean (TRUE or FALSE); if TRUE, skip raster predictions
#' @param abs.auc.diff boolean (TRUE or FALSE); if TRUE, take absolute value of
#' AUCdiff; default is TRUE
#' @param parallel boolean (TRUE or FALSE); if TRUE, run with parallel processing
#' @param numCores boolean (TRUE or FALSE); if TRUE, use specifed number of cores
#' for parallel processing
#' @param updateProgress boolean (TRUE or FALSE); if TRUE, use shiny progress
#' bar; only for use in shiny apps
#'
#' @return 
#'
#' @examples
#'
#' @importFrom magrittr %>%
#' @export 

ENMevaluate <- function(occs, envs = NULL, bg = NULL, occs.vals = NULL, bg.vals = NULL, 
                        mod.name, tune.args = NULL, other.args = NULL, categoricals = NULL, 
                        partitions = NULL, occs.folds = NULL, bg.folds = NULL, occs.ind = NULL, 
                        kfolds = NA, aggregation.factor = c(2, 2), n.bg = 10000, overlap = FALSE,   
                        doClamp = TRUE, skipRasters = FALSE, abs.auc.diff = TRUE, parallel = FALSE, 
                        numCores = NULL, updateProgress = FALSE) {
  
  # record start time
  start.time <- proc.time()
  # get model function from input name
  mod.fun <- get.mod.fun(mod.name)
  
  ########### #
  # CHECKS ####
  ########### #
  
  ## general parameter checks
  if(!(partitions %in% c("jackknife", "randomkfold", "block", "checkerboard1", 
                         "checkerboard2", "user", "independent", "none"))) {
    stop("Please make sure partition method is one of the available options.")
  }
  if(is.null(tune.args) & mod.name != "bioclim") {
    stop("Must specify tuning.args for all models except BIOCLIM.")
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
      stop("If inputting data without rasters (SWD), please specify both 
           occs.vals and bg.vals.")
    }
    occs.vals <- as.data.frame(occs.vals)
    if(is.null(bg.vals)) {
      stop("If inputting data without rasters (SWD), please specify both 
           occs.vals and bg.vals.")
    }
    bg.vals <- as.data.frame(bg.vals)
    if(is.null(bg)) {
      stop("If inputting data without rasters (SWD), please specify bg in 
           addition to occs.")
    }
    message("Data without rasters were input (SWD format), so no raster 
            predictions will be generated and AICc cannot be calculated.")
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
  # for occs.ind settings, partitions should be NULL
  if(partitions == "jackknife") {
    folds <- get.jackknife(occs, bg)
    message("Doing model evaluations with k-1 jackknife cross validation...")
  }
  if(partitions == "randomkfold") {
    folds <- get.randomkfold(occs, bg, kfolds)
    message(paste0("Doing model evaluations with random", kfolds, "-fold cross validation..."))
  }
  if(partitions == "block") {
    folds <- get.block(occs, bg)
    message("Doing model evaluations with spatial block (4-fold) cross validation...")
  }
  if(partitions == "checkerboard1") {
    if(is.null(envs)) stop("Cannot use checkerboard partitioning if envs is NULL.")
    folds <- get.checkerboard1(occs, envs, bg, aggregation.factor)
    message("Doing model evaluations with checkerboard (2-fold) cross validation...")
  }
  if(partitions == "checkerboard2") {
    if(is.null(envs)) stop("Cannot use checkerboard partitioning if envs is NULL.")
    folds <- get.checkerboard2(occs, envs, bg, aggregation.factor)
    message("Doing model evaluations with hierarchical checkerboard (4-fold) cross validation...")
  }
  if(partitions == "user") {
    folds <- list(occs.folds = occs.folds, bg.folds = bg.folds)
    userk <- length(unique(occs.folds))
    message(paste0("Doing model evaluations with user-defined ", userk, "-fold cross validation..."))
  }
  if(partitions == "independent") {
    if(is.null(occs.ind)) stop("Cannot use independent partitioning if occs.ind is NULL.")
    folds <- NULL
    message("Doing model evaluations with independent testing data...")
  }
  if(partitions == "none") {
    folds <- NULL
    message("Skipping model evaluations (only calculating AICc)...")
  }
  
  # unpack occs.folds and bg.folds
  occs.folds <- folds$occs.folds
  bg.folds <- folds$bg.folds
  
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
  
  ################ #
  # tuning 
  ################ #
  
  if(parallel == FALSE) {
    # regular implementation
    results <- list()
    n <- ifelse(nrow(tune.tbl) > 0, nrow(tune.tbl), 1)
    
    # set up the console progress bar
    pb <- txtProgressBar(0, n, style = 3)
    
    for(i in 1:n) {
      # and (optionally) the shiny progress bar (updateProgress)
      if(n > 1) {
        if(is.function(updateProgress)) {
          text <- paste0('Running ', paste(as.character(tune.tbl[i,]), collapse = ""), '...')
          updateProgress(detail = text)
        }
        setTxtProgressBar(pb, i)
      }
      results[[i]] <- cv.enm(occs.vals, bg.vals, occs.folds, bg.folds, envs, mod.fun, mod.name, 
                             partitions, tune.tbl[i,], other.args, categoricals, occs.ind, doClamp, 
                             skipRasters, abs.auc.diff)
    }
    close(pb)
  }else{
    # parallel implementation
    # set up parallel processing functionality
    allCores <- detectCores()
    if (is.null(numCores)) {
      numCores <- allCores
    }
    c1 <- makeCluster(numCores)
    registerDoParallel(c1)
    numCoresUsed <- getDoParWorkers()
    message(paste("Of", allCores, "total cores using", numCoresUsed))
    
    if(mod.name == "maxent.jar") pkgs <- c("dismo", "raster", "ENMeval", "rJava")
    if(mod.name == "maxnet") pkgs <- c("dismo", "raster", "ENMeval", "maxnet")
    if(mod.name == "brt") pkgs <- c("dismo", "raster", "ENMeval", "gbm")
    if(mod.name == "bioclim") pkgs <- c("dismo", "raster")
    
    message("Running in parallel...")
    n <- nrow(tune.tbl)
    results <- foreach(i = 1:n, .packages = pkgs) %dopar% {
      cv.enm(occs.vals, bg.vals, occs.folds, bg.folds, envs, mod.fun, mod.name, partitions, 
             tune.tbl[i,], other.args, categoricals, occs.ind, doClamp, skipRasters, abs.auc.diff)
    }
    stopCluster(c1)
  }
  
  ################# #
  # collate results 
  ################# #
  
  # gather all full models into list
  mod.full.all <- lapply(results, function(x) x$mod.full)
  # gather all training AUCs into vector
  auc.train.all <- sapply(results, function(x) x$auc.train)
  # gather all statistics into a data frame
  kstats.all <- lapply(results, function(x) x$kstats)
  if(skipRasters == FALSE & !is.null(envs)) {
    mod.full.pred.all <- raster::stack(sapply(results, function(x) x$mod.full.pred))
  }else{
    mod.full.pred.all <- raster::stack()
  }
  
  # make data frame of stats (for all folds)
  kstats.df <- dplyr::bind_rows(kstats.all)
  
  # define number of folds (the value of "k") as number of
  # rows in one of the model runs
  nk <- nrow(kstats.all[[1]])
  
  if(partitions != "none") {
    # define number of rows in tune.tbl
    n <- ifelse(nrow(tune.tbl) > 0, nrow(tune.tbl), 1)
    # define number of evaluation statistics
    nstat <- ncol(kstats.all[[1]])
    # define number of settings
    ns <- ncol(tune.tbl)
    # add in columns for tuning settings
    if(nrow(tune.tbl) > 0) {
      for(i in 1:ns) {
        kstats.df <- cbind(rep(tune.tbl[,i], each = nk), kstats.df)
        names(kstats.df)[1] <- names(tune.tbl)[i]
      }  
    }
    
    # make new column for fold number
    # kstats.df <- tibble::as_tibble(kstats.df)
    
    # get number of columns in kstats.df
    nc <- ncol(kstats.df)
    
    # if model settings for tuning were input, summarize by averaging all folds 
    # per model setting combination
    if(ns > 0) {
      kstats.df <- dplyr::group_by_at(kstats.df, 1:ns)
    }
    kstats.avg.df <- kstats.df %>% 
      dplyr::summarize(auc.test.mean = mean(auc.test),
                       auc.test.var = corrected.var(auc.test, nk),
                       auc.test.min = min(auc.test),
                       auc.diff.mean = mean(auc.diff),
                       auc.diff.var = corrected.var(auc.diff, nk),
                       auc.diff.max = max(auc.diff),
                       or.mtp.mean = mean(or.mtp),
                       or.mtp.var = var(or.mtp),
                       or.mtp.max = max(or.mtp),
                       or.10p.mean = mean(or.10p),
                       or.10p.var = var(or.10p),
                       or.10p.max = max(or.10p),
                       or.mss.mean = mean(or.mss),
                       or.mss.var = var(or.mss),
                       or.mss.max = max(or.mss)) %>%
      dplyr::ungroup()
    if(ns > 0) {
      stats.df <- cbind(kstats.avg.df[, 1:ns], auc.train = auc.train.all, 
                        kstats.avg.df[, seq(ns+1, ns+15)])  
    }else{
      stats.df <- cbind(auc.train = auc.train.all, kstats.avg.df) 
    }
  }else{
    tune.cols <- tune.tbl[order(tune.tbl[,1]),]
    stats.df <- cbind(tune.cols, auc.train = auc.train.all)
  }
  
  # calculate number of non-zero parameters in model
  nparams <- sapply(mod.full.all, function(x) no.params(x, mod.name))
  # calculate AICc for Maxent models
  if(mod.name %in% c("maxent.jar", "maxnet")) {
    stats.df <- cbind(stats.df, calc.aicc(nparams, occs, mod.full.pred.all))
  }else{
    warning(paste0("Not able to calculate AICc for ", mod.name, 
                   "... returning NAs"))
  }
  stats.df$nparam <- nparams
  stats.df <- tibble::as_tibble(stats.df)
  kstats.df <- tibble::as_tibble(kstats.df)
  
  res <- list(stats = stats.df, kstats = kstats.df, mods = mod.full.all,
              preds = mod.full.pred.all)
  
  if(is.null(occs.folds)) occs.folds <- 0
  if(is.null(bg.folds)) bg.folds <- 0
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
