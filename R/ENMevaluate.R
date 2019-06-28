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
#' @param user.enm ENMdetails object specified by the user; this model will be
#' used for the analysis, and is an alternative to specifying mod.name
#' @param tune.args named list of model settings to be tuned
#' @param other.args named list of any additional model arguments not specified 
#' for tuning
#' @param categoricals character vector of names of categorical 
#' environmental variables
#' @param partitions character of name of partitioning technique (see
#' \code{?partitions})
#' @param occ.grp numeric vector of partition group (fold) for each
#' occurrence locality, intended for user-defined partitions
#' @param bg.grp numeric vector of partition group (fold) for each background 
#' (or pseudo-absence) locality, intended for user-defined partitions
#' @param occs.ind matrix or data frame with two columns for longitude and latitude 
#' of occurrence localities, in that order, intended for independent evaluation;
#' when \code{partitions = "independent"}; these occurrences will be used only 
#' for evaluation, and not for model training, and thus no cross validation will 
#' be done
#' @param kfolds numeric for number of partition groups (grp), only for random
#' k-fold partitioning
#' @param aggregation.factor numeric vector with length 2 that specifies the
#' factors for aggregating \code{envs} in order to perform checkerboard
#' partitioning
#' @param n.bg numeric for number of random background (or pseudo-absence) points
#' to sample; necessary if \code{bg} is NULL
#' @param overlap boolean (TRUE or FALSE); if TRUE, calculate niche overlap 
#' statistics
#' @param overlapStat character; one or two (vector) niche overlap statistics:
#' choices are "D" and "I" -- see ?calc.niche.overlap for more details
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
#' @importFrom foreach %dopar%
#' @export 
#' 

ENMevaluate <- function(occs, envs = NULL, bg = NULL, occs.vals = NULL, bg.vals = NULL, 
                        tune.args = NULL, other.args = NULL, categoricals = NULL, mod.name,
                        user.enm = NULL,
                        partitions = NULL, occ.grp = NULL, bg.grp = NULL, occs.ind = NULL, 
                        kfolds = NA, aggregation.factor = c(2, 2), n.bg = 10000, overlap = FALSE, 
                        overlapStat = c("D", "I"), doClamp = TRUE, skipRasters = FALSE, 
                        abs.auc.diff = TRUE, parallel = FALSE, numCores = NULL, updateProgress = FALSE,
                        # legacy parameters
                        occ = NULL, env = NULL, bg.coords = NULL, RMvalues = NULL, fc = NULL,
                        algorithm = NULL, method = NULL, bin.output = NULL, rasterPreds = NULL,
                        clamp = NULL, progbar = NULL) {
  
  # legacy parameter handling so ENMevaluate doesn't break for older code
  all.legacy <- list(occ, env, bg.coords, RMvalues, fc, algorithm, method, bin.output, rasterPreds)
  if(sum(sapply(all.legacy, function(x) !is.null(x))) > 0) {
    message("Running ENMeval v1.0.0 with legacy parameters. These will be phased out in the next version.\n")
  }
  if(!is.null(occ)) occs <- occ
  if(!is.null(env)) envs <- env
  if(!is.null(bg.coords)) bg <- bg.coords
  if(!is.null(method)) partitions <- method
  if(!is.null(rasterPreds)) skipRasters <- rasterPreds
  if(!is.null(algorithm)) {
    mod.name <- algorithm
    tune.args <- list(fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
                      rm = seq(0.5, 4, 0.5))
  }
  if(!is.null(RMvalues)) tune.args$rm <- RMvalues
  if(!is.null(fc)) tune.args$fc <- fc
  
  if(is.null(mod.name)) {
    stop("Please select a model name.\n")
  }
  
  # record start time
  start.time <- proc.time()
  
  if(is.null(user.enm)) enm <- lookup.enm(mod.name)
  
  ########### #
  # CHECKS ####
  ########### #
  
  ## general parameter checks
  all.partitions <- c("jackknife", "randomkfold", "block", "checkerboard1", 
                      "checkerboard2", "user", "independent", "none")
  if(!(partitions %in% all.partitions)) {
    stop("Please enter an accepted partition method.\n")
  }
  if(is.null(tune.args) & mod.name != "bioclim") {
    stop("Please specify tuning.args.\n")
  }
  
  # print model-specific message
  msg <- enm@msgs(tune.args)
  message(paste("*** Running ENMevaluate using", msg, "***\n"))
  
  # make table for all tuning parameter combinations
  tune.tbl <- expand.grid(tune.args, stringsAsFactors = FALSE)
  
  ## data checks and formatting
  # if occs is combined occurrence and background with environmental
  # predictor values (SWD format)
  if(is.null(envs)) {
    if(is.null(occs.vals)) {
      stop("If inputting data without rasters (SWD), please specify both occs.vals and bg.vals.\n")
      occs.vals <- as.data.frame(occs.vals)
    }
    if(is.null(bg.vals)) {
      stop("If inputting data without rasters (SWD), please specify both occs.vals and bg.vals.\n")
      bg.vals <- as.data.frame(bg.vals)
    }
    if(is.null(bg)) {
      stop("If inputting data without rasters (SWD), please specify bg in addition to occs.\n")
    }
    warning("Data without rasters were input (SWD format), so no raster predictions will be generated and AICc cannot be calculated.\n")
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
    grp <- get.jackknife(occs, bg)
    parts.msg <- "Doing model evaluations with k-1 jackknife cross validation...\n"
  }
  if(partitions == "randomkfold") {
    grp <- get.randomkfold(occs, bg, kfolds)
    parts.msg <- paste0("Doing model evaluations with random", kfolds, "-fold cross validation...\n")
  }
  if(partitions == "block") {
    grp <- get.block(occs, bg)
    parts.msg <- "Doing model evaluations with spatial block (4-fold) cross validation...\n"
  }
  if(partitions == "checkerboard1") {
    if(is.null(envs)) stop("Cannot use checkerboard partitioning if envs is NULL.")
    grp <- get.checkerboard1(occs, envs, bg, aggregation.factor)
    parts.msg <- "Doing model evaluations with checkerboard (2-fold) cross validation...\n"
  }
  if(partitions == "checkerboard2") {
    if(is.null(envs)) stop("Cannot use checkerboard partitioning if envs is NULL.")
    grp <- get.checkerboard2(occs, envs, bg, aggregation.factor)
    parts.msg <- "Doing model evaluations with hierarchical checkerboard (4-fold) cross validation...\n"
  }
  if(partitions == "user") {
    grp <- list(occ.grp = occ.grp, bg.grp = bg.grp)
    userk <- length(unique(occ.grp))
    parts.msg <- paste0("Doing model evaluations with user-defined ", userk, "-fold cross validation...\n")
  }
  if(partitions == "independent") {
    if(is.null(occs.ind)) stop("Cannot use independent partitioning if occs.ind is NULL.")
    grp <- NULL
    parts.msg <- "Doing model evaluations with independent testing data...\n"
  }
  if(partitions == "none") {
    grp <- NULL
    parts.msg <- "Skipping model evaluations (only calculating AICc)...\n"
  }
  
  # unpack occ.grp and bg.grp
  occ.grp <- grp$occ.grp
  bg.grp <- grp$bg.grp
  
  ############# #
  # ANALYSIS ####
  ############# #
  
  # extract predictor variable values at coordinates for occs and bg
  if(!is.null(envs)) {
    occs.vals <- as.data.frame(raster::extract(envs, occs))
    bg.vals <- as.data.frame(raster::extract(envs, bg))  
  }
  
  # remove rows from occs, occs.vals, occ.grp with NA for any predictor variable
  occs.vals.na <- sum(rowSums(is.na(occs.vals)) > 0)
  bg.vals.na <- sum(rowSums(is.na(bg.vals)) > 0)
  if(occs.vals.na > 0) {
    i <- !apply(occs.vals, 1, anyNA)
    occs <- occs[i,]
    occ.grp <- occ.grp[i]
    occs.vals <- occs.vals[i,]
    warning(paste0("Occurrence records found (n = ", occs.vals.na, ") with NA for at least one predictor variable. Removed these from analysis, resulting in ", nrow(occs.vals), " occurrence records.\n"), immediate. = TRUE)
  }
  # do the same for associated bg variables
  if(bg.vals.na > 0) {
    i <- !apply(bg.vals, 1, anyNA)
    bg <- bg[i,]
    bg.grp <- bg.grp[i]
    bg.vals <- bg.vals[i,]
    warning(paste0("Background records found (n = ", bg.vals.na, ") with NA for at least one predictor variable. Removed these from analysis, resulting in ", nrow(bg.vals), " background records.\n"), immediate. = TRUE)
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
  
  message(parts.msg)
  
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
      results[[i]] <- cv.enm(occs.vals, bg.vals, occ.grp, bg.grp, envs, enm,
                             partitions, tune.tbl[i,], other.args, categoricals, 
                             occs.ind, doClamp, skipRasters, abs.auc.diff)
    }
    close(pb)
  }else{
    # parallel implementation
    # set up parallel processing functionality
    allCores <- parallel::detectCores()
    if (is.null(numCores)) {
      numCores <- allCores
    }
    cl <- parallel::makeCluster(numCores)
    doSNOW::registerDoSNOW(cl)
    numCoresUsed <- foreach::getDoParWorkers()
    message(paste0("Of ", allCores, " total cores using ", numCoresUsed, "..."))
    
    if(mod.name == "maxent.jar") pkgs <- c("dismo", "raster", "rJava")
    if(mod.name == "maxnet") pkgs <- c("dismo", "raster", "maxnet")
    if(mod.name == "brt") pkgs <- c("dismo", "raster", "gbm")
    if(mod.name == "bioclim") pkgs <- c("dismo", "raster")
    
    message("Running in parallel...")
    n <- ifelse(nrow(tune.tbl) > 0, nrow(tune.tbl), 1)
    pb <- txtProgressBar(0, n, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    results <- foreach::foreach(i = 1:n, .packages = pkgs, .options.snow = opts) %dopar% {

      cv.enm(occs.vals, bg.vals, occ.grp, bg.grp, envs, enm,
             partitions, tune.tbl[i,], other.args, categoricals, occs.ind, doClamp,
             skipRasters, abs.auc.diff)
    }
    close(pb)
    parallel::stopCluster(cl)
  }
  
  ################# #
  # collate results 
  ################# #
  
  if(nrow(tune.tbl) == 0) {
    # if not tuned settings, the "tune name" is the model name
    tune.names <- mod.name
  } else {
    # define tuned settings names and bind them to the tune table
    tune.names <- apply(tune.tbl, 1, function(x) paste(x, collapse = "_"))
    tune.tbl <- cbind(tune.tbl, tune.args = tune.names, stringsAsFactors = FALSE)
  }
  # gather all full models into list and name them
  mod.full.all <- lapply(results, function(x) x$mod.full)
  names(mod.full.all) <- tune.names
  # gather all training AUCs into vector
  auc.train.all <- sapply(results, function(x) x$train.AUC)
  # gather all statistics into a data frame
  kstats.all <- lapply(results, function(x) x$kstats)
  # gather all model prediction rasters into a stack and name them
  if(skipRasters == FALSE & !is.null(envs)) {
    mod.full.pred.all <- raster::stack(sapply(results, function(x) x$mod.full.pred))
    names(mod.full.pred.all) <- tune.names
  }else{
    mod.full.pred.all <- raster::stack()
  }
  
  # make data frame of stats (for all grp)
  kstats.df <- dplyr::bind_rows(kstats.all)
  
  # define number of grp (the value of "k") as number of
  # rows in one of the model runs
  nk <- nrow(kstats.all[[1]])
  
  if(partitions != "none") {
    # define number of rows in tune.tbl
    n <- ifelse(nrow(tune.tbl) > 0, nrow(tune.tbl), 1)
    # define number of evaluation statistics
    nstat <- ncol(kstats.all[[1]])
    # define number of settings (plus the tune.args field)
    ns <- ncol(tune.tbl)
    # add in columns for tuning settings
    if(nrow(tune.tbl) > 0) {
      kstats.df <- cbind(apply(tune.tbl, 2, rep, each = nk), kstats.df, stringsAsFactors = FALSE)
    }
    
    # get number of columns in kstats.df
    nc <- ncol(kstats.df)
    
    # if model settings for tuning were input, summarize by averaging all grp 
    # per model setting combination
    if(ns > 0) {
      kstats.df <- dplyr::group_by_at(kstats.df, 1:ns)
    }
    kstats.avg.df <- kstats.df %>% 
      dplyr::summarize(avg.test.AUC = mean(auc.test),
                       var.test.AUC = corrected.var(auc.test, nk),
                       min.test.AUC = min(auc.test),
                       avg.diff.AUC = mean(auc.diff),
                       var.diff.AUC = corrected.var(auc.diff, nk),
                       max.diff.AUC = max(auc.diff),
                       avg.test.orMTP = mean(or.mtp),
                       var.test.orMTP = var(or.mtp),
                       max.test.orMTP = max(or.mtp),
                       avg.test.or10pct = mean(or.10p),
                       var.test.or10pct = var(or.10p),
                       max.test.or10pct = max(or.10p)) %>%
      dplyr::ungroup()
    # reorder based on original order of tune names (summarize forces an alphanumeric reorder)
    kstats.avg.df <- kstats.avg.df[match(tune.names, kstats.avg.df$tune.args),]
    if(ns > 0) {
      stats.df <- cbind(kstats.avg.df[, 1:ns], auc.train = auc.train.all, 
                        kstats.avg.df[, seq(ns+1, ns+12)])  
    }else{
      stats.df <- cbind(auc.train = auc.train.all, kstats.avg.df) 
    }
  }else{
    tune.cols <- tune.tbl[order(tune.tbl[,1]),]
    stats.df <- cbind(tune.cols, auc.train = auc.train.all)
  }
  
  # calculate number of non-zero parameters in model
  nparams <- sapply(mod.full.all, enm@nparams)
  # calculate AICc for Maxent models
  if(mod.name %in% c("maxent.jar", "maxnet")) {
    stats.df <- cbind(stats.df, calc.aicc(nparams, occs, mod.full.pred.all))
  }else{
    warning(paste0("Not able to calculate AICc for ", mod.name, "... returning NAs.\n"))
  }
  stats.df$nparam <- nparams
  # stats.df <- tibble::as_tibble(stats.df)
  kstats.df <- as.data.frame(kstats.df)
  
  res <- list(stats = stats.df, kstats = kstats.df, mods = mod.full.all,
              preds = mod.full.pred.all)
  
  if(is.null(occ.grp)) occ.grp <- 0
  if(is.null(bg.grp)) bg.grp <- 0
  e <- ENMevaluation(algorithm = mod.name, tune.settings = tune.tbl,
                     results = res$stats, results.grp = res$kstats,
                     predictions = res$preds, models = res$mods, 
                     partition.method = partitions,
                     occ.pts = occs, occ.grp = occ.grp,
                     bg.pts = bg, bg.grp = bg.grp)
  
  # if niche overlap selected, calculate and add the resulting matrix to results
  nr <- raster::nlayers(e@predictions)
  if(overlap == TRUE) {
    if(nr == 0) {
      warning("Cannot calculate niche overlap without model prediction rasters.\n")
    }else if(nr == 1) {
      warning("Only 1 model prediction raster found. Need at least 2 rasters to calculate niche overlap. Increase number of tuning arguments and run again.\n") 
    }else{
      for(ovStat in overlapStat) {
        message(paste0("Calculating niche overlap for statistic ", ovStat, "...\n"))
        overlap.mat <- calc.niche.overlap(e@predictions, ovStat)
        e@overlap[[ovStat]] <- overlap.mat
      }
    }
  }
  
  # calculate time expended and print message
  timed <- proc.time() - start.time
  t.min <- floor(timed[3] / 60)
  t.sec <- timed[3] - (t.min * 60)
  message(paste("ENMevaluate completed in", t.min, "minutes", round(t.sec, 1), "seconds."))
  
  return(e)
}
