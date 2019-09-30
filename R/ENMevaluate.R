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
#' @param numCores numeric for number of cores to use for parallel processing
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
  
  # if occs is combined occurrence and background with environmental
  # predictor values (SWD format)
  if(!is.null(envs)) {
    # make sure envs is a RasterStack -- if RasterLayer, maxent.jar crashes
    envs <- raster::stack(envs)
    # if no background points specified, generate random ones
    if(is.null(bg)) bg <- dismo::randomPoints(envs, n = n.bg)
    # extract predictor variable values at coordinates for occs and bg
    occs.vals <- raster::extract(envs, occs)
    bg.vals <- raster::extract(envs, bg)
    # if envs is NULL and values are specified (SWD), make sure they are data frames
  }else{
    warning("Data without rasters were input (SWD format), so no raster predictions will be generated and AICc cannot be calculated.\n", immediate. = TRUE)
    if(is.null(occs.vals)) {
      stop("If inputting data without rasters (SWD), please specify both occs.vals and bg.vals.\n")
    }
    if(is.null(bg.vals)) {
      stop("If inputting data without rasters (SWD), please specify both occs.vals and bg.vals.\n")
    }
  }
  # make sure values are data frames
  occs.vals <- as.data.frame(occs.vals)
  bg.vals <- as.data.frame(bg.vals)
  if(ncol(occs.vals) == 1) names(occs.vals) <- names(envs)
  if(ncol(bg.vals) == 1) names(bg.vals) <- names(envs)
  
  # make sure occs and bg are data frames with identical column names
  if(all(names(occs) != names(bg))) {
    warning('Datasets "occs" and "bg" have different column names. Assuming the first column for each is longitude and the second is latitude...\n', immediate. = TRUE)
    occs <- data.frame(x = occs[,1], y = occs[,2])
    bg <- data.frame(x = bg[,1], y = bg[,2])
  }
  
  # partition occs based on selected partition method
  # for occs.ind settings, partitions should be NULL
  grps <- switch(partitions, 
                 jackknife = get.jackknife(occs, bg),
                 randomkfold = get.randomkfold(occs, bg, kfolds),
                 block = get.block(occs, bg),
                 checkerboard1 = get.checkerboard1(occs, envs, bg, aggregation.factor),
                 checkerboard2 = get.checkerboard2(occs, envs, bg, aggregation.factor),
                 user = list(occ.grp = occ.grp, bg.grp = bg.grp),
                 independent = NULL,
                 none = NULL)
  parts.msg <- switch(partitions,
                      jackknife = "Doing model evaluations with k-1 jackknife cross validation...\n",
                      randomkfold = paste0("Doing model evaluations with random ", kfolds, "-fold cross validation...\n"),
                      block = "Doing model evaluations with spatial block (4-fold) cross validation...\n",
                      checkerboard1 = "Doing model evaluations with checkerboard (2-fold) cross validation...\n",
                      checkerboard2 = "Doing model evaluations with hierarchical checkerboard (4-fold) cross validation...\n",
                      user = paste0("Doing model evaluations with user-defined ", length(unique(occ.grp)), "-fold cross validation...\n"),
                      independent = "Doing model evaluations with independent testing data...\n",
                      none = "Skipping model evaluations (only calculating AICc)...\n")
  message(parts.msg)
  
  # unpack partition values
  occs.grp <- grps$occ.grp
  bg.grp <- grps$bg.grp
  
  # remove rows from occs.vals and occ.grp with NA for any predictor variable
  occs.vals.naRows <- which(rowSums(is.na(occs.vals)) > 0)
  occs.num.NA <- length(occs.vals.naRows)
  if(occs.num.NA > 0) {
    warning(paste0("Occurrence records found (n = ", occs.num.NA, ") with NA for at least one predictor variable. Removing these from analysis...\n"), immediate. = TRUE)
    occs.vals <- occs.vals[-occs.vals.naRows,]
    occs.grp <- occs.grp[-occs.vals.naRows]
  }
  # remove rows from bg.vals and bg.grp with NA for any predictor variable
  bg.vals.naRows <- which(rowSums(is.na(bg.vals)) > 0)
  bg.num.NA <- length(bg.vals.naRows)
  if(bg.num.NA > 0) {
    warning(paste0("Background records found (n = ", bg.num.NA, ") with NA for at least one predictor variable. Removing these from analysis...\n"), immediate. = TRUE)
    bg.vals <- bg.vals[-bg.vals.naRows,, drop = FALSE]
    bg.grp <- bg.grp[-bg.vals.naRows]
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
  
  # choose a built-in ENMdetails object matching the input model name
  # unless the model is chosen by the user
  if(is.null(user.enm)) enm <- lookup.enm(mod.name)
  # print model-specific message
  msg <- enm@msgs(tune.args)
  message(paste("*** Running ENMeval v1.0.0 using", msg, "***\n"))
  
  # make table for all tuning parameter combinations
  tune.tbl <- expand.grid(tune.args, stringsAsFactors = FALSE)
  
  if(parallel) {
    results <- tune.parallel(occs.vals, bg.vals, occs.grp, bg.grp, envs, enm, 
                             partitions, tune.tbl, other.args, categoricals, 
                             occs.ind, doClamp, skipRasters, abs.auc.diff, numCores)  
  }else{
    results <- tune.regular(occs.vals, bg.vals, occs.grp, bg.grp, envs, enm, 
                            partitions, tune.tbl, other.args, categoricals, 
                            occs.ind, doClamp, skipRasters, abs.auc.diff, updateProgress)
  }
  
  ################# #
  # results 
  ################# #
  
  if(nrow(tune.tbl) == 0) {
    # if not tuned settings, the "tune name" is the model name
    tune.names <- mod.name
  }else{
    # define tuned settings names and bind them to the tune table
    tune.names <- apply(tune.tbl, 1, function(x) paste(x, collapse = "_"))
    tune.tbl <- dplyr::bind_cols(tune.tbl, tune.args = tune.names)
  }
  # gather all full models into list and name them
  mod.full.all <- lapply(results, function(x) x$mod.full)
  names(mod.full.all) <- tune.names
  # gather all training AUCs into vector
  auc.train.all <- sapply(results, function(x) x$train.AUC)
  # gather all model prediction rasters into a stack and name them
  # if skipRasters is TRUE or no envs, make an empty stack
  if(skipRasters == FALSE & !is.null(envs)) {
    mod.full.pred.all <- raster::stack(sapply(results, function(x) x$mod.full.pred))
    names(mod.full.pred.all) <- tune.names
  }else{
    mod.full.pred.all <- raster::stack()
  }
  # gather all k-fold statistics into a list of data frames,
  # (these are a single set of stats if no partitions were chosen)
  kstats.all <- lapply(results, function(x) x$kstats)
  # define number of grp (the value of "k") as number of
  # rows in one of the model runs
  nk <- nrow(kstats.all[[1]])
  # bind all kstats into a single data frame
  kstats.df <- dplyr::bind_rows(kstats.all)
  
  if(partitions != "none") {
    # if tune settings were specified and there is at least one partition,
    # calculate the kstats tbl
    # if(nrow(tune.tbl) > 0) {
    # define number of settings (plus the tune.args field)
    nset <- ncol(tune.tbl)
    # add in columns for tuning settings
    if(nrow(tune.tbl) > 0) {
      tune.tbl.reps <- as.data.frame(apply(tune.tbl, 2, rep, each = nk))
      kstats.df <- dplyr::bind_cols(tune.tbl.reps, kstats.df)
    }
    # if model settings for tuning were input, summarize by averaging all grp 
    # per model setting combination
    if(nset > 0) kstats.df <- dplyr::group_by_at(kstats.df, 1:nset)
    # set the variables to summarize (excludes the tuning settings and the "fold" column)
    colsExcl <- seq(-1, -(nset+1))
    vars.summarize <- names(kstats.df)[colsExcl]
    # if jackknife cross-validation (leave-one-out), correct variance for
    # non-independent samples (Shcheglovitova & Anderson 2013)
    if(partitions == "jackknife") {
      kstats.avg.df <- kstats.df %>% dplyr::summarize_at(vars.summarize, list(avg = mean, var = ~corrected.var(., nk), 
                                                                            min = min, max = max)) %>%
      dplyr::ungroup()
    }else{
      kstats.avg.df <- kstats.df %>% dplyr::summarize_at(vars.summarize, list(avg = mean, var = var, 
                                                                              min = min, max = max)) %>%
        dplyr::ungroup()
    }
    # change names
    colsExcl <- seq(-1, -nset)
    names(kstats.avg.df)[colsExcl] <- gsub("(.*)_([a-z]{3}$)", "\\1.\\2", names(kstats.avg.df)[colsExcl])
    
    # reorder based on original order of tune names (summarize forces an alphanumeric reorder)
    if(nrow(tune.tbl) > 0) kstats.avg.df <- kstats.avg.df[match(tune.names, kstats.avg.df$tune.args),]
    if(nset > 0) {
      if(partitions == "independent") {
        stats.df <- dplyr::bind_cols(kstats.df[, seq(1, nset)], auc.train = auc.train.all, kstats.df[, seq(-1, -(nset+1))])
      }else{
        stats.df <- dplyr::bind_cols(tune.tbl, auc.train = auc.train.all, kstats.avg.df[, seq(-1, -nset)])  
      }
    }else{
      stats.df <- dplyr::bind_cols(auc.train = auc.train.all, kstats.avg.df) 
    }
  }else{
    # if no partitions were specified, make the stats tbl without cross validation stats
    stats.df <- dplyr::bind_cols(tune.tbl, auc.train = auc.train.all)
  }
  
  # calculate number of non-zero parameters in model
  nparams <- sapply(mod.full.all, enm@nparams)
  # calculate AICc for Maxent models
  if(mod.name %in% c("maxent.jar", "maxnet")) {
    stats.df <- dplyr::bind_cols(stats.df, calc.aicc(nparams, occs, mod.full.pred.all))
  }else{
    warning(paste0("Not able to calculate AICc for ", mod.name, "... returning NAs.\n"))
  }
  stats.df$nparam <- nparams
  # stats.df <- tibble::as_tibble(stats.df)
  kstats.df <- as.data.frame(kstats.df)
  
  res <- list(stats = stats.df, kstats = kstats.df, mods = mod.full.all,
              preds = mod.full.pred.all)
  
  if(is.null(occs.grp)) occs.grp <- 0
  if(is.null(bg.grp)) bg.grp <- 0
  e <- ENMevaluation(algorithm = mod.name, tune.settings = tune.tbl,
                     results = as.data.frame(res$stats), results.grp = res$kstats,
                     predictions = res$preds, models = res$mods, 
                     partition.method = partitions,
                     occ.pts = occs, occ.grp = occs.grp,
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
