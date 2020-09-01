#' @title Tune ecological niche model (ENM) settings and calculate evaluation statistics
#' @description \code{ENMevaluate()} builds ecological niche models iteratively across a range of 
#' user-specified tuning settings. Users can choose to evaluate models with cross validation or a
#' full-withheld testing dataset. \code{ENMevaluate()} returns an \code{ENMevaluation} object with slots containing 
#' evaluation statistics for each combination of settings and for each cross validation fold therein, as
#' well as raster predictions for each model when raster data is input. The evaluation statistics in the 
#' results table should aid users in identifying model settings that balance fit and predictive ability.
#' 
#' @param occs matrix or data frame with three columns for taxon name, longitude, and latitude 
#' of occurrence localities, in that order; if specifying predictor variable values
#' assigned to presence/background localities (species with data "SWD" form), this table will also have 
#' one column for each predictor variable
#' @param envs Raster* object of environmental variables (must be in same geographic projection as occurrence data)
#' @param bg matrix or data frame with two columns for longitude and latitude of 
#' background (or pseudo-absence) localities, in that order; if specifying predictor variable values
#' assigned to presence/background localities (species with data "SWD" form), this table will also have 
#' one column for each predictor variable; if NULL, points will be randomly sampled across \code{envs} 
#' with the number specified by argument \code{n.bg}
#' @param tune.args named list of model settings to be tuned
#' @param other.args named list of any additional model arguments not specified for tuning
#' @param categoricals character vector of names of categorical environmental variables
#' @param algorithm character of the name of the chosen model
#' @param user.enm ENMdetails object specified by the user; this model will be
#' used for the analysis, and is an alternative to specifying algorithm
#' @param partitions character of name of partitioning technique (see \code{?partitions})
#' @param user.grp named list with occs.grp = vector of partition group (fold) for each
#' occurrence locality, intended for user-defined partitions, and bg.grp = same vector for 
#' background (or pseudo-absence) localities
#' @param occs.testing matrix or data frame with two columns for longitude and latitude 
#' of occurrence localities, in that order, intended for a testing dataset;
#' when \code{partitions = "testing"}; these occurrences will be used only 
#' for evaluation, and not for model training, and thus no cross validation will be done
#' @param kfolds numeric for number of partition groups (grp), only for random k-fold partitioning
#' @param aggregation.factor numeric vector with length 2 that specifies the factors for aggregating 
#' \code{envs} in order to perform checkerboard partitioning
#' @param n.bg numeric (default: 10000) for number of random background (or pseudo-absence) points to
#'  sample; necessary if \code{bg} is NULL
#' @param overlap boolean (TRUE or FALSE) which if TRUE, calculate niche overlap statistics
#' @param overlapStat character for one or two (vector) niche overlap statistics:
#' choices are "D" and "I" -- see ?calc.niche.overlap for more details
#' @param clamp boolean (TRUE or FALSE) which if TRUE, clamp model responses; currently only 
#' applicable for maxent.jar/maxnet models
#' @param pred.type character (default: "cloglog") that specifies which prediction type should be used to
#' generate prediction rasters for the ENMevaluation object; currently only applicable for maxent.jar/maxnet models
#' @param abs.auc.diff boolean (TRUE or FALSE) which if TRUE, take absolute value of AUCdiff; default is TRUE
#' @param validation.bg Character: either "full" to calculate AUC and CBI with respect to the full background, or
#' "partition" to calculate them with respect to the validation partition background 
#' @param user.val.grps matrix or data frame of user-defined test record coordinates and predictor variable values; this is mainly used
#' internally by ENMnullSims() to force each null model to evaluate with real test data
#' @param user.eval function specifying custom validation evaluation (see vignette for example)
#' @param rmm RMM object to be written to by ENMevaluate; if not specified, a new RMM object is output
#' @param parallel boolean (TRUE or FALSE) which if TRUE, run with parallel processing
#' @param numCores numeric for number of cores to use for parallel processing
#' @param parallelType character (default: "doSNOW") specifying either "doParallel" or "doSNOW"
#' @param updateProgress boolean (TRUE or FALSE) which if TRUE, use shiny progress bar; only for use in shiny apps
#'
#' @return 
#'
#' @examples
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#' @export 
#' 

ENMevaluate <- function(occs, envs = NULL, bg = NULL, tune.args = NULL, taxon.name = NULL, other.args = NULL, categoricals = NULL, algorithm = NULL,
                        user.enm = NULL, partitions = NULL, user.grp = NULL, occs.testing = NULL, 
                        kfolds = NA, aggregation.factor = c(2, 2), orientation = "lat_lon",
                        n.bg = 10000, overlap = FALSE, overlapStat = c("D", "I"), clamp = TRUE, pred.type = "cloglog", 
                        abs.auc.diff = TRUE, validation.bg = "full", user.val.grps = NULL, user.eval = NULL, rmm = NULL,
                        parallel = FALSE, numCores = NULL, parallelType = "doSNOW", updateProgress = FALSE, quiet = FALSE, 
                        # legacy parameters
                        occ = NULL, env = NULL, bg.coords = NULL, RMvalues = NULL, fc = NULL, occ.grp = NULL, bg.grp = NULL,
                        method = NULL, bin.output = NULL, rasterPreds = NULL, progbar = NULL) {
  
  # legacy argument handling so ENMevaluate doesn't break for older code
  all.legacy <- list(occ, env, bg.coords, RMvalues, fc, occ.grp, bg.grp, method, bin.output, rasterPreds)
  if(sum(sapply(all.legacy, function(x) !is.null(x))) > 0) {
    if(quiet != TRUE) message("* Running ENMeval v2.0.0 with legacy parameters. These will be phased out in the next version.")
  }
  if(!is.null(occ)) occs <- occ
  if(!is.null(env)) envs <- env
  if(!is.null(bg.coords)) bg <- bg.coords
  if(!is.null(method)) partitions <- method
  if(!is.null(rasterPreds)) {
    stop("Warning: This argument was deprecated. If you want to avoid generating model prediction rasters, include predictor variable values in the occurrence and background data frames (SWD format). See Details in ?ENMevaluate for more information.")
  }
  if(!is.null(RMvalues)) tune.args$rm <- RMvalues
  if(!is.null(fc)) tune.args$fc <- fc
  if(!is.null(occ.grp) & !is.null(bg.grp)) {
    user.grp <- list(occs.grp = occ.grp, bg.grp = bg.grp)
    if(quiet != TRUE) message('Warning: These arguments were deprecated and replaced with the argument "user.grp".')
  }
  if((!is.null(occ.grp) & is.null(bg.grp)) | (is.null(occ.grp) & !is.null(bg.grp))) {
    stop('For user partitions, please input both arguments "occ.grp" and "bg.grp". Warning: These are legacy parameters that were replaced with the argument "user.grp".')
  }
  
  if(is.null(algorithm) & is.null(user.enm)) {
    stop('* Please select a model name (argument "algorithm") or specify a user model (argument "user.enm").')
  }
  
  # record start time
  start.time <- proc.time()
  
  ######################## #
  # INITIAL DATA CHECKS ####
  ######################## #
  
  # coerce occs and bg to df
  occs <- as.data.frame(occs)
  if(!is.null(bg)) bg <- as.data.frame(bg)
  # extract species name and coordinates
  
  if(is.null(taxon.name)) {
    if(quiet != TRUE) message(paste0("*** Running initial checks... ***\n"))
  }else{
    if(quiet != TRUE) message(paste0("*** Running initial checks for ", taxon.name, " ... ***\n"))
  }
  
  ## general argument checks
  all.partitions <- c("jackknife", "randomkfold", "block", "checkerboard1", 
                      "checkerboard2", "user", "testing", "none")
  
  if(!(partitions %in% all.partitions)) {
    stop("* Please enter an accepted partition method.")
  }
  
  if(partitions == "testing" & is.null(occs.testing)) {
    stop("* If doing testing evaluations, please provide testing data (occs.testing).")
  }
  
  if((partitions == "checkerboard1" | partitions == "checkerboard2") & is.null(envs)) {
    stop('* For checkerboard partitioning, predictor variable rasters "envs" are required.')
  }
  
  if(partitions == "randomkfold") {
    if(is.null(kfolds)) {
      stop('* For random k-fold partitioning, a numeric, non-zero value of "kfolds" is required.')  
    }else{
      if(kfolds == 0) {
        stop('* For random k-fold partitioning, a numeric, non-zero value of "kfolds" is required.')  
      }
    }
  }
  
  if(is.null(tune.args) & overlap == TRUE) {
    if(quiet != TRUE) message('* As no tuning arguments were specified, turning off niche overlap.')
    overlap <- FALSE
  }
  
  # make sure occs and bg are data frames with identical column names
  if(all(names(occs) != names(bg)) & !is.null(bg)) {
    stop('* Datasets "occs" and "bg" have different column names. Please make them identical and try again.')
  }
  
  # if a vector of tuning arguments is numeric, make sure it is sorted (for results table and plotting)
  tune.args.num <- which((sapply(tune.args, class) %in% c("numeric", "integer")) & sapply(tune.args, length) > 1)
  for(i in tune.args.num) {
    tune.args[[i]] <- sort(tune.args[[i]])
  }
  
  # choose a built-in ENMdetails object matching the input model name
  # unless the model is chosen by the user
  if(is.null(user.enm)) {
    enm <- lookup.enm(algorithm)
  }else{
    enm <- user.enm
  }
  
  ########################################################### #
  # ASSEMBLE COORDINATES AND ENVIRONMENTAL VARIABLE VALUES ####
  ########################################################### #
  
  # if environmental rasters are input as predictor variables
  if(!is.null(envs)) {
    # make sure envs is a RasterStack -- if RasterLayer, maxent.jar crashes
    envs <- raster::stack(envs)
    envs.z <- raster::values(envs)
    envs.naMismatch <- sum(apply(envs.z, 1, function(x) !all(is.na(x)) & !all(!is.na(x))))
    if(envs.naMismatch > 0) {
      if(quiet != TRUE) message(paste0("* Found ", envs.naMismatch, " raster cells that were NA for one or more, but not all, predictor variables. Converting these cells to NA for all predictor variables."))
      envs.names <- names(envs)
      envs <- calc(envs, fun = function(x) if(sum(is.na(x)) > 0) x * NA else x)
      names(envs) <- envs.names
    }
    # if no background points specified, generate random ones
    if(is.null(bg)) {
      if(quiet != TRUE) message(paste0("* Randomly sampling ", n.bg, " background points ..."))
      bg <- as.data.frame(dismo::randomPoints(envs, n = n.bg))
      names(bg) <- names(occs)
    }
    
    # remove cell duplicates
    occs.cellNo <- raster::extract(envs, occs, cellnumbers = TRUE)
    occs.dups <- duplicated(occs.cellNo[,1])
    if(sum(occs.dups) > 0) if(quiet != TRUE) message(paste0("* Removed ", sum(occs.dups), " occurrence localities that shared the same grid cell as another."))
    occs <- occs[!occs.dups,]
    if(!is.null(user.grp)) user.grp$occs.grp <- user.grp$occs.grp[!occs.dups]
    
    # bind coordinates to predictor variable values for occs and bg
    occs.z <- raster::extract(envs, occs)
    bg.z <- raster::extract(envs, bg)
    occs <- cbind(occs, occs.z)
    bg <- cbind(bg, bg.z)
  }else{
    # if no bg included, stop
    if(is.null(bg)) {
      stop("* If inputting variable values without rasters, please make sure to input background coordinates with values as well as occurrences.")
    }
    # for occ and bg coordinates with environmental predictor values (SWD format)
    if(quiet != TRUE) message("* Variable values were input along with coordinates and not as raster data, so no raster predictions can be generated and AICc cannot be calculated for Maxent models.")
    # make sure both occ and bg have predictor variable values
    if(ncol(occs) < 3 | ncol(bg) < 3) stop("* If inputting variable values without rasters, please make sure these values are included in the occs and bg tables proceeding the coordinates.")
  }
  
  # if NA predictor variable values exist for occs or bg, remove these records and modify user.grp accordingly
  occs.z.na <- which(rowSums(is.na(occs)) > 0)
  if(length(occs.z.na) > 0) {
    if(quiet != TRUE) message(paste0("* Removed ", length(occs.z.na), " occurrence points with NA predictor variable values."))
    occs <- occs[-occs.z.na,]
    if(!is.null(user.grp)) user.grp$occs.grp <- user.grp$occs.grp[-occs.z.na]
  }
  
  bg.z.na <- which(rowSums(is.na(bg)) > 0)
  if(length(bg.z.na) > 0) {
    if(quiet != TRUE) message(paste0("* Removed ", length(bg.z.na), " background points with NA predictor variable values."))
    bg <- bg[-bg.z.na,]
    if(!is.null(user.grp)) user.grp$bg.grp <- user.grp$bg.grp[-bg.z.na]
  }
  
  # make main df with coordinates and predictor variable values
  d <- rbind(occs, bg)
  
  # add presence-background identifier for occs and bg
  d$pb <- c(rep(1, nrow(occs)), rep(0, nrow(bg)))
  
  # if user-defined partitions, assign grp variable before filtering out records with NA predictor variable values
  # for all other partitioning methods, grp assignments occur after filtering
  if(!is.null(user.grp)) {
    d[d$pb == 1, "grp"] <- as.numeric(as.character(user.grp$occs.grp))
    d[d$pb == 0, "grp"] <- as.numeric(as.character(user.grp$bg.grp))
  }
  
  ################################# #
  # ASSIGN CATEGORICAL VARIABLES ####
  ################################# #
  
  # convert fields for categorical data to factor class
  if(!is.null(categoricals)) {
    for(i in 1:length(categoricals)) {
      if(enm@name == "bioclim") {
        if(quiet != TRUE) message("* As specified model is BIOCLIM, removing categorical variables.")
        d[, categoricals[i]] <- NULL
        envs <- raster::dropLayer(envs, categoricals)
      }else{
        if(quiet != TRUE) message(paste0("* Assigning variable ", categoricals[i], " to categorical ..."))
        d[, categoricals[i]] <- as.factor(d[, categoricals[i]])  
      }
    }
  }
  
  ###################### #
  # ASSIGN PARTITIONS ####
  ###################### #
  
  # unpack occs and bg records for partitioning
  d.occs <- d %>% dplyr::filter(pb == 1) %>% dplyr::select(1:2)
  d.bg <- d %>% dplyr::filter(pb == 0) %>% dplyr::select(1:2)
  
  # if occs.testing input, coerce partitions to 'testing'
  if(!is.null(occs.testing) & partitions != "testing") {
    partitions <- "testing"
  }
  
  # partition occs based on selected partition method
  grps <- switch(partitions, 
                 jackknife = get.jackknife(d.occs, d.bg),
                 randomkfold = get.randomkfold(d.occs, d.bg, kfolds),
                 block = get.block(d.occs, d.bg, orientation),
                 checkerboard1 = get.checkerboard1(d.occs, envs, d.bg, aggregation.factor),
                 checkerboard2 = get.checkerboard2(d.occs, envs, d.bg, aggregation.factor),
                 user = NULL,
                 testing = list(occs.grp = rep(2, nrow(d.occs)), bg.grp = rep(0, nrow(d.bg))),
                 none = list(occs.grp = rep(0, nrow(d.occs)), bg.grp = rep(0, nrow(d.bg))))
  
  
  # choose a user message reporting on partition choice
  parts.message <- switch(partitions,
                      jackknife = "* Model evaluations with k-1 jackknife (leave-one-out) cross validation...",
                      randomkfold = paste0("* Model evaluations with random ", kfolds, "-fold cross validation..."),
                      block =  paste0("* Model evaluations with spatial block (4-fold) cross validation and ", orientation, " orientation..."),
                      checkerboard1 = "* Model evaluations with checkerboard (2-fold) cross validation...",
                      checkerboard2 = "* Model evaluations with hierarchical checkerboard (4-fold) cross validation...",
                      user = paste0("* Model evaluations with user-defined ", length(unique(user.grp$occs.grp)), "-fold cross validation..."),
                      testing = "* Model evaluations with testing data...",
                      none = "* Skipping model evaluations (only calculating full model statistics)...")
  if(quiet != TRUE) message(parts.message)
  
  # record partition settings
  parts.settings <- switch(partitions,
                           randomkfold = list(kfolds = kfolds),
                           checkerboard1 = list(aggregation.factor = aggregation.factor),
                           checkerboard2 = list(aggregation.factor = aggregation.factor),
                           user = list(kfolds = length(unique(user.grp$occs.grp))))
  if(is.null(parts.settings)) parts.settings <- list()
  
  # if jackknife partitioning, do not calculate CBI because there are too few validation occurrences
  # per partition (n=1) to have a meaningful result
  if(partitions == "jackknife") cbi.cv <- FALSE else cbi.cv <- TRUE
  
  # add these values as the 'grp' column
  if(!is.null(grps)) d$grp <- factor(c(grps$occs.grp, grps$bg.grp))
  
  ############################################ #
  # ADD TESTING DATA (IF INPUT) ####
  ############################################ #
  
  # add testing data to main df if provided
  if(partitions == "testing") {
    occs.testing.z <- as.data.frame(raster::extract(envs, occs.testing))
    occs.testing.z <- cbind(occs.testing, occs.testing.z)
    occs.testing.z$pb <- 1
    # the grp here is 1 so that the first cv iteration will evaluate the full dataset on the testing data
    # and the second iteration is not performed
    occs.testing.z$grp <- 1
    user.val.grps <- occs.testing.z
    # change the factor levels to accomodate grp 1 (originally it only has grp 2 for occs and grp 0 for bg)
    # d$grp <- factor(d$grp, levels = 0:2)
    # and then add the testing data with grp value 1
    # d <- rbind(d, occs.testing.z)
  }
  
  ################# #
  # MODEL TUNING #### 
  ################# #
  
  # print model-specific message
  if(is.null(taxon.name)) {
    if(quiet != TRUE) message(paste("\n*** Running ENMeval v2.0.0 with", enm@msgs(tune.args), "***\n"))
  }else{
    if(quiet != TRUE) message(paste("\n*** Running ENMeval v2.0.0 for", taxon.name, "with", enm@msgs(tune.args), "***\n"))
  }
  
  # make table for all tuning parameter combinations
  tune.tbl <- expand.grid(tune.args, stringsAsFactors = FALSE)
  # make tune.tbl NULL, not an empty table, if no settings are specified
  # this makes it easier to use tune.i as a parameter in function calls
  # when tune.args does not exist
  if(nrow(tune.tbl) == 0) tune.tbl <- NULL
  
  # put all settings into list
  other.settings <- list(other.args = other.args, clamp = clamp, pred.type = pred.type, 
                         abs.auc.diff = abs.auc.diff, validation.bg = validation.bg, cbi.cv = cbi.cv)
  
  if(parallel) {
    results <- tune.parallel(d, envs, enm, partitions, tune.tbl, other.settings, user.val.grps, numCores, parallelType, user.eval, quiet)  
  }else{
    results <- tune.regular(d, envs, enm, partitions, tune.tbl, other.settings, user.val.grps, updateProgress, user.eval, quiet)
  }
  
  ##################### #
  # ASSEMBLE RESULTS #### 
  ##################### #
  
  # flatten all training statistics data frames from results list into a single data frame
  train.stats.all <- dplyr::bind_rows(lapply(results, function(x) x$train.stats))
  # flatten all validation statistics data frames from results list into a single data frame
  # (these are no validation stats if no partitions were chosen)
  val.stats.all <- dplyr::bind_rows(lapply(results, function(x) x$cv.stats))
  
  if(is.null(tune.tbl)) {
    # if not tuned settings, the "tune name" is the model name
    tune.names <- enm@name
  }else{
    # define tuned settings names and bind them to the tune table
    tune.tbl <- dplyr::mutate_all(tune.tbl, as.factor)
    tune.names <- train.stats.all$tune.args
    tune.tbl$tune.args <- factor(tune.names, levels = tune.names)
  }
  # gather all full models into list and name them
  mod.full.all <- lapply(results, function(x) x$mod.full)
  names(mod.full.all) <- tune.names
  
  # gather all model prediction rasters into a stack and name them
  # if envs is null, make an empty stack
  if(!is.null(envs)) {
    mod.full.pred.all <- raster::stack(sapply(results, function(x) x$mod.full.pred))
    names(mod.full.pred.all) <- tune.names
  }else{
    mod.full.pred.all <- raster::stack()
  }
  
  # define number of grp (the value of "k") for occurrences
  # k is 1 for partition "testing"
  # k is 0 for partitions "none" and "user"
  occGrps <- unique(d[d$pb == 1, "grp"])
  if(length(occGrps) == 1 & 0 %in% occGrps) {
    nk <- 0
  }else{
    nk <- length(occGrps)
  }
  
  # if partitions were specified
  if(nk > 0) {
    # define number of settings (plus the tune.args field)
    nset <- ifelse(!is.null(tune.tbl), ncol(tune.tbl), 0)
    
    # if jackknife cross-validation (leave-one-out), correct variance for
    # non-independent samples (Shcheglovitova & Anderson 2013)
    if(partitions == "jackknife") {
      sum.list <- list(avg = mean, sd = ~sqrt(corrected.var(., nk)))
    }else{
      sum.list <- list(avg = mean, sd = sd)
    } 
    
    # if there is one partition, or if using an testing dataset, do not take summary statistics
    if(nk == 1 | partitions == "testing") sum.list <- list(function(x) {x})
    
    # if tune.tbl exists, make tune.args column a factor to keep order after using dplyr functions
    if(!is.null(tune.tbl)) val.stats.all$tune.args <- factor(val.stats.all$tune.args, levels = tune.names)
    
    # calculate summary statistics
    cv.stats.sum <- val.stats.all %>% 
      dplyr::group_by(tune.args) %>%
      dplyr::select(-fold) %>% 
      dplyr::summarize_all(sum.list) %>%
      dplyr::ungroup() 
    
    # change names (replace _ with .)
    names(cv.stats.sum) <- gsub("(.*)_(.*)", "\\1.\\2", names(cv.stats.sum))
    # order columns alphabetically
    cv.stats.sum <- cv.stats.sum[, order(colnames(cv.stats.sum))]
    
    # if tune.tbl exists
    if(!is.null(tune.tbl)) {
      # make tune.args column in training stats factor too for smooth join
      train.stats.all$tune.args <- factor(train.stats.all$tune.args, levels = tune.names)
      # eval.stats is the join of tune.tbl, training stats, and cv stats
      eval.stats <- tune.tbl %>% dplyr::left_join(train.stats.all, by = "tune.args") %>%
        dplyr::left_join(cv.stats.sum, by = "tune.args")
    }else{
      # if tune.tbl does not exist, eval.stats is the binding of training stats to cv stats
      train.stats.all$tune.args <- NULL
      cv.stats.sum$tune.args <- NULL
      eval.stats <- dplyr::bind_cols(train.stats.all, cv.stats.sum)
    }
  }else{
    # make tune.args column in training stats factor too for smooth join
    train.stats.all$tune.args <- factor(train.stats.all$tune.args, levels = tune.names)
    # if no partitions assigned, eval.stats is the join of tune.tbl to training stats
    eval.stats <- dplyr::left_join(tune.tbl, train.stats.all, by = "tune.args")
  }
  
  # calculate number of non-zero parameters in model
  ncoefs <- sapply(mod.full.all, enm@ncoefs)
  
  # calculate AICc
  if((enm@name == "maxnet" | enm@name == "maxent.jar")) {
    pred.type.raw <- switch(enm@name, maxnet = "exponential", maxent.jar = "raw")
    aic.settings <- other.settings
    aic.settings$pred.type <- pred.type.raw
    if(!is.null(envs)) {
      pred.all.raw <- raster::stack(lapply(mod.full.all, enm@predict, envs, aic.settings))
    }else{
      pred.all.raw <- NULL
    }
    occs.pred.raw <- dplyr::bind_rows(lapply(mod.full.all, enm@predict, occs[,-c(1,2)], aic.settings))
    aic <- aic.maxent(occs.pred.raw, ncoefs, pred.all.raw)
    eval.stats <- dplyr::bind_cols(eval.stats, aic)
  }
  
  # add ncoef column
  eval.stats$ncoef <- ncoefs
  
  # assign all groups to 1 if no partitions were selected
  # this avoids putting a NULL in the object slot
  if(nk == 0) d$grp <- 1
  if(is.null(taxon.name)) taxon.name <- ""
  if(is.null(tune.tbl)) tune.tbl <- data.frame()
  
  # get variable importance for all models
  varimp.all <- lapply(mod.full.all, enm@varimp)
  
  # assemble the ENMevaluation object
  e <- ENMevaluation(algorithm = enm@name, tune.settings = tune.tbl,
                     results = eval.stats, results.partitions = val.stats.all,
                     predictions = mod.full.pred.all, models = mod.full.all, 
                     variable.importance = varimp.all,
                     partition.method = partitions, partition.settings = parts.settings,
                     other.settings = other.settings, taxon.name = taxon.name,
                     occs = d[d$pb == 1, 1:(ncol(d)-2)], occs.grp = factor(d[d$pb == 1, "grp"]),
                     bg = d[d$pb == 0, 1:(ncol(d)-2)], bg.grp = factor(d[d$pb == 0, "grp"]),
                     rmm = list())
  
  # add the rangeModelMetadata object to the ENMevaluation object
  # write to existing RMM if input by user
  e@rmm <- buildRMM(e, envs, rmm)
  
  # if niche overlap selected, calculate and add the resulting matrix to results
  if(overlap == TRUE) {
    nr <- raster::nlayers(e@predictions)
    if(nr == 0) {
      warning("Cannot calculate niche overlap without model prediction rasters.")
    }else if(nr == 1) {
      warning("Only 1 model prediction raster found. Need at least 2 rasters to calculate niche overlap. Increase number of tuning arguments and run again.") 
    }else{
      for(ovStat in overlapStat) {
        if(quiet != TRUE) message(paste0("Calculating niche overlap for statistic ", ovStat, "..."))
        # turn negative values to 0 for niche overlap calculations
        predictions.noNegs <- calc(e@predictions, function(x) {x[x<0] <- 0; x})
        overlap.mat <- calc.niche.overlap(predictions.noNegs, ovStat, quiet)
        e@overlap[[ovStat]] <- overlap.mat
      }
    }
  }
  
  # calculate time expended and print message
  timed <- proc.time() - start.time
  t.min <- floor(timed[3] / 60)
  t.sec <- timed[3] - (t.min * 60)
  if(quiet != TRUE) message(paste("ENMevaluate completed in", t.min, "minutes", round(t.sec, 1), "seconds."))
  
  return(e)
}
