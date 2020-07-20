#' @title Tune ecological niche model (ENM) settings and calculate evaluation statistics
#' @description \code{ENMevaluate()} builds ecological niche models iteratively across a range of 
#' user-specified tuning settings. Users can choose to tune with cross validation or an independent 
#' occurrence dataset. \code{ENMevaluate()} returns an \code{ENMevaluation} object with slots containing 
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
#' with the number specified by parameter \code{n.bg}
#' @param tune.args named list of model settings to be tuned
#' @param other.args named list of any additional model arguments not specified for tuning
#' @param categoricals character vector of names of categorical environmental variables
#' @param mod.name character of the name of the chosen model
#' @param user.enm ENMdetails object specified by the user; this model will be
#' used for the analysis, and is an alternative to specifying mod.name
#' @param partitions character of name of partitioning technique (see \code{?partitions})
#' @param user.grp named list with occ.grp = vector of partition group (fold) for each
#' occurrence locality, intended for user-defined partitions, and bg.grp = same vector for 
#' background (or pseudo-absence) localities
#' @param occs.ind matrix or data frame with two columns for longitude and latitude 
#' of occurrence localities, in that order, intended for independent evaluation;
#' when \code{partitions = "independent"}; these occurrences will be used only 
#' for evaluation, and not for model training, and thus no cross validation will be done
#' @param kfolds numeric for number of partition groups (grp), only for random k-fold partitioning
#' @param aggregation.factor numeric vector with length 2 that specifies the factors for aggregating 
#' \code{envs} in order to perform checkerboard partitioning
#' @param n.bg numeric (default: 10000) for number of random background (or pseudo-absence) points to
#'  sample; necessary if \code{bg} is NULL
#' @param overlap boolean (TRUE or FALSE) which if TRUE, calculate niche overlap statistics
#' @param overlapStat character for one or two (vector) niche overlap statistics:
#' choices are "D" and "I" -- see ?calc.niche.overlap for more details
#' @param doClamp boolean (TRUE or FALSE) which if TRUE, clamp model responses; currently only 
#' applicable for maxent.jar/maxnet models
#' @param pred.type character (default: "cloglog") that specifies which prediction type should be used to
#' generate prediction rasters for the ENMevaluation object; currently only applicable for maxent.jar/maxnet models
#' @param cbi.eval character (default: "envs") specifying which should be used to calculate the expected frequency
#' of occurrences for the Continuous Boyce Index: "envs" for the predictions over the entire predictor variable rasters
#' (which necessitates creating a new prediction raster over the full extent for every partition), and "bg" for the 
#' predictions at all localities (training and testing occurrences and backgrounds)
#' @param abs.auc.diff boolean (TRUE or FALSE) which if TRUE, take absolute value of AUCdiff; default is TRUE
#' @param user.test.grps matrix or data frame of user-defined test record coordinates and predictor variable values; this is mainly used
#' internally by ENMnullSims() to force each null model to evaluate with real test data
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

ENMevaluate <- function(occs, envs = NULL, bg = NULL, tune.args = NULL, taxon.name = NULL, other.args = NULL, categoricals = NULL, mod.name = NULL,
                        user.enm = NULL, partitions = NULL, user.grp = NULL, occs.ind = NULL, kfolds = NA, aggregation.factor = c(2, 2), 
                        n.bg = 10000, overlap = FALSE, overlapStat = c("D", "I"), doClamp = TRUE, pred.type = "cloglog", cbi.eval = "envs", 
                        abs.auc.diff = TRUE, user.test.grps = NULL,
                        parallel = FALSE, numCores = NULL, parallelType = "doSNOW", updateProgress = FALSE, 
                        # legacy parameters
                        occ = NULL, env = NULL, bg.coords = NULL, RMvalues = NULL, fc = NULL,
                        algorithm = NULL, method = NULL, bin.output = NULL, rasterPreds = NULL,
                        clamp = NULL, progbar = NULL, occ.grp = NULL, bg.grp = NULL) {
  
  # legacy parameter handling so ENMevaluate doesn't break for older code
  all.legacy <- list(occ, env, bg.coords, RMvalues, fc, algorithm, method, bin.output, rasterPreds, occ.grp,
                     bg.grp, clamp)
  if(sum(sapply(all.legacy, function(x) !is.null(x))) > 0) {
    message("* Running ENMeval v1.0.0 with legacy parameters. These will be phased out in the next version.")
  }
  if(!is.null(occ)) occs <- occ
  if(!is.null(env)) envs <- env
  if (!is.null(bg.coords)) {
    bg <- bg.coords
    names(bg) <- names(occs)
  }
  if(!is.null(method)) partitions <- method
  if(!is.null(rasterPreds)) {
    stop("This parameter was deprecated. If you want to avoid generating model prediction rasters, include predictor variable values in the occurrence and background data frames (SWD format). See Details in ?ENMevaluate for more information.")
  }
  if(!is.null(algorithm)) {
    mod.name <- algorithm
    tune.args <- list(fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
                      rm = seq(0.5, 4, 0.5))
  }
  if(!is.null(RMvalues)) tune.args$rm <- RMvalues
  if(!is.null(fc)) tune.args$fc <- fc
  if(!is.null(occ.grp) & !is.null(bg.grp)) user.grp <- list(occ.grp = occ.grp, bg.grp = bg.grp)
  
  if (!is.null(occ.grp) & !is.null(bg.grp)) {
    user.grp <- list(occ.grp = occ.grp, bg.grp = bg.grp)
  }
  
  if(!is.null(clamp)) doClamp <- clamp
  
  if(is.null(mod.name) & is.null(user.enm)) {
    stop("* Please select a model name (mod.name) or specify a user model (user.enm).")
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
    message(paste0("*** Running initial checks... ***\n"))
  }else{
    message(paste0("*** Running initial checks for ", taxon.name, " ... ***\n"))
  }
  
  ## general parameter checks
  all.partitions <- c("jackknife", "randomkfold", "block", "checkerboard1", 
                      "checkerboard2", "user", "independent", "none")
  
  if(!(partitions %in% all.partitions)) {
    stop("* Please enter an accepted partition method.")
  }
  
  if(partitions == "independent" & is.null(occs.ind)) {
    stop("* If doing independent evaluations, please provide independent testing data (occs.ind).")
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
    message('* As no tuning arguments were specified, turning off niche overlap.')
    overlap <- FALSE
  }
  
  # make sure occs and bg are data frames with identical column names
  if(all(names(occs) != names(bg)) & !is.null(bg)) {
    stop('* Datasets "occs" and "bg" have different column names. Please make them identical and try again.')
  }
  
  if(is.null(envs) & cbi.eval == "envs") {
    warning('* Setting "cbi.eval = envs" requires input of "envs" (predictor variables in raster form).')
    cbi.eval <- "bg"
  }
  
  # if a vector of tuning arguments is numeric, make sure it is sorted (for results table and plotting)
  tune.args.num <- which((sapply(tune.args, class) %in% c("numeric", "integer")) & sapply(tune.args, length) > 1)
  for(i in tune.args.num) {
    tune.args[[i]] <- sort(tune.args[[i]])
  }
  
  # choose a built-in ENMdetails object matching the input model name
  # unless the model is chosen by the user
  if(is.null(user.enm)) {
    enm <- lookup.enm(mod.name)
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
      message(paste0("* Found ", envs.naMismatch, " raster cells that were NA for one or more, but not all, predictor variables. Converting these cells to NA for all predictor variables."))
      envs.names <- names(envs)
      envs <- calc(envs, fun = function(x) if(sum(is.na(x)) > 0) x * NA else x)
      names(envs) <- envs.names
    }
    # if no background points specified, generate random ones
    if(is.null(bg)) {
      message(paste0("* Randomly sampling ", n.bg, " background points ..."))
      bg <- as.data.frame(dismo::randomPoints(envs, n = n.bg))
      names(bg) <- names(occs)
    }
    
    # remove cell duplicates
    occs.cellNo <- raster::extract(envs, occs, cellnumbers = TRUE)
    occs.dups <- duplicated(occs.cellNo[,1])
    if(sum(occs.dups) > 0) message(paste0("* Removed ", sum(occs.dups), " occurrence localities that shared the same grid cell as another."))
    occs <- occs[!occs.dups,]
    if(!is.null(user.grp)) user.grp$occ.grp <- user.grp$occ.grp[!occs.dups]
    
    # bind coordinates to predictor variable values for occs and bg
    occs.vals <- raster::extract(envs, occs)
    bg.vals <- raster::extract(envs, bg)
    occs <- cbind(occs, occs.vals)
    bg <- cbind(bg, bg.vals)
  }else{
    # if no bg included, stop
    if(is.null(bg)) {
      stop("* If inputting variable values without rasters, please make sure to input background coordinates with values as well as occurrences.")
    }
    # for occ and bg coordinates with environmental predictor values (SWD format)
    message("* Variable values were input along with coordinates and not as raster data, so no raster predictions can be generated and AICc cannot be calculated for Maxent models.")
    # make sure both occ and bg have predictor variable values
    if(ncol(occs) < 3 | ncol(bg) < 3) stop("* If inputting variable values without rasters, please make sure these values are included in the occs and bg tables proceeding the coordinates.")
  }
  
  # if NA predictor variable values exist for occs or bg, remove these records and modify user.grp accordingly
  occs.vals.na <- which(rowSums(is.na(occs)) > 0)
  if(length(occs.vals.na) > 0) {
    message(paste0("* Removed ", length(occs.vals.na), " occurrence points with NA predictor variable values."))
    occs <- occs[-occs.vals.na,]
    if(!is.null(user.grp)) user.grp$occ.grp <- user.grp$occ.grp[-occs.vals.na]
  }
  
  bg.vals.na <- which(rowSums(is.na(bg)) > 0)
  if(length(bg.vals.na) > 0) {
    message(paste0("* Removed ", length(bg.vals.na), " background points with NA predictor variable values."))
    bg <- bg[-bg.vals.na,]
    if(!is.null(user.grp)) user.grp$bg.grp <- user.grp$bg.grp[-bg.vals.na]
  }
  
  # make main df with coordinates and predictor variable values
  d <- rbind(occs, bg)
  
  # add presence-background identifier for occs and bg
  d$pb <- c(rep(1, nrow(occs)), rep(0, nrow(bg)))
  
  # if user-defined partitions, assign grp variable before filtering out records with NA predictor variable values
  # for all other partitioning methods, grp assignments occur after filtering
  if(!is.null(user.grp)) {
    d[d$pb == 1, "grp"] <- as.numeric(as.character(user.grp$occ.grp))
    d[d$pb == 0, "grp"] <- as.numeric(as.character(user.grp$bg.grp))
  }
  
  ################################# #
  # ASSIGN CATEGORICAL VARIABLES ####
  ################################# #
  
  # convert fields for categorical data to factor class
  if(!is.null(categoricals)) {
    for(i in 1:length(categoricals)) {
      if(enm@name == "bioclim") {
        message("* As specified model is BIOCLIM, removing categorical variables.")
        d[, categoricals[i]] <- NULL
      }else{
        message(paste0("* Assigning variable ", categoricals[i], " to categorical ..."))
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
  
  # if occs.ind input, coerce partitions to 'independent'
  if(!is.null(occs.ind) & partitions != "independent") {
    partitions <- "independent"
  }
  
  # partition occs based on selected partition method
  grps <- switch(partitions, 
                 jackknife = get.jackknife(d.occs, d.bg),
                 randomkfold = get.randomkfold(d.occs, d.bg, kfolds),
                 block = get.block(d.occs, d.bg),
                 checkerboard1 = get.checkerboard1(d.occs, envs, d.bg, aggregation.factor),
                 checkerboard2 = get.checkerboard2(d.occs, envs, d.bg, aggregation.factor),
                 user = NULL,
                 independent = list(occ.grp = rep(2, nrow(d.occs)), bg.grp = rep(0, nrow(d.bg))),
                 none = list(occ.grp = rep(0, nrow(d.occs)), bg.grp = rep(0, nrow(d.bg))))
  
  
  # choose a user message reporting on partition choice
  parts.message <- switch(partitions,
                      jackknife = "* Model evaluations with k-1 jackknife (leave-one-out) cross validation...",
                      randomkfold = paste0("* Model evaluations with random ", kfolds, "-fold cross validation..."),
                      block = "* Model evaluations with spatial block (4-fold) cross validation...",
                      checkerboard1 = "* Model evaluations with checkerboard (2-fold) cross validation...",
                      checkerboard2 = "* Model evaluations with hierarchical checkerboard (4-fold) cross validation...",
                      user = paste0("* Model evaluations with user-defined ", length(unique(user.grp$occ.grp)), "-fold cross validation..."),
                      independent = "* Model evaluations with independent testing data...",
                      none = "* Skipping model evaluations (only calculating full model statistics)...")
  message(parts.message)
  
  # record partition settings
  parts.settings <- switch(partitions,
                           randomkfold = list(kfolds = kfolds),
                           checkerboard1 = list(aggregation.factor = aggregation.factor),
                           checkerboard2 = list(aggregation.factor = aggregation.factor),
                           user = list(kfolds = length(unique(user.grp$occ.grp))))
  if(is.null(parts.settings)) parts.settings <- list()
  
  # add these values as the 'grp' column
  if(!is.null(grps)) d$grp <- factor(c(grps$occ.grp, grps$bg.grp))
  
  ############################################ #
  # ADD INDEPENDENT TESTING DATA (IF INPUT) ####
  ############################################ #
  
  # add independent tesing data to main df if provided
  if(partitions == "independent") {
    occs.ind.vals <- as.data.frame(raster::extract(envs, occs.ind))
    occs.ind.vals <- cbind(occs.ind, occs.ind.vals)
    occs.ind.vals$pb <- 1
    # the grp here is 1 so that the first cv iteration will evaluate the full dataset on the independent data
    # and the second iteration is not performed
    occs.ind.vals$grp <- 1
    user.test.grps <- occs.ind.vals
    # change the factor levels to accomodate grp 1 (originally it only has grp 2 for occs and grp 0 for bg)
    # d$grp <- factor(d$grp, levels = 0:2)
    # and then add the independent testing data with grp value 1
    # d <- rbind(d, occs.ind.vals)
  }
  
  ##################################### #
  # TURN ON/OFF CBI.TEST CALCULATION ####
  ##################################### #
  
  # for 1) spatial cross validation and 2) jackknife, calculating the continuous Boyce Index on testing data is problematic, as
  # 1) the full study area must be considered, and 2) too few test records are considered, so currently we turn it off
  bg.grp.vals <- unique(d[d$pb==0,"grp"]) == 0
  if(!all(bg.grp.vals) == TRUE | partitions == "jackknife") {
    message("* Turning off test evaluation for Continuous Boyce Index (CBI), as there is no current implementation for jackknife or partitioned background cross-validation (which includes spatial partitioning).")
    cbi.cv <- FALSE
  }else{
    cbi.cv <- TRUE
  }
  
  ################# #
  # MODEL TUNING #### 
  ################# #
  
  # print model-specific message
  if(is.null(taxon.name)) {
    message(paste("\n*** Running ENMeval v1.0.0 with", enm@msgs(tune.args), "***\n"))
  }else{
    message(paste("\n*** Running ENMeval v1.0.0 for", taxon.name, "with", enm@msgs(tune.args), "***\n"))
  }
  
  
  # make table for all tuning parameter combinations
  tune.tbl <- expand.grid(tune.args, stringsAsFactors = FALSE)
  
  # put all settings into list
  other.settings <- list(other.args = other.args, doClamp = doClamp, pred.type = pred.type, 
                         abs.auc.diff = abs.auc.diff, cbi.cv = cbi.cv, cbi.eval = cbi.eval)
  
  if(parallel) {
    results <- tune.parallel(d, envs, enm, partitions, tune.tbl, other.settings, user.test.grps, numCores, parallelType)  
  }else{
    results <- tune.regular(d, envs, enm, partitions, tune.tbl, other.settings, user.test.grps, updateProgress)
  }
  
  ##################### #
  # ASSEMBLE RESULTS #### 
  ##################### #
  
  if(nrow(tune.tbl) == 0) {
    # if not tuned settings, the "tune name" is the model name
    tune.names <- enm@name
  }else{
    # define tuned settings names and bind them to the tune table
    tune.tbl <- dplyr::mutate_all(tune.tbl, as.factor)
    tune.names <- apply(tune.tbl, 1, function(x) paste(x, collapse = "_"))
    tune.tbl$tune.args <- factor(tune.names, levels = tune.names)
  }
  # gather all full models into list and name them
  mod.full.all <- lapply(results, function(x) x$mod.full)
  names(mod.full.all) <- tune.names
  # gather all training AUCs into vector
  train.stats.all <- dplyr::bind_rows(lapply(results, function(x) x$train.stats))
  # gather all model prediction rasters into a stack and name them
  # if envs is null, make an empty stack
  if(!is.null(envs)) {
    mod.full.pred.all <- raster::stack(sapply(results, function(x) x$mod.full.pred))
    names(mod.full.pred.all) <- tune.names
  }else{
    mod.full.pred.all <- raster::stack()
  }
  # gather all k-fold statistics into a list of data frames,
  # (these are a single set of stats if no partitions were chosen)
  cv.stats.all <- dplyr::bind_rows(lapply(results, function(x) x$cv.stats))
  
  # define number of grp (the value of "k") for occurrences
  # k is 1 for partition "independent"
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
    nset <- ncol(tune.tbl)
    
    # if jackknife cross-validation (leave-one-out), correct variance for
    # non-independent samples (Shcheglovitova & Anderson 2013)
    if(partitions == "jackknife") {
      sum.list <- list(avg = mean, sd = ~sqrt(corrected.var(., nk)))
    }else{
      sum.list <- list(avg = mean, sd = sd)
    } 
    
    # if there is one partition, or if using an independent evaluation dataset, do not take summary statistics
    if(nk == 1 | partitions == "independent") sum.list <- list(function(x) {x})
    
    # if tune.tbl exists, make tune.args column a factor to keep order after using dplyr functions
    if(nrow(tune.tbl) > 0) cv.stats.all$tune.args <- factor(cv.stats.all$tune.args, levels = tune.names)
    
    # calculate summary statistics
    cv.stats.sum <- cv.stats.all %>% 
      dplyr::group_by(tune.args) %>%
      dplyr::select(-fold) %>% 
      dplyr::summarize_all(sum.list) %>%
      dplyr::ungroup() 
    
    # change names (replace _ with .)
    names(cv.stats.sum) <- gsub("(.*)_(.*)", "\\1.\\2", names(cv.stats.sum))
    # order columns alphabetically
    cv.stats.sum <- cv.stats.sum[, order(colnames(cv.stats.sum))]
    
    # if tune.tbl exists
    if(nrow(tune.tbl) > 0) {
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
  nparams <- sapply(mod.full.all, enm@nparams)
  
  # calculate AICc
  if((enm@name == "maxnet" | enm@name == "maxent.jar") & !is.null(envs)) {
    pred.type.raw <- switch(enm@name, maxnet = "exponential", maxent.jar = "raw")
    aic.settings <- other.settings
    aic.settings$pred.type <- pred.type.raw
    pred.all.raw <- raster::stack(lapply(mod.full.all, function(x) enm@pred(x, envs, aic.settings)))
    occs.pred.raw <- raster::extract(pred.all.raw, occs[,1:2])
    aic <- enm@aic(occs.pred.raw, nparams, pred.all.raw)
    eval.stats <- dplyr::bind_cols(eval.stats, aic)
  }
  
  # add nparam column
  eval.stats$nparam <- nparams
  
  # assign all groups to 1 if no partitions were selected
  # this avoids putting a NULL in the object slot
  if(nk == 0) d$grp <- 1
  if(is.null(taxon.name)) taxon.name <- ""
  
  e <- ENMevaluation(algorithm = enm@name, tune.settings = tune.tbl,
                     results = eval.stats, results.grp = cv.stats.all,
                     predictions = mod.full.pred.all, models = mod.full.all, 
                     partition.method = partitions, partition.settings = parts.settings,
                     other.settings = other.settings, taxon.name = taxon.name,
                     occs = d[d$pb == 1, 1:(ncol(d)-2)], occ.grp = factor(d[d$pb == 1, "grp"]),
                     bg = d[d$pb == 0, 1:(ncol(d)-2)], bg.grp = factor(d[d$pb == 0, "grp"]),
                     rmm = list())
  # add the rangeModelMetadata object to the ENMevaluation object
  e@rmm <- buildRMM(e, envs)
  
  # if niche overlap selected, calculate and add the resulting matrix to results
  if(overlap == TRUE) {
    nr <- raster::nlayers(e@predictions)
    if(nr == 0) {
      warning("Cannot calculate niche overlap without model prediction rasters.")
    }else if(nr == 1) {
      warning("Only 1 model prediction raster found. Need at least 2 rasters to calculate niche overlap. Increase number of tuning arguments and run again.") 
    }else{
      for(ovStat in overlapStat) {
        message(paste0("Calculating niche overlap for statistic ", ovStat, "..."))
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
