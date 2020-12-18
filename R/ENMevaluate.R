#' @title Tune ecological niche model (ENM) settings and calculate evaluation statistics
#' @description \code{ENMevaluate()} builds ecological niche models iteratively across a range of 
#' user-specified tuning settings. Users can choose to evaluate models with cross validation or a
#' full-withheld testing dataset. \code{ENMevaluate()} returns an \code{ENMevaluation} object with slots containing 
#' evaluation statistics for each combination of settings and for each cross validation fold therein, as
#' well as raster predictions for each model when raster data is input. The evaluation statistics in the 
#' results table should aid users in identifying model settings that balance fit and predictive ability.
#' @details There are a few methodological details in the implementation of ENMeval 2.0 that are important to mention.
#' They are also discussed briefly in ?other.settings and ?ENMnulls.
#' 
#' 1. By default, validation AUC is calculated with respect to the full background (training + validation).
#' This approach follows Radosavljevic & Anderson (2014).This setting can be changed by assigning 
#' other.settings$validation.bg to "partition", which will calculate AUC with respect 
#' to the validation background only. The default value for other.settings$validation.bg is "full".
#' 
#' 2. The continuous Boyce index is not calculated with respect to the RasterStack delineating the study extent,
#' but instead to the background records. This decision was made to simplify the code and improve running time. 
#' If the background records are a good representation of the study extent, there should not be much difference
#' between this and the raster approach.
#' 
#' 3. Null occurrences for null ENMs are sample randomly from the background records, not the RasterStack
#' delineating the study extent. This decision was made for similar reasons to that for CBI, and as above,
#' there should be little difference as long as the background records represent the study extent well.
#' 
#' @param occs matrix / data frame: occurrence records with two columns for longitude and latitude 
#' of occurrence localities, in that order -- if specifying predictor variable values
#' assigned to presence/background localities (species with data "SWD" form), this table will also have 
#' one column for each predictor variable
#' @param envs RasterStack: environmental predictor variables (must be in same geographic projection as occurrence data)
#' @param bg matrix / data frame: background records with two columns for longitude and latitude of 
#' background (or pseudo-absence) localities, in that order; if specifying predictor variable values
#' assigned to presence/background localities (species with data "SWD" form), this table will also have 
#' one column for each predictor variable; if NULL, points will be randomly sampled across \code{envs} 
#' with the number specified by argument \code{n.bg}
#' @param tune.args named list: model settings to be tuned (i.e., list(fc = c("L","Q"), rm = 1:3))
#' @param partitions character: name of partitioning technique (see \code{?partitions})
#' @param algorithm character: name of the algorithm used to build models -- one of "maxnet" or
#' "maxent.jar", else the name from a custom ENMdetails implementation
#' @param partition.settings named list: settings specific to certain partitions -- see ?partition.settings
#' @param other.settings named list: other settings for analysis -- see ?other.settings
#' @param categoricals character vector: name or names of categorical environmental variables -- if not specified,
#' all predictor variables will be treated as continuous unless they are factors (if categorical variables
#' are already factors, specifying names of such variables in this argument is not needed)
#' @param doClamp boolean: if TRUE (default), model prediction extrapolations will be restricted to the upper and lower
#' bounds of the predictor variables -- this avoids extreme predictions for non-analog environments, but
#' if extrapolation is a study aim, this should be set to FALSE
#' @param clamp.directions named list: specifies the direction ("left" for minimum, "right" for maximum) 
#' of clamping for predictor variables -- (e.g.) list(left = c("bio1","bio5"), right = c("bio10","bio15"))
#' @param user.enm ENMdetails object: an alternative to specifying an algorithm, users can insert a custom
#' ENMdetails object to build models with
#' @param user.grp named list: specifies user-defined partition groups, where occs.grp = vector of partition group (fold) for each
#' occurrence locality, intended for user-defined partitions, and bg.grp = same vector for background (or pseudo-absence) localities
#' @param occs.testing matrix / data frame: a full withheld testing dataset with two columns for longitude and latitude 
#' of occurrence localities, in that order -- when \code{partitions = "testing"}; these occurrences will be used only 
#' for evaluation but not for model training, and thus no cross validation will be performed
#' @param taxon.name character: name of the focal species or taxon -- used primarily for annotating
#' the ENMevaluation object and output metadata (rmm), but not necessary
#' @param n.bg numeric: if background records not already provided, this specifies the 
#' number of background (or pseudo-absence) points to randomly sample over envs raster (default: 10000)
#' @param overlap boolean: if TRUE, calculate niche overlap statistics
#' @param overlapStat character: niche overlap statistics to be calculated -- 
#' "D" (Schoener's D) and or "I" (Hellinger's I) -- see ?calc.niche.overlap for more details
#' @param user.val.grps matrix / data frame: user-defined validation record coordinates and predictor variable values -- 
#' this is used internally by ENMnulls() to force each null model to evaluate with empirical validation data
#' @param user.eval function: specify custom validation evaluation (see vignette for example)
#' @param rmm rangeModelMetadata object: if specified, ENMevaluate() will write metadata details for the analysis into
#' this object, but if not, a new rangeModelMetadata object will be generated and written to
#' @param parallel boolean: if TRUE, run with parallel processing
#' @param numCores numeric: number of cores to use for parallel processing; if NULL, all available cores will be used
#' @param parallelType character:: either "doParallel" or "doSNOW" (default: "doSNOW") 
#' @param updateProgress boolean: if TRUE, use shiny progress bar; only for use in shiny apps
#' @param quiet boolean: if TRUE, silence all function messages (but not errors)
#' @param legacy.arguments these are included to avoid unnecessary errors for older scripts, but in a later version
#' these arguments will be permanently deprecated
#' 
#' @return An ENMevaluation object. See ?ENMevaluation for details.
#'
#' @examples
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#' @export 
#' 

ENMevaluate <- function(occs, envs = NULL, bg = NULL, tune.args = NULL, partitions = NULL, algorithm = NULL, 
                        partition.settings = list(orientation = "lat_lon", aggregation.factor = 2, kfolds = 5), 
                        other.settings = list(abs.auc.diff = TRUE, pred.type = "cloglog", validation.bg = "full"), 
                        categoricals = NULL, doClamp = TRUE, clamp.directions = NULL,
                        user.enm = NULL, user.grp = NULL, occs.testing = NULL, taxon.name = NULL, 
                        n.bg = 10000, overlap = FALSE, overlapStat = c("D", "I"), 
                        user.val.grps = NULL, user.eval = NULL, rmm = NULL,
                        parallel = FALSE, numCores = NULL, parallelType = "doSNOW", updateProgress = FALSE, quiet = FALSE, 
                        # legacy arguments
                        occ = NULL, env = NULL, bg.coords = NULL, RMvalues = NULL, fc = NULL, occ.grp = NULL, bg.grp = NULL,
                        method = NULL, bin.output = NULL, rasterPreds = NULL, clamp = NULL, progbar = NULL) {
  
  # legacy argument handling so ENMevaluate doesn't break for older code
  all.legacy <- list(occ, env, bg.coords, RMvalues, fc, occ.grp, bg.grp, method, bin.output, rasterPreds)
  if(sum(sapply(all.legacy, function(x) !is.null(x))) > 0) {
    if(quiet != TRUE) message("* Running ENMeval v2.0.0 with legacy arguments. These will be phased out in the next version.")
  }
  if(!is.null(occ)) occs <- occ
  if(!is.null(env)) envs <- env
  if(!is.null(bg.coords)) bg <- bg.coords
  if(!is.null(method)) partitions <- method
  if(!is.null(clamp)) doClamp <- clamp
  if(!is.null(rasterPreds)) {
    stop("This argument was deprecated. If you want to avoid generating model prediction rasters, include predictor variable values in the occurrence and background data frames (SWD format). See Details in ?ENMevaluate for more information.")
  }
  if(!is.null(RMvalues)) tune.args$rm <- RMvalues
  if(!is.null(fc)) tune.args$fc <- fc
  if(!is.null(occ.grp) & !is.null(bg.grp)) {
    user.grp <- list(occs.grp = occ.grp, bg.grp = bg.grp)
    if(quiet != TRUE) message('Warning: These arguments were deprecated and replaced with the argument "user.grp".')
  }
  if((!is.null(occ.grp) & is.null(bg.grp)) | (is.null(occ.grp) & !is.null(bg.grp))) {
    stop('For user partitions, please input both arguments "occ.grp" and "bg.grp". Warning: These are legacy arguments that were replaced with the argument "user.grp".')
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
    if(is.null(partition.settings$kfolds)) {
      stop('* For random k-fold partitioning, a numeric, non-zero value of "kfolds" is required.')  
    }else{
      if(partition.settings$kfolds == 0) {
        stop('* For random k-fold partitioning, a numeric, non-zero value of "kfolds" is required.')  
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
  
  if(is.null(tune.args) & overlap == TRUE) {
    if(quiet != TRUE) message('* As no tuning arguments were specified, turning off niche overlap.')
    overlap <- FALSE
  }
  
  # make sure occs and bg are data frames with identical column names
  if(all(names(occs) != names(bg)) & !is.null(bg)) {
    stop('* Datasets "occs" and "bg" have different column names. Please make them identical and try again.')
  }
  
  if(doClamp == FALSE & !is.null(clamp.directions)) {
    stop("If specifying clamp directions, please make doClamp = TRUE.")
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
      envs <- raster::stack(calc(envs, fun = function(x) if(sum(is.na(x)) > 0) x * NA else x))
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
  
  ############################################ #
  # ADD TESTING DATA (IF INPUT) ####
  ############################################ #
  
  # add testing data to main df if provided
  if(partitions == "testing") {
    if(!is.null(envs)) {
      occs.testing.z <- as.data.frame(raster::extract(envs, occs.testing))
      occs.testing.z <- cbind(occs.testing, occs.testing.z)
    }else{
      occs.testing.z <- occs.testing
    }
  }else{
    occs.testing.z <- NULL
  }
  
  ################################# #
  # ASSIGN CATEGORICAL VARIABLES ####
  ################################# #
  
  # find factor rasters or columns and identify them as categoricals
  if(!is.null(envs)) {
    categoricals <- unique(c(categoricals, names(envs)[which(raster::is.factor(envs))]))
  }else{
    categoricals <- unique(c(categoricals, names(occs)[which(sapply(occs, is.factor))]))
  }
  if(length(categoricals) == 0) categoricals <- NULL
  
  # if categoricals argument was specified, convert these columns to factor class
  if(!is.null(categoricals)) {
    for(i in 1:length(categoricals)) {
      if(quiet != TRUE) message(paste0("* Assigning variable ", categoricals[i], " to categorical ..."))
      d[, categoricals[i]] <- as.factor(d[, categoricals[i]])
      if(!is.null(user.val.grps)) user.val.grps[, categoricals[i]] <- factor(user.val.grps[, categoricals[i]], levels = levels(d[, categoricals[i]]))
      if(!is.null(occs.testing.z)) occs.testing.z[, categoricals[i]] <- factor(occs.testing.z[, categoricals[i]], levels = levels(d[, categoricals[i]]))
    }
  }
  
  # drop categoricals designation in other.settings to feed into other functions
  other.settings$categoricals <- categoricals
  
  ###################### #
  # ASSIGN PARTITIONS ####
  ###################### #
  
  # unpack occs and bg records for partitioning
  d.occs <- d %>% dplyr::filter(pb == 1) %>% dplyr::select(1:2)
  d.bg <- d %>% dplyr::filter(pb == 0) %>% dplyr::select(1:2)
  
  # partition occs based on selected partition method
  grps <- switch(partitions, 
                 jackknife = get.jackknife(d.occs, d.bg),
                 randomkfold = get.randomkfold(d.occs, d.bg, partition.settings$kfolds),
                 block = get.block(d.occs, d.bg, partition.settings$orientation),
                 checkerboard1 = get.checkerboard1(d.occs, envs, d.bg, partition.settings$aggregation.factor),
                 checkerboard2 = get.checkerboard2(d.occs, envs, d.bg, partition.settings$aggregation.factor),
                 user = NULL,
                 testing = list(occs.grp = rep(0, nrow(d.occs)), bg.grp = rep(0, nrow(d.bg))),
                 none = list(occs.grp = rep(0, nrow(d.occs)), bg.grp = rep(0, nrow(d.bg))))
  
  
  # choose a user message reporting on partition choice
  parts.message <- switch(partitions,
                          jackknife = "* Model evaluations with k-1 jackknife (leave-one-out) cross validation...",
                          randomkfold = paste0("* Model evaluations with random ", partition.settings$kfolds, "-fold cross validation..."),
                          block =  paste0("* Model evaluations with spatial block (4-fold) cross validation and ", partition.settings$orientation, " orientation..."),
                          checkerboard1 = "* Model evaluations with checkerboard (2-fold) cross validation...",
                          checkerboard2 = "* Model evaluations with hierarchical checkerboard (4-fold) cross validation...",
                          user = paste0("* Model evaluations with user-defined ", length(unique(user.grp$occs.grp)), "-fold cross validation..."),
                          testing = "* Model evaluations with testing data...",
                          none = "* Skipping model evaluations (only calculating full model statistics)...")
  if(quiet != TRUE) message(parts.message)
  
  # record user partition settings
  if(partitions == "user") partition.settings <- c(partition.settings, kfolds = length(unique(user.grp$occs.grp)))
  
  # if jackknife partitioning, do not calculate CBI because there are too few validation occurrences
  # per partition (n=1) to have a meaningful result
  if(partitions == "jackknife") other.settings$cbi.cv <- FALSE else other.settings$cbi.cv <- TRUE
  
  # add these values as the 'grp' column
  if(!is.null(grps)) d$grp <- factor(c(grps$occs.grp, grps$bg.grp))
  
  ################# #
  # MESSAGE
  ################# #
  # print model-specific message
  if(is.null(taxon.name)) {
    if(quiet != TRUE) message(paste("\n*** Running ENMeval v2.0.0 with", enm@msgs(tune.args, other.settings), "***\n"))
  }else{
    if(quiet != TRUE) message(paste("\n*** Running ENMeval v2.0.0 for", taxon.name, "with", enm@msgs(tune.args, other.settings), "***\n"))
  }
  
  ################# #
  # CLAMPING ####
  ################# #
  if(doClamp == TRUE) {
    if(is.null(envs)) {
      if(!quiet) message("Warning: cannot make model extrapolations (and therefore clamp) without predictor variable rasters.")
    }else{
      if(is.null(clamp.directions)) {
        clamp.envs <- names(envs)[!names(envs) %in% categoricals]
        clamp.directions <- list(left = clamp.envs, right = clamp.envs)
      }
      envs <- clamp.vars(predictors = envs, p.z = rbind(occs.z, bg.z), 
                         left = clamp.directions$left, right = clamp.directions$right, 
                         categoricals = categoricals)  
      if(quiet != TRUE) message("* Clamping predictor variable rasters...")
    }
  }
  
  ################# #
  # MODEL TUNING #### 
  ################# #
  
  # make table for all tuning parameter combinations
  tune.tbl <- expand.grid(tune.args, stringsAsFactors = FALSE) %>% tibble::as_tibble()
  # make tune.tbl NULL, not an empty table, if no settings are specified
  # this makes it easier to use tune.i as a parameter in function calls
  # when tune.args does not exist
  if(nrow(tune.tbl) == 0) tune.tbl <- NULL
  
  if(parallel) {
    results <- tune.parallel(d, envs, enm, partitions, tune.tbl, other.settings, partition.settings, 
                             user.val.grps, occs.testing.z, numCores, parallelType, user.eval, quiet)  
  }else{
    results <- tune.regular(d, envs, enm, partitions, tune.tbl, other.settings, partition.settings,
                            user.val.grps, occs.testing.z, updateProgress, user.eval, quiet)
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
  if(partitions %in% c("testing", "none")) {
    nk <- 0
  }else{
    nk <- length(unique(d[d$pb == 1, "grp"]))
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
    if(nrow(val.stats.all) > 0) eval.stats <- dplyr::left_join(eval.stats, val.stats.all, by = "tune.args")
    if("fold" %in% names(eval.stats)) eval.stats <- eval.stats %>% select(-fold)
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
  
  if(is.null(taxon.name)) taxon.name <- ""
  if(is.null(tune.tbl)) tune.tbl <- data.frame()
  if(is.null(occs.testing.z)) occs.testing.z <- data.frame()
  if(partitions != "block") partition.settings$orientation <- NULL
  if(partitions != "checkerboard1" | partitions != "checkerboard2") partition.settings$aggregation.factor <- NULL
  if(partitions != "randomkfold") partition.settings$kfolds <- NULL
  if(is.null(partition.settings) | length(partition.settings) == 0) partition.settings <- list()
  if(is.null(clamp.directions)) clamp.directions <- list()
  
  # get variable importance for all models
  varimp.all <- lapply(mod.full.all, enm@varimp)
  
  # assemble the ENMevaluation object
  e <- ENMevaluation(algorithm = enm@name, tune.settings = as.data.frame(tune.tbl),
                     results = as.data.frame(eval.stats), results.partitions = val.stats.all,
                     predictions = mod.full.pred.all, models = mod.full.all, 
                     variable.importance = varimp.all,
                     partition.method = partitions, partition.settings = partition.settings,
                     other.settings = other.settings, doClamp = doClamp, clamp.directions = clamp.directions, 
                     taxon.name = taxon.name,
                     occs = d[d$pb == 1, 1:(ncol(d)-2)], occs.testing = occs.testing.z, occs.grp = factor(d[d$pb == 1, "grp"]),
                     bg = d[d$pb == 0, 1:(ncol(d)-2)], bg.grp = factor(d[d$pb == 0, "grp"]),
                     rmm = list())
  
  # add the rangeModelMetadata object to the ENMevaluation object
  # write to existing RMM if input by user
  e@rmm <- buildRMM(e, envs, rmm)
  
  # if niche overlap selected, calculate and add the resulting matrix to results
  if(overlap == TRUE) {
    nr <- raster::nlayers(e@predictions)
    if(nr == 0) {
      if(quiet != TRUE) message("Warning: calculate niche overlap without model prediction rasters.")
    }else if(nr == 1) {
      if(quiet != TRUE) message("Warning: only 1 model prediction raster found. Need at least 2 rasters to calculate niche overlap. Increase number of tuning arguments and run again.") 
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
