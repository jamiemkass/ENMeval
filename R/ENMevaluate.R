#' @title Tune ecological niche model (ENM) settings and calculate evaluation statistics
#' @description \code{ENMevaluate()} is the primary function for the \pkg{ENMeval} package. This
#' function builds ecological niche models iteratively across a range of user-specified tuning 
#' settings. Users can choose to evaluate models with cross validation or a full-withheld testing 
#' dataset. \code{ENMevaluate()} returns an \code{ENMevaluation} object with slots containing 
#' evaluation statistics for each combination of settings and for each cross validation fold therein, as
#' well as raster predictions for each model when raster data is input. The evaluation statistics in the 
#' results table should aid users in identifying model settings that balance fit and predictive ability. See
#' the extensive vignette for fully worked examples: 
#' <https://jamiemkass.github.io/ENMeval/articles/ENMeval-2.0-vignette.html>.
#' 
#' @param occs matrix / data frame: occurrence records with two columns for longitude and latitude 
#' of occurrence localities, in that order. If specifying predictor variable values
#' assigned to presence/background localities (without inputting raster data), this table should also have 
#' one column for each predictor variable. See Note for important distinctions between running the function
#' with and without rasters.
#' @param envs SpatRaster: environmental predictor variables. These should be in same geographic projection as occurrence data.
#' @param bg matrix / data frame: background records with two columns for longitude and latitude of 
#' background (or pseudo-absence) localities, in that order. If NULL, points will be randomly sampled across \code{envs} 
#' with the number specified by argument \code{n.bg}. If specifying predictor variable values
#' assigned to presence/background localities (without inputting raster data), this table should also have 
#' one column for each predictor variable. See Details for important distinctions between running the function
#' with and without rasters.
#' @param tune.args named list: model settings to be tuned (i.e., for Maxent models:  \code{list(fc = c("L","Q"), rm = 1:3)})
#' @param partitions character: name of partitioning technique. Currently available options are
#' the nonspatial partitions "randomkfold" and "jackknife", the spatial partitions "block" and
#' "checkerboard", "testing" for partitioning with fully withheld data (see 
#' argument occs.testing), the "user" option (see argument user.grp), and "none" for no partitioning 
#' (see \code{?partitions} for details).
#' @param algorithm character: name of the algorithm used to build models. Currently one of "maxnet",
#' "maxent.jar", or "bioclim", else the name from a custom ENMdetails implementation.
#' @param partition.settings named list: used to specify certain settings for partitioning schema.
#' See Details and ?partitions for descriptions of these settings.
#' @param other.settings named list: used to specify extra settings for the analysis. 
#' All of these settings have internal defaults, so if they are not specified the analysis will be run 
#' with default settings. See Details for descriptions of these settings, including how to specify arguments
#' for maxent.jar.
#' @param categoricals character vector: name or names of categorical environmental variables. If not specified,
#' all predictor variables will be treated as continuous unless they are factors. If categorical variables
#' are already factors, specifying names of such variables in this argument is not needed.
#' @param doClamp boolean: if TRUE (default), model prediction extrapolations will be restricted to the upper and lower
#' bounds of the predictor variables. Clamping avoids extreme predictions for environment values outside
#' the range of the training data. If free extrapolation is a study aim, this should be set to FALSE, but
#' for most applications leaving this at the default of TRUE is advisable to avoid unrealistic predictions. 
#' When predictor variables are input, they are clamped internally before making model predictions when clamping is on.
#' When no predictor variables are input and data frames of coordinates and variable values are used instead (SWD format),
#' validation data is clamped before making model predictions when clamping is on.
#' @param raster.preds boolean: if TRUE (default), return model prediction rasters. If this is FALSE,
#' the predictions slot in the ENMevaluation object will be empty, which is the same as if no raster 
#' data is input. You can still make model prediction rasters using the model objects in the models slot
#' with the predict() function.
#' @param clamp.directions named list: specifies the direction ("left" for minimum, "right" for maximum) 
#' of clamping for predictor variables -- (e.g., \code{list(left = c("bio1","bio5"), right = c("bio10","bio15"))}).
#' @param user.enm ENMdetails object: a custom ENMdetails object used to build models. 
#' This is an alternative to specifying \code{algorithm} with a character string.
#' @param user.grp named list: specifies user-defined partition groups, where \code{occs.grp} = vector of partition group 
#' (fold) for each occurrence locality, intended for user-defined partitions, and \code{bg.grp} = same vector for 
#' background (or pseudo-absence) localities.
#' @param occs.testing matrix / data frame: a fully withheld testing dataset with two columns for longitude and latitude 
#' of occurrence localities, in that order when \code{partitions = "testing"}. These occurrences will be used only 
#' for evaluation but not for model training, and thus no cross validation will be performed.
#' @param n.bg numeric: the number of background (or pseudo-absence) points to randomly sample over the environmental  
#' raster data (default: 10000) if background records were not already provided.
#' @param overlap boolean: if TRUE, calculate range overlap statistics (Warren \emph{et al.} 2008).
#' @param overlapStat character: range overlap statistics to be calculated -- 
#' "D" (Schoener's D) and or "I" (Hellinger's I) -- see ?calc.niche.overlap for more details.
#' @param user.val.grps matrix / data frame: user-defined validation record coordinates and predictor variable values. 
#' This is used internally by \code{ENMnulls()} to force each null model to evaluate with empirical validation data,
#' and does not have any current use when running \code{ENMevaluate()} independently.
#' @param user.eval function: custom function for specifying performance metrics not included in \pkg{ENMeval}.
#' The function must first be defined and then input as the argument \code{user.eval}. 
#' This function should have a single argument called \code{vars}, which is a list that includes different data 
#' that can be used to calculate the metric. See Details below and the vignette for a worked example.
#' @param rmm rangeModelMetadata object: if specified, \code{ENMevaluate()} will write metadata details for the analysis into
#' this object, but if not, a new \code{rangeModelMetadata} object will be generated and included in the output
#' \code{ENMevaluation} object.
#' @param parallel boolean: if TRUE, run with parallel processing.
#' @param numCores numeric: number of cores to use for parallel processing. If NULL, all available cores will be used.
#' @param updateProgress boolean: if TRUE, use shiny progress bar. This is only for use in shiny apps.
#' @param quiet boolean: if TRUE, silence all function messages (but not errors).
#' 
#' @details There are a few methodological details in the implementation of ENMeval >=2.0.0 that are important to mention.
#' There is also a brief discussion of some points relevant to null models in ?ENMnulls.
#' 
#' 1. By default, validation AUC is calculated with respect to the full background (training + validation).
#' This approach follows Radosavljevic & Anderson (2014).This setting can be changed by assigning 
#' other.settings$validation.bg to "partition", which will calculate AUC with respect 
#' to the validation background only. The default value for other.settings$validation.bg is "full". 
#' NOTE: When examining validation AUC and other discrimination metrics, the "full" option will likely 
#' result in higher performance than for the "partition" option because more and varied background data 
#' should lead to higher discriminatory power for the model. Users should thus make sure they are correctly
#' interpreting the evaluation results.
#' 
#' 2. The continuous Boyce index (always) and AICc (when no raster is provided) are not calculated using 
#' the predicted values of the SpatRaster delineating the full study extent, but instead using the predicted
#' values for the background records. This decision to use the background only for calculating the continuous 
#' Boyce index was made to simplify the code and improve running time. The decision for AICc was made in order
#' to allow AICc calculations for datasets that do not include raster data. See ?calc.aicc for more details,
#' and for caveats when calculating AICc without raster data (mainly, that if the background does not 
#' adequately represent the occurrence records, users should use the raster approach, for reasons explained
#' in the calc.aicc documentation). For both metrics, if the background records are a good representation 
#' of the study extent, there should not be much difference between this approach using the background 
#' data and the approach that uses rasters.
#' 
#' 3. When running \code{ENMevaluate()} without raster data, and instead adding the environmental predictor values
#' to the occurrence and background data tables, users may notice some differences in the results. Occurrence records
#' that share a raster grid cell are automatically removed when raster data is provided, but without raster data
#' this functionality cannot operate, and thus any such duplicate occurrence records can remain in the training data.
#' The Java implementation of Maxent (maxent.jar implemented with \code{MaxEnt()} from the R package \code{predicts}) should automatically remove these records, but the R implementation 
#' \code{maxnet} does not, and the \code{envelope()} function from the R package \code{predicts} does not as well. Therefore,  
#' it is up to the user to remove such records before running \code{ENMevaluate()} when raster data are not included.
#' 
#' Below are descriptions of the parameters used in the other.settings, partition.settings, and user.eval arguments.
#' 
#' For other.settings, the options are:\cr*
#' path - character: the folder path designating where maxent.jar files should be saved\cr*
#' removeduplicates - boolean: whether or not to remove grid-cell duplicates for occurrences 
#' (this controls behavior for maxent.jar and ENMeval)\cr*
#' addsamplestobackground - boolean: whether or not to add occurrences to the background
#' when modeling with maxnet -- the default is TRUE.\cr*
#' abs.auc.diff - boolean: if TRUE, take absolute value of AUCdiff (default: TRUE)\cr*
#' pred.type - character: specifies which prediction type should be used to generate maxnet or 
#' maxent.jar prediction rasters (default: "cloglog").\cr*
#' validation.bg - character: either "full" to calculate training and validation AUC and CBI 
#' for cross-validation with respect to the full background (default), or "partition" (meant for 
#' spatial partitions only) to calculate each with respect to the partitioned background only 
#' (i.e., training occurrences are compared to training background, and validation occurrences 
#' compared to validation background).\cr*
#' other.args - named list: any additional model arguments not specified for tuning; this can
#' include arguments for maxent.jar, which are described in the software's Help file, such as
#' "jackknife=TRUE" for a variable importance jackknife plot or "responsecurves=TRUE" for response
#' curve plots -- note the the "path" must be specified (see above).\cr
#' 
#' For partition.settings, the current options are:\cr*
#' orientation - character: one of "lat_lon" (default), "lon_lat", "lat_lat", or "lon_lon" (required for block partition).\cr* 
#' aggregation.factor - numeric vector: one or two numbers specifying the factor with which to aggregate the envs (default: 2)
#' raster to assign partitions (required for the checkerboard partitions).\cr*
#' kfolds - numeric: the number of folds (i.e., partitions) for random partitions (default: 5).\cr
#' 
#' For the block partition, the orientation specifications are abbreviations for "latitude" and "longitude", 
#' and they determine the order and orientations with which the block partitioning function creates the partition groups. 
#' For example, "lat_lon" will split the occurrence localities first by latitude, then by longitude. For the checkerboard 
#' partitions, the aggregation factor specifies how much to aggregate the existing cells in the envs raster
#' to make new spatial partitions. For example, 'basic' checkerboard with an aggregation factor value of 2 will make squares  
#' 4 times larger than the input rasters and assign occurrence and background records to partition groups based on which square they fall in. 
#' Using two aggregation factors makes the checkerboard partitions hierarchical, where squares are first aggregated to define groups as in the 'basic' checkerboard, but a 
#' second aggregation is then made to separate the resulting two bins into four bins (see ?partitions for more details).
#' 
#' For user.eval, the variables you have access to in order to run your custom function are below. 
#' See the vignette for a worked example.\cr*
#' enm - ENMdetails object\cr*
#' occs.train.z - data frame: predictor variable values for training occurrences\cr*
#' occs.val.z - data frame: predictor variable values for validation occurrences\cr*
#' bg.train.z - data frame: predictor variable values for training background\cr*
#' bg.val.z - data frame: predictor variable values for validation background\cr*
#' mod.k - Model object for current partition (k)\cr*
#' nk - numeric: number of folds (i.e., partitions)\cr*
#' other.settings - named list: other settings specified in ENMevaluate()\cr*
#' partitions - character: name of the partition method (e.g., "block")\cr*
#' occs.train.pred - numeric: predictions made by mod.k for training occurrences\cr*
#' occs.val.pred - numeric: predictions made by mod.k for validation occurrences\cr*
#' bg.train.pred - numeric: predictions made by mod.k for training background\cr*
#' bg.val.pred - numeric: predictions made by mod.k for validation background
#' 
#' @references 
#' 
#' Muscarella, R., Galante, P. J., Soley-Guardia, M., Boria, R. A., Kass, J. M., Uriarte, M., & Anderson, R. P. (2014). ENMeval: An R package for conducting spatially independent evaluations and estimating optimal model complexity for Maxent ecological niche models. \emph{Methods in Ecology and Evolution}, \bold{5}: 1198-1205. \doi{10.1111/2041-210X.12261}
#' 
#' Warren, D. L., Glor, R. E., Turelli, M. & Funk, D. (2008) Environmental niche equivalency versus conservatism: quantitative approaches to niche evolution. \emph{Evolution}, \bold{62}: 2868-2883. \doi{10.1111/j.1558-5646.2008.00482.x}
#' 
#' @return An ENMevaluation object. See ?ENMevaluation for details and description of the columns
#' in the results table.
#'
#' @importFrom foreach %dopar%
#' @importFrom grDevices rainbow
#' @importFrom methods new slot validObject
#' @importFrom stats pnorm quantile runif sd quantile
#' @importFrom utils citation combn packageVersion setTxtProgressBar txtProgressBar
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(ENMeval)
#' 
#' occs <- read.csv(file.path(system.file(package="predicts"), 
#' "/ex/bradypus.csv"))[,2:3]
#' envs <- rast(list.files(path=paste(system.file(package="predicts"), 
#' "/ex", sep=""), pattern="tif$", full.names=TRUE))
#' occs.xyEnv <- cbind(occs, extract(envs, occs, ID = FALSE))
#' occs.xyEnv$biome <- factor(occs.xyEnv$biome)
#' bg <- as.data.frame(predicts::backgroundSample(envs, n = 10000))
#' names(bg) <- names(occs)
#' bg.xyEnv <- cbind(bg, extract(envs, bg, ID = FALSE))
#' bg.xyEnv$biome <- factor(bg.xyEnv$biome)
#' 
#' # set other.settings -- pred.type is only for Maxent models
#' os <- list(abs.auc.diff = FALSE, pred.type = "cloglog", 
#' validation.bg = "partition")
#' # set partition.settings -- here's an example for the block method
#' # see Details for the required settings for other partition methods
#' ps <- list(orientation = "lat_lat")
#' 
#' # here's a run with maxnet -- note the tune.args for feature classes (fc)
#' # and regularization multipliers (rm), as well as the designation of the
#' # categorical variable we are using (this can be a vector if multiple
#' # categorical variables are used)
#' e.maxnet <- ENMevaluate(occs, envs, bg, 
#' tune.args = list(fc = c("L","LQ","LQH","H"), rm = 1:5), 
#' partitions = "block", other.settings = os, partition.settings = ps,
#' algorithm = "maxnet", categoricals = "biome", overlap = TRUE)

#' # you can run the same model without rasters using the xyEnv data
#' e.maxnet <- ENMevaluate(occs.xyEnv, envs = NULL, bg.xyEnv, 
#' tune.args = list(fc = c("L","LQ","LQH","H"), rm = 1:5), 
#' partitions = "block", other.settings = os, partition.settings = ps,
#' algorithm = "maxnet", categoricals = "biome", overlap = TRUE)
#' 
#' # print the tuning results
#' eval.results(e.maxnet)
#' 
#' # you can plot the marginal response curves of a maxnet object with plot(), 
#' # and you can also extract the data for plotting to make your own custom plots
#' mods.maxnet <- eval.models(e.maxnet)
#' m <- mods.maxnet$fc.LQH_rm.2
#' plot(m, type = "cloglog")
#' rcurve_data <- maxnet::response.plot(m, "bio1", type = "cloglog", plot = FALSE)
#' 
#' # there is currently no native function to make raster model predictions for
#' # maxnet models, but ENMeval can be used to make them like this:
#' # here's an example where we make a prediction based on the L2 model
#' # (feature class: Linear, regularization multiplier: 2) for our envs data
#' pred.LQH2 <- maxnet.predictRaster(m, envs)
#' plot(pred.L2)
#' 
#' # here's a run with maxent.jar -- note that if the R package rJava cannot 
#' # install or load, or if you have other issues with Java on your computer, 
#' # maxent.jar will not function
#' e.maxent.jar <- ENMevaluate(occs, envs, bg, 
#' tune.args = list(fc = c("L","LQ","LQH","H"), rm = 1:5), 
#' partitions = "block", other.settings = os, partition.settings = ps,
#' algorithm = "maxent.jar", categoricals = "biome", overlap = TRUE)
#' 
#' # here's a run of maxent.jar with a path specified for saving the html and 
#' # plot files -- you can also turn on jackknife variable importance or 
#' # response curves, etc., to have these plots saved there
#' e.maxent.jar <- ENMevaluate(occs, envs, bg, 
#' tune.args = list(fc = c("L","LQ","LQH","H"), rm = 1:5), 
#' partitions = "block", partition.settings = ps,
#' algorithm = "maxent.jar", categoricals = "biome", overlap = TRUE,
#' other.settings = list(path = "analyses/mxnt_results", 
#' other.args = c("jackknife=TRUE", "responsecurves=TRUE")))
#' 
#' # print the tuning results
#' eval.results(e.maxent.jar)
#' 
#' # raster predictions can be made for maxent.jar models with predicts or 
#' # ENMeval
#' mods.maxent.jar <- eval.models(e.maxent.jar)
#' pred.L2 <- predict(mods.maxent.jar$fc.L_rm.2, envs, 
#' args = "outputform=cloglog")
#' pred.L2 <- maxnet.predictRaster(mods.maxent.jar$fc.L_rm.2, envs, os)
#' plot(pred.L2)
#' 
#' # this will give you the percent contribution (not deterministic) and
#' # permutation importance (deterministic) values of variable importance for
#' # Maxent models, and it only works with maxent.jar
#' eval.variable.importance(e.maxent.jar)
#' 
#' # here's a run with BIOCLIM. Note that we need to remove the categorical
#' # variable here because this algorithm only takes continuous variables. We 
#' # also should point out that the way BIOCLIM is tuned is by comparing 
#' # performance for different ways to make predictions (as opposed to comparing 
#' # performance for models fit in different ways like for maxnet or maxent.jar). 
#' # Namely, BIOCLIM can ignore different tails of the distribution when making 
#' # predictions, and this is what is tuned in ENMevaluate (see 
#' # ?predicts::envelope).
# e.bioclim <- ENMevaluate(occs, envs[[-9]], bg,
# tune.args = list(tails = c("low", "high", "both")),
# partitions = "block", other.settings = os, partition.settings = ps,
# algorithm = "bioclim", overlap = TRUE)
#' 
#' # print the tuning results
#' eval.results(e.bioclim)
#' # make raster predictions with predicts or ENMeval
#' mods.bioclim <- eval.models(e.bioclim)
#' # note: the models for low, high, and both are actually all the same, and
#' # the only difference for tuning is how they are predicted during
#' # cross-validation
#' pred.both <- predict(mods.bioclim$tails.both, envs, tails = "both")
#' plot(pred.both)
#' 
#' # please see the vignette for more examples of model tuning, 
#' # partitioning, plotting functions, and null models
#' # https://jamiemkass.github.io/ENMeval/articles/ENMeval-2.0-vignette.html
#' }
#' 
#' @export 

ENMevaluate <- function(occs, envs = NULL, bg = NULL, tune.args = NULL, 
                        partitions = NULL, algorithm = NULL, 
                        partition.settings = NULL, 
                        other.settings = list(), categoricals = NULL, 
                        doClamp = TRUE, raster.preds = TRUE,
                        clamp.directions = NULL, user.enm = NULL, 
                        user.grp = NULL, occs.testing = NULL, 
                        n.bg = 10000, overlap = FALSE, 
                        overlapStat = c("D", "I"), user.val.grps = NULL, 
                        user.eval = NULL, rmm = NULL, parallel = FALSE, 
                        numCores = NULL, updateProgress = FALSE) {
  
  
  # record start time
  start.time <- proc.time()
  
  # check if ecospat is installed, and if not, prevent CBI calculations
  if((requireNamespace("ecospat", quietly = TRUE))) {
    ecospat.use <- TRUE
  }else{
    message("Package ecospat is not installed, so Continuous Boyce Index (CBI) cannot be calculated.")
    ecospat.use <- FALSE
  }
  # add whether to use ecospat to other.settings to avoid multiple requires
  other.settings <- c(other.settings, ecospat.use = ecospat.use)
  
  ######################## #
  # INITIAL DATA CHECKS ####
  ######################## #
  
  # coerce occs and bg to df
  occs <- as.data.frame(occs)
  if(!is.null(bg)) bg <- as.data.frame(bg)
  
  # add defaults to other.settings if not entered by user
  other.settings <- checks.other.settings(other.settings)
  
  # make sure taxon name column is not included
  if(inherits(occs[,1], "character") | inherits(bg[,1], "character")) 
    stop("* The two first columns of occs should be longitude and latitude.")
  
  # initial message
  message("*** Running initial checks... ***\n")
  
  # common partitions input errors
  checks.partitions(partitions, partition.settings, occs.testing)
  
  # incompatibilities between envs and algorithms
  checks.envs(envs, algorithm)
  
  # if a vector of tuning arguments is numeric, make sure it is sorted 
  # (for results table and plotting)
  tune.args <- checks.tuning(tune.args)
  
  # turn off overlap is no tuning specified  
  if(is.null(tune.args)) {
    overlap <- FALSE
  }
  
  # make sure occs and bg are data frames with identical column names
  if(all(names(occs) != names(bg)) & !is.null(bg)) {
    stop('Datasets "occs" and "bg" have different column names. Please make them identical and try again.')
  }
  
  if(doClamp == FALSE & !is.null(clamp.directions)) {
    stop("If specifying clamp directions, please make doClamp = TRUE.")
  }
  
  # make sure that validation.bg is only "bg" if the partitions are spatial
  if(!(partitions %in% c("block","checkerboard","user")) & 
     other.settings$validation.bg == "partition") {
    stop('If using non-spatial partitions, please set validation.bg to "full". The "partition" option only makes sense when partitions represent different regions of the study extent. See ?ENMevaluate for details.')
  }else if(partitions == "user" & other.settings$validation.bg == "partition") {
    message('* Please make sure that the user-specified partitions are spatial, else validation.bg should be set to "full". The "partition" option only makes sense when partitions represent different regions of the study extent. See ?ENMevaluate for details.')
  }
  
  # choose a built-in ENMdetails object matching the input model name
  # unless the model is chosen by the user
  if(is.null(user.enm)) {
    enm <- lookup.enm(algorithm)
  }else{
    enm <- user.enm
  }
  
  # get unique values of user partition groups to make sure they all remain 
  # after occurrence data processing
  if(!is.null(user.grp)) {
    user.nk <- length(unique(c(user.grp$occs.grp, user.grp$bg.grp)))
    partition.settings$kfolds <- user.nk
  }
  
  # algorithm-specific errors
  enm@errors(occs, envs, bg, tune.args, partitions, algorithm, 
             partition.settings, other.settings, 
             categoricals, doClamp, clamp.directions)
  
  
  ############################################################## #
  # ASSEMBLE OCCS / BG DATA AND ENVIRONMENTAL VARIABLE VALUES ####
  ############################################################## #
  
  # keep records of which points were removed during cleaning
  occs.rem <- numeric(nrow(occs))
  
  # if environmental rasters are input as predictor variables
  if(!is.null(envs)) {
    # if no background points specified, generate random ones
    if(is.null(bg)) {
      message(paste0("* Randomly sampling ", n.bg, " background points ..."))
      bg <- terra::spatSample(envs, size = n.bg, xy = TRUE, na.rm = TRUE,
                              values = FALSE) |> as.data.frame()
      names(bg) <- names(occs)
    }
    
    # flag occurrence record cell duplicates
    if(other.settings$removeduplicates == TRUE) {
      # extract cell number values from occs to determine duplicates
      occs.cellNo <- terra::extract(envs, occs, cells = TRUE, ID = FALSE)
      occs.dups <- which(duplicated(occs.cellNo[,"cell"]) == 1)
      if(sum(occs.dups) > 0) {
        occs.rem[occs.dups] <- 1
        message(paste0("* Removing ", length(occs.dups), " occurrence record grid cell duplicates."))
      }
    }
    
    # extract env values from occs and bg
    occs.z <- terra::extract(envs, occs, ID = FALSE)  
    bg.z <- terra::extract(envs, bg, ID = FALSE)
    
    # bind coordinates to predictor variable values for occs and bg
    # occs <- cbind(occs, occs.z)
    # bg <- cbind(bg, bg.z)
    
    # if no environmental rasters input (SWD)
  }else{
    # assign occs/bg to coordinates and occs.z/bg.z to variable values
    occs.z <- occs[,-1:-2]
    occs <- occs[,1:2]
    bg.z <- bg[,-1:-2]
    bg <- bg[,1:2]
    # if no bg included, cannot continue because no raster for random sampling
    if(is.null(bg)) stop("* If using species with data (SWD) format, input coordinates for both occurrence and background records.")
    # make sure both occ and bg have predictor variable values
    if(ncol(occs) < 3 | ncol(bg) < 3) stop("* If using species with data (SWD) format, input variable values in occs and bg tables proceeding the coordinates.")
    message("* Using species with data (SWD) format, so raster predictions cannot be generated. Also, AICc is calculated with background data for Maxent models.")
  }
  
  # keep records of which points were removed during cleaning
  bg.rem <- numeric(nrow(bg))
  
  # if NA predictor variable values exist for occs or bg, flag these records 
  occs.na <- which(rowSums(is.na(occs.z)) > 0)
  if(length(occs.na) > 0) {
    # which many occs are newly filtered?
    occs.na.rem <- occs.na[which(occs.rem[occs.na] == 0)]
    message(paste0("* Removing ", length(occs.na.rem), " occurrence records with NA predictor variable values."))
    occs.rem[occs.na.rem] <- 1
  }
  bg.na <- which(rowSums(is.na(bg.z)) > 0)
  if(length(bg.na) > 0) {
    message(paste0("* Removing ", length(bg.na), " background records with NA predictor variable values."))
    bg.rem[bg.na] <- 1
  }
  
  # remove flagged occs and bg records
  occs <- occs[which(occs.rem == 0),]
  occs.z <- occs.z[which(occs.rem == 0),]
  bg <- bg[which(bg.rem == 0),]
  bg.z <- bg.z[which(bg.rem == 0),]
  
  # take out removed occs and bg from user.grp
  if(!is.null(user.grp)) {
    user.grp$occs.grp <- user.grp$occs.grp[which(occs.rem == 0)]
    user.nk.occs <- length(unique(user.grp$occs.grp))
    if(user.nk > user.nk.occs) stop("Occurrence record cleaning removed one or more user partition groups.")
    user.grp$bg.grp <- user.grp$bg.grp[which(bg.rem == 0)]
    user.nk.bg <- length(unique(user.grp$bg.grp))
    if(user.nk > user.nk.bg) stop("Background record cleaning removed one or more user partition groups.")
  }
  
  ############################################ #
  # ADD TESTING DATA (IF INPUT) ####
  ############################################ #
  
  # add testing data to main df if provided
  if(partitions == "testing") {
    if(!is.null(envs)) {
      occs.testing.z <- as.data.frame(terra::extract(envs, occs.testing, ID = FALSE))
    }else{
      occs.testing.z <- occs.testing[,-1:-2]
    }
  }else{
    occs.testing.z <- NULL
  }
  
  ################################# #
  # CHECK CATEGORICAL VARIABLES ####
  ################################# #
  
  checks.cats(occs.z, envs, categoricals)
  
  # put categoricals designation in other.settings to feed into other functions
  other.settings$categoricals <- categoricals
  
  ################# #
  # CLAMPING ####
  ################# #
  if(doClamp == TRUE) {
    # when predictor variable rasters are input
    if(!is.null(envs)) {
      # assign both clamp directions to all variables if none are set
      if(is.null(clamp.directions)) {
        clamp.directions$left <- names(envs)
        clamp.directions$right <- names(envs)
      }
      # record clamp directions in other.settings
      other.settings$clamp.directions <- clamp.directions
      # run function to clamp predictor variable rasters
      envs <- clamp.vars(orig.vals = envs, ref.vals = rbind(occs.z, bg.z), 
                         left = clamp.directions$left, 
                         right = clamp.directions$right, 
                         categoricals = categoricals)
      message("* Clamping predictor variable rasters...")
    }else{
      # if no predictor variable rasters are input, assign both clamp directions
      # to all variable names (columns besides the first two, which should be
      # coordinates) if none are set
      # during cross-validation, validation data will be clamped using the
      # clamp.vars() function
      if(is.null(clamp.directions)) {
        clamp.directions$left <- names(occs.z)
        clamp.directions$right <- names(occs.z)
      }
    }
  }
  # record clamping choice as FALSE in other.settings regardless of doClamp
  # selection -- this is because the clamping is done on the predictor rasters 
  # or values directly, so that internally all model predictions are made with 
  # clamping off
  # when enm.predict() functions are run externally, users can specify
  # other.settings$doClamp to turn on clamping functionality
  other.settings$doClamp <- FALSE
  
  ###################### #
  # ASSIGN PARTITIONS ####
  ###################### #
  
  # partition occs / bg based on selected partition method
  grps <- assign.partitions(occs, bg, envs, partitions, partition.settings)
  
  # choose a user message reporting on partition choice
  message(assign.partition.msg(partitions, partition.settings))
  
  # do not calculate CBI if one or more partitions has only 1 occurrence record
  if(partitions == "jackknife") {
    other.settings$cbi.cv <- FALSE
    message("For jackknife partitioning, not calculating CBI because each partition has too few occurrences for this metric.")
  }
  if(partitions == "user") {
    tb <- table(user.grp$occs.grp)
    user.grp.1 <- names(tb)[which(tb == 1)]
    if(length(user.grp.1 > 0)) other.settings$cbi.cv <- FALSE
    message("For user partitioning, not calculating CBI because partition ", user.grp.1, 
        "has too few occurrences for this metric.")
  }
  
  ################# #
  # MESSAGE
  ################# #
  message(paste("\n*** Running ENMeval v2.0.5 with", enm@msgs(tune.args, other.settings), "***\n"))
  
  ################# #
  # MODEL TUNING #### 
  ################# #
  
  # make table for all tuning parameter combinations
  tune.tbl <- expand.grid(tune.args, stringsAsFactors = FALSE) |> tibble::as_tibble()
  # make tune.tbl NULL, not an empty table, if no settings are specified
  # this makes it easier to use tune.i as a parameter in function calls
  # when tune.args does not exist
  if(nrow(tune.tbl) == 0) tune.tbl <- NULL
  
  # run parallel or regular
  if(parallel == TRUE) {
    results <- tuneParallel(occs.z, bg.z, grps, enm, partitions, tune.tbl, 
                            doClamp, other.settings, partition.settings, 
                            user.val.grps, occs.testing.z, numCores, user.eval, 
                            algorithm, updateProgress)
  }else{
    results <- tune(occs.z, bg.z, grps, enm, partitions, tune.tbl, doClamp, other.settings, 
                    partition.settings, user.val.grps, occs.testing.z, 
                    numCores, user.eval, algorithm, updateProgress)
  }
  
  ##################### #
  # ASSEMBLE RESULTS #### 
  ##################### #
  
  # flatten all training statistics data frames from results list into a single data frame
  train.stats.all <- dplyr::bind_rows(lapply(results, function(x) x$train.stats))
  # flatten all validation statistics data frames from results list into a single data frame
  # (these are no validation stats if no partitions were chosen)
  val.stats.all <- dplyr::bind_rows(lapply(results, function(x) x$cv.stats))
  
  if(!is.null(tune.tbl)) {
    # define tuned settings names and bind them to the tune table
    tune.names <- train.stats.all$tune.args
    tune.tbl$tune.args <- factor(tune.names, levels = tune.names)
  }else{
    # if not tuned settings, the "tune name" is the model name
    tune.names <- enm@name
  }
  # gather all full models into list and name them
  mod.full.all <- lapply(results, function(x) x$mod.full)
  names(mod.full.all) <- tune.names
  
  # gather all model prediction rasters into a stack and name them
  # if envs is null, make an empty stack
  if(!is.null(envs) & raster.preds == TRUE) {
    f <- function(x) enm@predict(x$mod.full, envs, other.settings)
    message("Making model prediction rasters...")
    mod.full.pred.all <- terra::rast(lapply(results, f))
    names(mod.full.pred.all) <- tune.names
  }else{
    mod.full.pred.all <- terra::rast()
  }
  
  ## assemble model evaluation statistics table
  if(partitions == "none") {
    # if no partitions were specified
    eval.stats <- assemble.stats.noParts(train.stats.all, val.stats.all, tune.names, tune.tbl)
  }else{
    # if partitions were specified
    nk <- length(unique(grps$occs.grp))
    eval.stats <- assemble.stats.parts(nk, train.stats.all, val.stats.all, tune.names, tune.tbl)
  }
  
  # calculate number of non-zero parameters in model
  ncoefs <- sapply(mod.full.all, enm@ncoefs)
  
  # calculate AICc and bind to eval.stats
  # currently, only works for Maxent models
  aic <- calcAIC(occs.z, envs, enm, mod.full.all, ncoefs, other.settings)
  eval.stats <- dplyr::bind_cols(eval.stats, aic)
  
  # add ncoef column
  eval.stats$ncoef <- ncoefs
  
  # clean up variables for making ENMevaluation object
  if(is.null(tune.tbl)) tune.tbl <- data.frame()
  if(is.null(occs.testing.z)) occs.testing.z <- data.frame()
  if(partitions != "block") partition.settings$orientation <- NULL
  if(partitions != "checkerboard") partition.settings$aggregation.factor <- NULL
  if(partitions != "randomkfold") partition.settings$kfolds <- NULL
  if(is.null(partition.settings) | length(partition.settings) == 0) partition.settings <- list()
  if(is.null(clamp.directions)) clamp.directions <- list()
  
  # get variable importance for all models
  variable.importance.all <- lapply(mod.full.all, enm@variable.importance)
  
  # overwrite the doClamp = FALSE used internally with the one entered by user
  other.settings$doClamp <- doClamp
  
  # if range overlap selected, calculate and add the resulting matrix to results
  if(overlap == TRUE) {
    ov.ls <- assemble.overlap(mod.full.pred.all, overlapStat)
  }else{
    ov.ls <- NULL
  }
  
  # assemble the ENMevaluation object
  e <- ENMevaluation(algorithm = enm@name, tune.settings = as.data.frame(tune.tbl),
                     results = as.data.frame(eval.stats), 
                     results.partitions = val.stats.all,
                     predictions = mod.full.pred.all, models = mod.full.all, 
                     variable.importance = variable.importance.all,
                     partition.method = partitions, 
                     partition.settings = partition.settings,
                     other.settings = other.settings, doClamp = doClamp, 
                     clamp.directions = clamp.directions, occs = occs.z, 
                     occs.testing = occs.testing.z, 
                     occs.grp = factor(grps$occs.grp), bg = bg.z, 
                     bg.grp = factor(grps$bg.grp), overlap = ov.ls, rmm = list())
  
  # add the rangeModelMetadata object to the ENMevaluation object
  # write to existing RMM if input by user
  e@rmm <- buildRMM(e, envs, rmm)
  
  
  # calculate time expended and print message
  timed <- proc.time() - start.time
  t.min <- floor(timed[3] / 60)
  t.sec <- timed[3] - (t.min * 60)
  message(paste("ENMevaluate completed in", t.min, "minutes", round(t.sec, 1), "seconds."))
  
  return(e)
}
