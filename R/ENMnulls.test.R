#' @title Compare model accuracy metrics of Ecological Niche Models (ENMs) built with different set of predictors.
#' #' @description \code{ENMnulls.test()} Iteratively builds null ENMs for "k" sets of user-specified model
#' settings bsaed on "k" input ENMevaluation objects, from which all other analysis 
#' settings are extracted.Summary statistics of the performance metrics for the null ENMs 
#' are taken (averages and standard deviations) and effect sizes and p-values are calculated by 
#' comparing these summary statistics to the empirical values of the performance metrics 
#' (i.e., from the model built with the empirical data). This is an extension of {ENMnulls()} for comparisons
#' of two or more predictor variable sets.
#' 
#' @param e.list: list of ENMevaluation objects to be compared
#' @param mod.settings.list named list: model settings corresponding to ENMevaluation objects in e.list 
#' that specify the settings used for building null models
#' @param no.iter numeric: number of null model iterations.
#' @param eval.stats character: model accuarcy metrics that will be used to calculate null model statistics. 
#' Can be “auc.val”, “auc.diff”, “cbi.val”, “or.mtp”, “or.10p”.
#' @param user.enm ENMdetails object: if implementing a user-specified model. 
#' @user.eval.type character: if implementing a user-specified model -- either "knonspatial", "kspatial", 
#' "testing" or "none".
#' @param userStats.signs named list: user-defined evaluation statistics attributed with either 1 or -1 
#' to designate whether the expected difference between empirical and null models is positive or negative; 
#' this is used to calculate the p-value of the z-score when comparing two predictor variable sets. Default is NULL.
#' @param removeMxTemp boolean: if TRUE, delete all temporary data generated when using maxent.jar for modeling
#' @param parallel boolean: if TRUE, use parallel processing.
#' @param numCores numeric: number of cores to use for parallel processing; if NULL, all available cores will be used.
#' @param parallelType character: either "doParallel" or "doSNOW" (default: "doSNOW").
#' @param quiet boolean: if TRUE, silence all function messages (but not errors).
#' 
#' #' @details This null ENM technique extends the implementation in Bohl \emph{et al.} (2019)and Kass \emph{et al.} (2020),
#' which follows the original methodology of Raes & ter Steege (2007). Here we evaluate if observed differences in accuracy metric values 
#' (e.g., AUC, omission rates, CBI) of empirical models built with different sets of predictor variable are greater than expected 
#' at random. This is done by building the null distributions of the difference in accuracy metrics
#  employing the same withheld validation data used to evaluate the empirical models. Please see the vignette for a brief example: <
#' 
#' This function avoids using raster data to speed up each iteration, and instead samples null occurrences from the 
#' partitioned background records. Thus, you should avoid running this when your background records are not well 
#' sampled across the study extent, as this limits the extent that null occurrences can be sampled from.
#' 
#' @references 
#' Bohl, C. L., Kass, J. M., & Anderson, R. P. (2019). A new null model approach to quantify performance and significance for ecological niche models of species distributions. \emph{Journal of Biogeography}, \bold{46}: 1101-1111. \url{https://doi.org/10.1111/jbi.13573}
#' 
#' Kass, J. M., Anderson, R. P., Espinosa-Lucas, A., Juárez-Jaimes, V., Martínez-Salas, E., Botello, F.,  Tavera, G., Flores-Martínez, J. J., & Sánchez-Cordero, V. (2020). Biotic predictors with phenological information improve range estimates for migrating monarch butterflies in Mexico. \emph{Ecography}, \bold{43}: 341-352. \url{https://doi.org/10.1111/ecog.04886}
#'
#' Raes, N., & ter Steege, H. (2007). A null-model for significance testing of presence-only species distribution models. \emph{Ecography}, \bold{30}: 727-736. \url{https://doi.org/10.1111/j.2007.0906-7590.05041.x} 
#' 
#' #' @return An \code{ENMnull}An ENMnull object with slots containing evaluation summary statistics 
#' for the null models and their cross-validation results, as well as differences in results between the 
#' empirical and null models. This comparison table includes T-statistics for pairwise comparisons (T-test) 
#' and F-statistic (ANOVA) of these differences and their associated p-values (under a normal distribution). 

ENMnulls.test <- 	function(e.list, mod.settings.list, no.iter, 
                           eval.stats = c("auc.val", "auc.diff","cbi.val","or.mtp","or.10p"),
                           user.enm = NULL, user.eval.type = NULL, userStats.signs = NULL,
                           removeMxTemp = TRUE, parallel = FALSE, numCores = NULL, 
                           parallelType = "doSnow", quiet = FALSE){
  #loading dependencies
  require(doParallel)
  require(doSNOW)
  ## checks
  #more than one input enm
  if(length(e.list) == 1){
    stop("Please input more than one ENM to run comparisons.")
  }
  
  # model settings are all single entries for each enm treatment
  for(k in 1:length(mod.settings.list)){
    if(!all(sapply(mod.settings.list[[k]], length) == 1)){
      stop("Please input a single set of model settings.")
    }
  }
  # model settings are correct for input algorithm and are entered in the right order -- 
  #if not, put them in the right order, else indexing models later will fail because the model 
  # name will be incorrect
  for(k in 1:length(e.list)){
    if(e.list[[k]]@algorithm %in% c("maxent.jar", "maxnet")){
        if(length(mod.settings.list[[k]]) != 2){
          stop("Please input two complexity settings (fc [feature classes] and rm [regularization
           multipliers]) for mod.settings for maxent.jar and maxnet models.")
        }
      if(all(names(mod.settings.list[[k]]) %in% c("fc", "rm"))) {
        if(!all(names(mod.settings.list[[k]]) == c("fc", "rm"))) {
          mod.settings.list[[k]] <- mod.settings.list[[k]][c("fc", "rm")]
        }
      }else{
        stop('Please input only "fc" (feature classes) and "rm" (regularization multipliers) for
           mod.settings for maxent.jar and maxnet models.')
        
      }
  }else if(e.list[[1]]@algorithm == "bioclim") {
    if(length(mod.settings.list[[k]]) != 1) {
      stop("Please input one complexity setting (tails) for mod.settings for BIOCLIM models.")
    }
    if(!all(names(mod.settings.list[[k]]) == "tails")) {
      stop('Please input only "tails" for mod.settings for BIOCLIM models.')
    }
    }
  }
  
    # assign evaluation type based on partition method
    if(is.null(user.eval.type)) {
      eval.type <- switch(e.list[[1]]@partition.method,
                          randomkfold = "knonspatial",
                          jackknife = "knonspatial",
                          block = "kspatial",
                          checkerboard1 = "kspatial",
                          checkerboard2 = "kspatial",
                          testing = "testing",
                          none = "none")  
    }else{
      eval.type <- user.eval.type
    }
  # assign directionality of sign for evaluation stats in post-hoc tests
  if(!is.null(userStats.signs)){
  signs <- c(list("auc.val" = 1, "auc.train" = 1, "cbi.val" = 1, "cbi.train" = 1,
                  "auc.diff" = -1, "or.10p" = -1, "or.mtp" = -1), userStats.signs)
  }
  
  # record start time
  start.time <- proc.time()
  
  ##############################
  ## 1. Create null occurrences
  ##############################
  #Create list of "i" sets of null "n" occurrences for each ENMevaluation object in e.list
  
  # assign the number of cross validation iterations
  # each ENMevaluation object in e.list should be built using the same parameters
  nk <- max(as.numeric(as.character(e.list[[1]]@occs.grp)))
  
  # get number of occurrence points by partition
  occs.grp.tbl.list <- lapply(e.list, function(e){table(e@occs.grp)})
  
  # if more than one background partition exists, assume spatial CV and
  # keep existing partitions
  
  #Create list of background SWD dataframes for each ENM treatment
  null.samps.list <- lapply(e.list, function(e){
    #Get occ points env values from ENMevaluation object
    null.samps <- cbind(rbind(e@occs, e@bg), grp = c(e@occs.grp, e@bg.grp))
    #Get bg points env values from ENMevaluation object
    return(null.samps)
  })
  
  for(k in 1:length(e.list)){
    if(e.list[[k]]@algorithm == "maxent.jar") {
      # create temp directory to store maxent.jar output, for potential removal later
      tmpdir <- paste(tempdir(), runif(1,0,1), sep = "/")
      dir.create(tmpdir, showWarnings = TRUE, recursive = FALSE)
    }  
  }
  
  # assign user algorithm if provided
  if(!is.null(user.enm)) {
    for(k in 1:length(e.list)){
      e.list[[k]]@algorithm <- user.enm
    }
  }
  
  #########################################
  ## 2. Specify empirical model statistics 
  #########################################
  # Specify empirical model statistics for each ENM treatment
  
  emp.mod.res.list<- mapply(function(e, m){
    mod.tune.args <- paste(names(m), m, collapse = "_", sep = ".")
    emp.mod <- e@models[[mod.tune.args]]
    emp.mod.res <- e@results %>% dplyr::filter(tune.args == mod.tune.args)
  },e.list, mod.settings.list)
  
  #########################################
  ## 3. Build null models
  #########################################           
  #Iteratively apply ENMulls for each ENMevaluation object in e.list
  
  if(quiet == FALSE) message(paste("Building and evaluating null ENMs with", no.iter, "iterations..."))
  if(quiet == FALSE) pb <- txtProgressBar(0, no.iter, style = 3)
  
  # set up parallel processing functionality
  if(parallel == TRUE) {
    allCores <- parallel::detectCores()
    if (is.null(numCores)) {
      numCores <- allCores
    }
    cl <- parallel::makeCluster(numCores, setup_strategy = "sequential")
    if(quiet != TRUE) progress <- function(n) setTxtProgressBar(pb, n)  
    
    if(parallelType == "doParallel") {
      doParallel::registerDoParallel(cl)
      opts <- NULL
    } else if(parallelType == "doSNOW") {
      doSNOW::registerDoSNOW(cl)
      if(quiet != TRUE) opts <- list(progress=progress) else opts <- NULL
    }
    numCoresUsed <- foreach::getDoParWorkers()
    if(quiet != TRUE) message(paste0("\nOf ", allCores, " total cores using ", numCoresUsed, "..."))
    if(quiet != TRUE) message(paste0("Running in parallel using ", parallelType, "..."))  
  }
  
  #Specify clamping directions
  clamp.directions.list <- lapply(e.list, function(e){
    if(length(e@clamp.directions) == 0){
      clamp.directions.i <- NULL
    }else{
      clamp.directions.i <- e@clamp.directions  
    }
  })
  
  
  # define function to run null model for iteration i
  null_i <- function(i) {
    null.occs.ik <- list()
    if(eval.type == "kspatial") {
      # randomly sample the same number of training occs over each k partition
      # of envs; if kspatial evaluation, only sample over the current spatial
      # partition of envs.z
      for(k in 1:nk) {
        # sample null occurrences only from
        # the records in partition group k
        null.samps.k <- null.samps %>% dplyr::filter(grp == k)
        # randomly sample n null occurrences, where n equals the number
        # of empirical occurrence in partition group k
        samp.k <- sample(1:nrow(null.samps.k), occs.grp.tbl[k])
        null.occs.ik[[k]] <- null.samps.k[samp.k, ]
      }
    }else if(eval.type == "knonspatial") {
      for(k in 1:nk) {
        # randomly sample n null occurrences, where n equals the number
        # of empirical occurrence in partition group k
        samp.k <- sample(1:nrow(null.samps), occs.grp.tbl[k])
        null.occs.ik[[k]] <- null.samps[samp.k, ]
      }
    }else if(eval.type %in% c("testing", "none")) {
      samp.test <- sample(1:nrow(null.samps), occs.grp.tbl)
      null.occs.ik[[1]] <- null.samps[samp.test, ]
    }
    
    # bind rows together to make full null occurrence dataset
    null.occs.i.df <- dplyr::bind_rows(null.occs.ik)
    if(eval.type == "knonspatial") {
      if(e@partition.method == "randomkfold") null.occs.i.df$grp <- get.randomkfold(null.occs.i.df, e@bg, kfolds = e@partition.settings$kfolds)$occs.grp
      if(e@partition.method == "jackknife") null.occs.i.df$grp <- get.jackknife(null.occs.i.df, e@bg)$occs.grp
    }
    null.occs.i.z <- null.occs.i.df %>% dplyr::select(-grp)
    
    # shortcuts
    categoricals <- names(which(sapply(e@occs, is.factor)))
    if(length(categoricals) == 0) categoricals <- NULL
    
    if(eval.type %in% c("testing", "none")) {
      partitions <- eval.type
      user.grp <- NULL
      user.val.grps <- NULL
    }else{
      # assign the null occurrence partitions as user partition settings, but
      # keep the empirical model background partitions
      user.grp <- list(occs.grp = null.occs.i.df$grp, bg.grp = e@bg.grp)
      # assign user validation partitions to those used in the empirical model
      user.val.grps <- cbind(e@occs, grp = e@occs.grp)
      partitions <- "user"
    }
    
    # check if ecospat is installed, and if not, prevent CBI calculations
    #if(requireNamespace("ecospat", quietly = TRUE)) {
    #  e@other.settings$ecospat.use <- TRUE
    #}else{
    #  e@other.settings$ecospat.use <- FALSE
    #}
    
    args.i <- list(occs = null.occs.i.z, bg = e@bg, tune.args = mod.settings, categoricals = categoricals, partitions = partitions,
                   algorithm = e@algorithm, other.settings = e@other.settings, partition.settings = e@partition.settings,
                   occs.testing = e@occs.testing, user.val.grps = user.val.grps, user.grp = user.grp, 
                   doClamp = e@doClamp, clamp.directions = clamp.directions.i, quiet = TRUE)
    
    null.e.i <- tryCatch({
      do.call(ENMevaluate, args.i)  
    }, error = function(cond) {
      if(quiet != TRUE) message(paste0("\n", cond, "\n"))
      # Choose a return value in case of error
      return(NULL)
    })
    
    if(is.null(null.e.i)) {
      results.na <- e@results[1,] %>% dplyr::mutate(dplyr::across(auc.train:ncoef, ~NA))
      mod.settings.i <- paste(names(mod.settings), mod.settings, collapse = "_", sep = ".")
      if(nrow(e@results.partitions) > 0) {
        results.partitions.na <- e@results.partitions %>% dplyr::filter(tune.args == mod.settings.i) %>% dplyr::mutate(dplyr::across(3:ncol(.), ~NA)) %>% dplyr::mutate(iter = i)  
      }else{
        results.partitions.na <- e@results.partitions
      }
      
      out <- list(results = results.na, results.partitions = results.partitions.na)
    }else{
      out <- list(results = null.e.i@results, 
                  results.partitions = null.e.i@results.partitions %>% dplyr::mutate(iter = i) %>% dplyr::select(iter, dplyr::everything()))  
      # restore NA row if partition evaluation is missing (model was NULL)
      if(eval.type != "testing") {
        allParts <- unique(user.grp$occs.grp) %in% out$results.partitions$fold
        if(!all(allParts)) {
          inds <- which(allParts == FALSE)
          newrow <- out$results.partitions[1,]
          newrow[,4:ncol(newrow)] <- NA
          for(ind in inds) {
            out$results.partitions <- dplyr::bind_rows(out$results.partitions, newrow %>% dplyr::mutate(fold = ind))  
          }
          out$results.partitions <- dplyr::arrange(out$results.partitions, fold)
        }  
      }
    }
    
    
    return(out)
  }
  
  #Run null models
  if(parallel == TRUE) {
    outs <- foreach::foreach(e = e.list, null.samps = null.samps.list, 
                             occs.grp.tbl = occs.grp.tbl.list, 
                             mod.settings = mod.settings.list,
                             clamp.directions.i = clamp.directions.list) %do% {
                               null <- foreach::foreach(i = 1:no.iter, .options.snow = opts, .packages = c("dplyr", "ENMeval")) %dopar% {
                                 null_i(i)}
                               return(null)
                             }
    return(outs)
  }
      
  }else{
    outs <- foreach::foreach(e = e.list, null.samps = null.samps.list, 
                             occs.grp.tbl = occs.grp.tbl.list, 
                             mod.settings = mod.settings.list,
                             clamp.directions.i = clamp.directions.list) %do% {
                               null <- foreach::foreach(i = 1:no.iter) %dopar% {
                                 null_i(i)}
                               return(null)
                               if(quiet == FALSE) setTxtProgressBar(pb, i)
                             }
    }
#START HERE
  if(quiet != TRUE) close(pb)
  if(parallel == TRUE) parallel::stopCluster(cl)
  
  #########################################
  ## 4. Run statistical tests
  #########################################
  #Extract relevant model accuracy metrics
  #Run statistical tests
  
  #########################################
  ## 5. Run post-hoc tests
  #########################################
  #Run post-hoc test (one-tailed z-test) for pairwise differences
}
