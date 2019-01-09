#' 

ENMevaluate <- function(occs, envs, bg = NULL, mod.settings, categoricals = NULL, partitions = NULL,
                        algorithm, n.bg = 10000, occs.folds = NULL, bg.folds = NULL,
                        overlap = FALSE, aggregation.factor = c(2, 2), kfolds = NA, bin.output = FALSE,
                        clamp = TRUE, rasterPreds = TRUE, parallel = FALSE, numCores = NULL, 
                        progbar = TRUE, updateProgress = FALSE, ...) {

  # record start time
  start.time <- proc.time()
  
  # general parameter checks
  if(is.null(partitions)) {
    stop("Please specify a partition method for cross validation.")
  }
  if(!(partitions %in% c("jackknife", "randomkfold", "block", "checkerboard1", "checkerboard2", "user"))) {
    stop("Please make sure partition method is one of the available options.")
  }
  if(progbar==TRUE & is.function(updateProgress)) {
    stop('Cannot specify both "progbar" and "updateProgress". Please turn one off.')
  }
  
  # mod.settings checks and algorithm-specific set-up
  if(algorithm %in% c("maxent.jar", "maxnet")) {
    if(!("rm" %in% names(mod.settings)) | !("fc" %in% names(mod.settings))) {
      stop("For Maxent, please specify both 'rm' and 'fc' settings. See ?mod.settings for help.")
    }else{
      if(!is.numeric(mod.settings[["rm"]])) {
        stop("Please input numeric values for 'rm' settings for Maxent.")
      }
      all.fc <- unlist(sapply(1:5, function(x) apply(combn(c("L","Q","H","P","T"), x), 2, function(y) paste(y, collapse = ""))))
      if(any(!mod.settings[["fc"]] %in% all.fc)) {
        stop("Please input accepted values for 'fc' settings for Maxent.")
      }
    }
    
    if(algorithm == 'maxent.jar') {
      user.args <- list(...)
      maxent.args <- c("addsamplestobackground", "addallsamplestobackground", "allowpartialdata", 
                       "beta_threshold", "beta_categorical", "beta_lqp", "beta_hinge", "convergencethreshold",
                       "defaultprevalence", "extrapolate", "fadebyclamping", "jackknife", "maximumbackground", 
                       "maximumiterations", "removeduplicates")
      if(!all(names(user.args) %in% maxent.args)) {
        stop("One or more input Maxent settings are not implemented in ENMeval or are misspelled.")
      }
      # construct user message with version info
      algorithm.ver <- paste("Maxent", maxentJARversion(), "via dismo", packageVersion('dismo'))
    }
    
    if(algorithm == 'maxnet') {
      # construct user message with version info
      algorithm.ver <- paste("maxnet", packageVersion('maxnet'))
    }
  }
  
  if(algorithm == 'brt') {
    if(all("n.trees" %in% names(mod.settings), "interaction.depth" %in% names(mod.settings), "shrinkage" %in% names(mod.settings))) {
      stop("BRT settings must include 'n.trees', 'interaction.depth', and 'shrinkage'.")
      # construct user message with version info
      algorithm.ver <- paste("gbm", packageVersion('gbm'))
    }
  }
  
  # handle user.args
  if(length(user.args) == 0) {
    user.args <- NULL
  }else{
    user.args <- paste(names(user.args), unlist(user.args), sep='=')
  }
  
  
  
  args <- make.args(mod.settings, algorithm)
  args.lab <- make.args(mod.settings, algorithm, labels = TRUE)
  
  message(paste("*** Running ENMevaluate using", algorithm.ver, "***"))

  ## data checks and formatting
  
  # if no background points specified, generate random ones
  if(is.null(bg)) {
    bg <- randomPoints(envs, n = n.bg)
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
    ptn.msg <- "checkerboard (2-fold)"
    folds <- get.checkerboard1(occs, envs, bg, aggregation.factor)
  }
  if(partitions == "checkerboard2") {
    ptn.msg <- "hierarchical checkerboard (4-fold)"
    folds <- get.checkerboard2(occs, envs, bg, aggregation.factor)
  }
  if(partitions == "user") {
    ptn.msg <- "user-defined k-fold"
    folds <- list(occs.folds = occs.folds, bg.folds = bg.folds)
  }
  
  message(paste("Doing evaluations with", ptn.msg, "cross validation..."))
  
  # run internal tuning function
  results <- tuning(occs, envs, bg, folds, algorithm,
                    args, args.lab, categoricals, aggregation.factor,
                    kfolds, bin.output, clamp, alg, rasterPreds, parallel, 
                    numCores, progbar, updateProgress, user.args)
  
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
  
  return(results)
}
