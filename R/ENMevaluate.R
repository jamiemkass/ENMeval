#' 

ENMevaluate <- function(occs, envs, bg = NULL, mod.settings, categoricals = NULL, partitions = NULL,
                        algorithm, n.bg = 10000, occs.folds = NULL, bg.folds = NULL,
                        overlap = FALSE, aggregation.factor = c(2, 2), kfolds = NA, bin.output = FALSE,
                        clamp = TRUE, skipRasters = FALSE, parallel = FALSE, numCores = NULL, 
                        updateProgress = FALSE, ...) {

  # record start time
  start.time <- proc.time()
  
  # general parameter checks
  if(is.null(partitions)) {
    stop("Please specify a partition method for cross validation.")
  }
  if(!(partitions %in% c("jackknife", "randomkfold", "block", "checkerboard1", "checkerboard2", "user"))) {
    stop("Please make sure partition method is one of the available options.")
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
  
  # unpack occs.folds and bg.folds
  occs.folds <- folds$occs.folds
  bg.folds <- folds$bg.folds
  
  message(paste("Doing evaluations with", ptn.msg, "cross validation..."))
  
  ###############
  # ANALYSIS ####
  ###############
  
  # extract predictor variable values at coordinates for occs and bg
  occs.vals <- as.data.frame(extract(envs, occs))
  bg.vals <- as.data.frame(extract(envs, bg))
  
  # remove rows from occs, occs.vals, occs.folds with NA for any predictor variable
  occs.vals.na <- sum(rowSums(is.na(occs.vals)) > 0)
  bg.vals.na <- sum(rowSums(is.na(bg.vals)) > 0)
  if(occs.vals.na > 0) {
    i <- !apply(occs.vals, 1, anyNA)
    occs <- occs[i,]
    occs.folds <- occs.folds[i]
    occs.vals <- occs.vals[i,]
    message(paste("There were", occs.vals.na, "occurrence records with NA for at least
                  one predictor variable. Removed these from analysis,
                  resulting in", nrow(occs.vals), "occurrence records."))
  }
  # do the same for associated bg variables
  if(bg.vals.na > 0) {
    i <- !apply(bg.vals, 1, anyNA)
    occs <- occs[i,]
    bg.folds <- bg.folds[i]
    bg.vals <- bg.vals[i,]
    message(paste("There were", bg.vals.na, "background records with NA for at least
                  one predictor variable. Removed these from analysis,
                  resulting in", nrow(bg.vals), "background records."))
  }
  
  # convert fields for categorical data to factor class
  if(!is.null(categoricals)) {
    for (i in 1:length(categoricals)) {
      occs.vals[, categoricals[i]] <- as.factor(occs.vals[, categoricals[i]])
      bg.vals[, categoricals[i]] <- as.factor(bg.vals[, categoricals[i]])
    }
  }
  
  # if using parallel processing...
  if(parallel == TRUE) {
    # set up parallel processing functionality
    allCores <- detectCores()
    if (is.null(numCores)) {
      numCores <- allCores
    }
    c1 <- makeCluster(numCores)
    registerDoParallel(c1)
    numCoresUsed <- getDoParWorkers()
    message(paste("Of", allCores, "total cores using", numCoresUsed))
    
    message("Running in parallel...")
    if (algorithm == 'maxnet') {
      tune.results <- foreach(i = seq_len(length(args)),
                     .packages = c("dismo", "raster", "ENMeval", "maxnet")) %dopar% {
                       modelTune.maxnet(occs.vals, bg.vals, envs, folds, args[[i]],
                                        rasterPreds, clamp)
                     }
    }
    if (algorithm == 'maxent.jar') {
      tune.results <- foreach(i = seq_len(length(args)),
                     .packages = c("dismo", "raster", "ENMeval", "rJava")) %dopar% {
                       modelTune.maxentJar(occs.vals, bg.vals, envs, folds, args[[i]],
                                           userArgs, rasterPreds, clamp)
                     }
    }
    stopCluster(c1)
  } 
  # if not using parallel processing...
  if(parallel == FALSE) {
    tune.results <- list()
    # set up the console progress bar
    pb <- txtProgressBar(0, length(args), style = 3)
    
    for(i in 1:length(args)) {
      # and (optionally) the shiny progress bar (updateProgress)
      if(length(args) > 1) {
        if(is.function(updateProgress)) {
          text <- paste0('Running ', args.lab[[1]][i], args.lab[[2]][i], '...')
          updateProgress(detail = text)
        }
        setTxtProgressBar(pb, i)
      }
      if (algorithm == 'maxnet') {
        tune.results[[i]] <- modelTune.maxnet(occs.vals, bg.vals, envs, occs.folds, bg.folds, args[[i]],
                                     rasterPreds, clamp)
      } else if (algorithm == 'maxent.jar') {
        tune.results[[i]] <- modelTune.maxentJar(occs.vals, bg.vals, envs, nk, folds, args[[i]],
                                        userArgs, rasterPreds, clamp)
      }
    }
    close(pb)  
  }
  
  # gather all full models into list
  mod.full.all <- lapply(tune.results, function(x) x$mod.full)
  # gather all statistics into a data frame
  stats.all <- lapply(tune.results, function(x) x$stats)
  if(skipRasters == FALSE) {
    mod.full.pred.all <- stack(sapply(tune.results, function(x) x$mod.full.pred))
  } else {
    mod.full.pred.all <- stack()
  }
  
  # define a corrected variance function
  corrected.var <- function(x, nk){
    rowSums((x - rowMeans(x))^2) * ((nk-1)/nk)
  }
  
  # make data frame of stats (for all folds)
  stats.all.df <- do.call("rbind", stats.all)
  # make new column for fold number
  stats.all.df$fold <- rep(1:4, length(args))
  # concatenate fc and rm to make settings column
  stats.all.df$fc <- unlist(lapply(args.lab[[1]], function(x) rep(x, 4)))
  stats.all.df$rm <- unlist(lapply(args.lab[[2]], function(x) rep(x, 4)))
  args.lab.concat <- paste0(args.lab[[1]], args.lab[[2]])
  stats.all.df$settings <- unlist(lapply(args.lab.concat, function(x) rep(x, 4)))
  stats.all.df <- stats.all.df %>% dplyr::select(settings, fc, rm, fold, auc.test, auc.diff, or.min, or.10)

  stats.sum.df <- stats.all.df %>% dplyr::group_by(settings) %>% summarize(mean(auc.test), mean(auc.diff), mean(or.min), mean(or.10))
  
  
  
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
