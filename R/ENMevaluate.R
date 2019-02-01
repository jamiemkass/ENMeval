#' @export 

ENMevaluate <- function(occs, envs, bg = NULL, mod.fun, tune.args, other.args, categoricals = NULL, 
                        partitions = NULL, occs.folds = NULL, bg.folds = NULL, n.bg = 10000, 
                        overlap = FALSE, aggregation.factor = c(2, 2), kfolds = NA, bin.output = FALSE,
                        doClamp = TRUE, skipRasters = FALSE, abs.auc.diff = TRUE, parallel = FALSE, 
                        numCores = NULL, updateProgress = FALSE) {

  # record start time
  start.time <- proc.time()
  # get model function's name
  mod.fun.name <- as.character(substitute(mod.fun))[3]
  
  ########### #
  # CHECKS ####
  ########### #
  
  # general parameter checks
  if(is.null(partitions)) {
    stop("Please specify a partition method for cross validation.")
  }
  if(!(partitions %in% c("jackknife", "randomkfold", "block", "checkerboard1", "checkerboard2", "user"))) {
    stop("Please make sure partition method is one of the available options.")
  }
  
  # tune.args checks and algorithm-specific set-up
  if(mod.fun.name %in% c("maxent", "maxnet")) {
    if(!("rm" %in% names(tune.args)) | !("fc" %in% names(tune.args))) {
      stop("For Maxent, please specify both 'rm' and 'fc' settings. See ?tune.args for help.")
    }else{
      if(!is.numeric(tune.args[["rm"]])) {
        stop("Please input numeric values for 'rm' settings for Maxent.")
      }
      all.fc <- unlist(sapply(1:5, function(x) apply(combn(c("L","Q","H","P","T"), x), 2, function(y) paste(y, collapse = ""))))
      if(any(!tune.args[["fc"]] %in% all.fc)) {
        stop("Please input accepted values for 'fc' settings for Maxent.")
      }
    }
    
    if(mod.fun.name == 'maxent') {
      # maxent.args <- c("fc", "rm", "addsamplestobackground", "addallsamplestobackground", "allowpartialdata", 
      #                  "beta_threshold", "beta_categorical", "beta_lqp", "beta_hinge", "convergencethreshold",
      #                  "defaultprevalence", "extrapolate", "fadebyclamping", "jackknife", "maximumbackground", 
      #                  "maximumiterations", "removeduplicates")
      # if(!all(names(tune.args) %in% maxent.args)) {
        # stop("One or more input Maxent settings are not implemented in ENMeval or are misspelled.")
      # }
      # construct user message with version info
      algorithm.ver <- paste("Maxent", maxentJARversion(), "via dismo", packageVersion('dismo'))
    }
    
    if(mod.fun.name == 'maxnet') {
      # construct user message with version info
      algorithm.ver <- paste("maxnet", packageVersion('maxnet'))
    }
  }
  
  if(mod.fun.name == 'brt') {
    if(all("n.trees" %in% names(tune.args), "interaction.depth" %in% names(tune.args), "shrinkage" %in% names(tune.args))) {
      stop("BRT settings must include 'n.trees', 'interaction.depth', and 'shrinkage'.")
      # construct user message with version info
      algorithm.ver <- paste("gbm", packageVersion('gbm'))
    }
  }
  
  tune.tbl <- expand.grid(tune.args, stringsAsFactors = FALSE)
  
  message(paste("*** Running ENMevaluate using", algorithm.ver, "***"))

  ## data checks and formatting
  
  # if no background points specified, generate random ones
  if(is.null(bg)) {
    bg <- dismo::randomPoints(envs, n = n.bg)
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
  
  ############# #
  # ANALYSIS ####
  ############# #
  
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
    for(i in 1:length(categoricals)) {
      occs.vals[, categoricals[i]] <- as.factor(occs.vals[, categoricals[i]])
      bg.vals[, categoricals[i]] <- as.factor(bg.vals[, categoricals[i]])
    }
  }
  
  if(parallel == TRUE) {
    results <- tune.enm.parallel(occs.vals, bg.vals, occs.folds, bg.folds, envs, 
                                 mod.fun, tune.tbl, other.args, categoricals, 
                                 doClamp, skipRasters, abs.auc.diff)
  } 
  else {
    results <- tune.enm(occs.vals, bg.vals, occs.folds, bg.folds, envs, 
                        mod.fun, tune.tbl, other.args, categoricals, 
                        doClamp, skipRasters, abs.auc.diff, updateProgress)
  }
  
  # gather all full models into list
  mod.full.all <- lapply(results, function(x) x$mod.full)
  # gather all statistics into a data frame
  kstats.all <- lapply(results, function(x) x$kstats)
  if(skipRasters == FALSE) {
    mod.full.pred.all <- raster::stack(sapply(results, function(x) x$mod.full.pred))
  } else {
    mod.full.pred.all <- raster::stack()
  }
  
  # define a corrected variance function
  corrected.var <- function(x, nk){
    rowSums((x - rowMeans(x))^2) * ((nk-1)/nk)
  }
  
  # make data frame of stats (for all folds)
  kstats.df <- do.call("rbind", kstats.all)
  # define number of folds (the value of "k")
  nk <- length(unique(occs.folds))
  # concatenate fc and rm to make settings column
  for(i in 1:ncol(tune.tbl)) {
    kstats.df <- cbind(kstats.df, rep(tune.tbl[,i], each = nk))
    names(kstats.df)[ncol(kstats.df)] <- names(tune.tbl)[i]
  }
  # make new column for fold number
  kstats.df$fold <- rep(1:nk, nrow(tune.tbl))
  # get number of columns in kstats.df
  nc <- ncol(kstats.df)
  # rearrange the columns
  kstats.df <- kstats.df[, c(seq(nc-ncol(tune.tbl), nc), seq_len(nc-ncol(tune.tbl)-1))]
  
  stats.sum.df <- kstats.df %>% dplyr::group_by(fc, rm) %>% dplyr::summarize(mean(auc.test), mean(auc.diff), mean(or.min), mean(or.10))
  
  
  
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
