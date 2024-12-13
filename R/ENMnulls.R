#' @title Generate null ecological niche models (ENMs) and compare null with 
#' empirical performance metrics
#' @description \code{ENMnulls()} iteratively builds null ENMs for a single set 
#' of user-specified model settings based on an input ENMevaluation object, 
#' from which all other analysis settings are extracted. Summary statistics of 
#' the performance metrics for the null ENMs are taken (averages and standard 
#' deviations) and effect sizes and \emph{p}-values are calculated by comparing 
#' these summary statistics to the empirical values of the performance metrics 
#' (i.e., from the model built with the empirical data). See the references 
#' below for more details on this method.
#' 
#' @param e ENMevaluation object
#' @param mod.settings named list: one set of model settings with which to 
#' build null ENMs.
#' @param no.iter numeric: number of null model iterations.
#' @param eval.stats character vector: the performance metrics that will be 
#' used to calculate null model statistics.
#' @param user.enm ENMdetails object: if implementing a user-specified model.
#' @param user.eval function: custom function for specifying performance metrics not included in \pkg{ENMeval}.
#' The function must first be defined and then input as the argument \code{user.eval}.
#' @param user.eval.type character: if implementing a user-specified model, 
#' specify here which evaluation type to use -- either "knonspatial", 
#' "kspatial", "testing", or "none".
#' @param userStats.signs named list: user-defined evaluation statistics
#' attributed with either 1 or -1 to designate whether the expected difference 
#' between empirical and null models is positive or negative; this is used to 
#' calculate the p-value of the z-score. For example, for AUC, the difference 
#' should be positive (the empirical model should have a higher score), whereas 
#' for omission rate it should be negative (the empirical model should have a 
#' lower score).
#' @param removeMxTemp boolean: if TRUE, delete all temporary data generated 
#' when using maxent.jar for modeling.
#' @param parallel boolean: if TRUE, use parallel processing.
#' @param numCores numeric: number of cores to use for parallel processing; 
#' if NULL, all available cores will be used.
#' @param parallelType character:: either "doParallel" or "doSNOW" 
#' (default: "doSNOW").
#' @param quiet boolean: if TRUE, silence all function messages 
#' (but not errors).
#' 
#' @details This null ENM technique is based on the implementation in Bohl 
#' \emph{et al.} (2019), which follows the original methodology of Raes & ter 
#' Steege (2007) but makes an important modification: instead of evaluating 
#' each null model on random validation data, here we evaluate the null models 
#' on the same withheld validation data used to evaluate the empirical model. 
#' Bohl \emph{et al.} (2019) demonstrates this approach using a single defined 
#' withheld partition group, but Kass \emph{et al.} (2020) extended it to use 
#' spatial partitions by drawing null occurrences from the area of the predictor 
#' raster data defining each partition. Please see the vignette for a brief 
#' example.
#' 
#' This function avoids using raster data to speed up each iteration, and 
#' instead samples null occurrences from the partitioned background records. 
#' Thus, you should avoid running this when your background records are not well 
#' sampled across the study extent, as this limits the extent that null 
#' occurrences can be sampled from.
#' 
#' @references 
#' Bohl, C. L., Kass, J. M., & Anderson, R. P. (2019). A new null model approach to quantify performance and significance for ecological niche models of species distributions. \emph{Journal of Biogeography}, \bold{46}: 1101-1111. \doi{10.1111/jbi.13573}
#' 
#' Kass, J. M., Anderson, R. P., Espinosa-Lucas, A., Juárez-Jaimes, V., Martínez-Salas, E., Botello, F.,  Tavera, G., Flores-Martínez, J. J., & Sánchez-Cordero, V. (2020). Biotic predictors with phenological information improve range estimates for migrating monarch butterflies in Mexico. \emph{Ecography}, \bold{43}: 341-352. \doi{10.1111/ecog.04886}
#' 
#' Raes, N., & ter Steege, H. (2007). A null-model for significance testing of presence-only species distribution models. \emph{Ecography}, \bold{30}: 727-736. \doi{10.1111/j.2007.0906-7590.05041.x}
#' 
#' @return An \code{ENMnull} object with slots containing evaluation summary statistics for the null models 
#' and their cross-validation results, as well as differences in results between the empirical and null models. 
#' This comparison table includes z-scores of these differences and their associated p-values (under a normal distribution).
#' See ?ENMnull for more details.
#' 
#' 
#' @examples
#' \dontrun{
#' library(ENMeval)
#' 
#' # first, let's tune some models
#' occs <- read.csv(file.path(system.file(package="predicts"), 
#' "/ex/bradypus.csv"))[,2:3]
#' envs <- rast(list.files(path=paste(system.file(package="predicts"), 
#'                                    "/ex", sep=""), pattern="tif$", 
#'                                    full.names=TRUE))
#' bg <- as.data.frame(predicts::backgroundSample(envs, n = 10000))
#' names(bg) <- names(occs)
#' 
#' ps <- list(orientation = "lat_lat")
#' 
#' # as an example, let's use two user-specified evaluation metrics
#' conf.and.cons <- function(vars) {
#'   observations <- c(
#'     rep(x = 1, times = length(vars$occs.train.pred)),
#'     rep(x = 0, times = length(vars$bg.train.pred)),
#'     rep(x = 1, times = length(vars$occs.val.pred)),
#'     rep(x = 0, times = length(vars$bg.val.pred))
#'   )
#'   predictions <- c(vars$occs.train.pred, vars$bg.train.pred, 
#'   vars$occs.val.pred, vars$bg.val.pred)
#'   evaluation_mask <- c(
#'     rep(x = FALSE, times = length(vars$occs.train.pred) + 
#'     length(vars$bg.train.pred)),
#'     rep(x = TRUE, times = length(vars$occs.val.pred) + 
#'     length(vars$bg.val.pred))
#'   )
#'   measures <- confcons::measures(observations = observations, 
#'   predictions = predictions, 
#'   evaluation_mask = evaluation_mask, df = TRUE)[, c("CPP_eval", "DCPP")]
#'   colnames(measures) <- c("confidence", "consistency")
#'   return(measures)
#' }
#' 
#' e <- ENMevaluate(occs, envs, bg, 
#'                  tune.args = list(fc = c("L","LQ","LQH"), rm = 2:4), 
#'                  partitions = "block", partition.settings = ps, 
#'                  algorithm = "maxnet", categoricals = "biome",
#'                  user.eval  = conf.and.cons, parallel = TRUE)
#' 
#' d <- eval.results(e)
#' 
#' # here, we will choose an optimal model based on validation CBI, but you can
#' # choose yourself what evaluation statistics to use
#' opt <- d |> filter(cbi.val.avg == max(cbi.val.avg))
#' 
#' # now we can run our null models, and we can specify to include estimates for
#' # our user-specified variables too, but we need to make sure we note what 
#' # sign we expect these statistics to be 
#' # NOTE: you should use at least 100 iterations in practice -- this is just an
#' # example
#' nulls <- ENMnulls(e, 
#'                   mod.settings = list(fc = opt$fc, rm = opt$rm),
#'                   no.iter = 10, 
#'                   user.eval = conf.and.cons,
#'                   eval.stats = c("cbi.val", "confidence", "consistency"),
#'                   userStats.signs = c("confidence" = 1, "consistency" = 1))
#' 
#' # here are the results of all the null iterations
#' null.results(nulls)
#' # and here are the comparisons between the null and empirical values for
#' # the evaluation statistics, including the z-score and p-value
#' # for more details, see Bohl et al. 2019
#' null.emp.results(nulls)
#' }
#' 
#' @export
#'

# for split evaluation, label training occs "1" and testing evaluation occs "2" in partitions
ENMnulls <- function(e, mod.settings, no.iter, eval.stats = c("auc.val","auc.diff","cbi.val","or.mtp","or.10p"),
                     user.enm = NULL, user.eval = NULL, user.eval.type = NULL, userStats.signs = NULL, 
                     removeMxTemp = TRUE, parallel = FALSE, numCores = NULL, parallelType = "doSNOW", quiet = FALSE) {
  
  # assign evaluation type based on partition method
  if(is.null(user.eval.type)) {
    eval.type <- switch(e@partition.method,
                        randomkfold = "knonspatial",
                        jackknife = "knonspatial",
                        block = "kspatial",
                        checkerboard = "kspatial",
                        testing = "testing",
                        none = "none")  
  }else{
    eval.type <- user.eval.type
  }
  
  ## checks
  # model settings are all single entries
  if(!all(sapply(mod.settings, length) == 1)) stop("Please input a single set of model settings.")
  # model settings are correct for input algorithm and are entered in the right order -- 
  #if not, put them in the right order, else indexing models later will fail because the model 
  # name will be incorrect
  if(e@algorithm %in% c("maxent.jar", "maxnet")) {
    if(length(mod.settings) != 2) {
      stop("Please input two complexity settings (fc [feature classes] and rm [regularization
           multipliers]) for mod.settings for maxent.jar and maxnet models.")
    }
    if(all(names(mod.settings) %in% c("fc", "rm"))) {
      if(!all(names(mod.settings) == c("fc", "rm"))) {
        mod.settings <- mod.settings[c("fc", "rm")]
      }
    }else{
      stop('Please input only "fc" (feature classes) and "rm" (regularization multipliers) for
           mod.settings for maxent.jar and maxnet models.')
    }
  }else if(e@algorithm == "bioclim") {
    if(length(mod.settings) != 1) {
      stop("Please input one complexity setting (tails) for mod.settings for BIOCLIM models.")
    }
    if(!all(names(mod.settings) == "tails")) {
      stop('Please input only "tails" for mod.settings for BIOCLIM models.')
    }
  }
  
  
  # assign directionality of sign for evaluation stats
  signs <- c(list("auc.val" = 1, "auc.train" = 1, "cbi.val" = 1, "cbi.train" = 1,
                  "auc.diff" = -1, "or.10p" = -1, "or.mtp" = -1), userStats.signs)
  
  # record start time
  start.time <- proc.time()
  
  # assign the number of cross validation iterations
  nk <- max(as.numeric(as.character(e@occs.grp)))
  
  # get number of occurrence points by partition
  occs.grp.tbl <- table(e@occs.grp)
  
  # if more than one background partition exists, assume spatial CV and
  # keep existing partitions
  null.samps <- cbind(rbind(e@occs, e@bg), grp = c(e@occs.grp, e@bg.grp))
  
  if(e@algorithm == "maxent.jar") {
    # create temp directory to store maxent.jar output, for potential removal later
    tmpdir <- paste(tempdir(), runif(1,0,1), sep = "/")
    dir.create(tmpdir, showWarnings = TRUE, recursive = FALSE)
  }
  
  # assign user algorithm if provided
  if(!is.null(user.enm)) {
    e@algorithm <- user.enm
  }
  
  ############################## #
  # specify empirical model statistics ####
  ############################## #
  
  mod.tune.args  <- paste(names(mod.settings), mod.settings, collapse = "_", sep = ".")
  emp.mod <- e@models[[mod.tune.args]]
  emp.mod.res <- e@results |> dplyr::filter(tune.args == mod.tune.args)
  
  ############################## #
  # build null models ####
  ############################## #
  
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
  
  if(length(e@clamp.directions) == 0) clamp.directions.i <- NULL else clamp.directions.i <- e@clamp.directions
  
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
        null.samps.k <- null.samps |> dplyr::filter(grp == k)
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
    null.occs.i.z <- null.occs.i.df |> dplyr::select(-grp)
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
    if((requireNamespace("ecospat", quietly = TRUE))) {
      e@other.settings$ecospat.use <- TRUE
    }else{
      e@other.settings$ecospat.use <- FALSE
    }
    
    args.i <- list(occs = null.occs.i.z, bg = e@bg, tune.args = mod.settings, categoricals = categoricals, partitions = partitions,
                   algorithm = e@algorithm, other.settings = e@other.settings, partition.settings = e@partition.settings,
                   occs.testing = e@occs.testing, user.val.grps = user.val.grps, user.grp = user.grp, user.eval = user.eval,
                   doClamp = e@doClamp, clamp.directions = clamp.directions.i, quiet = TRUE)
    
    null.e.i <- tryCatch({
      do.call(ENMevaluate, args.i)  
    }, error = function(cond) {
      if(quiet != TRUE) message(paste0("\n", cond, "\n"))
      # Choose a return value in case of error
      return(NULL)
    })
    
    if(is.null(null.e.i)) {
      results.na <- emp.mod.res |> dplyr::mutate(dplyr::across(auc.train:ncoef, ~NA))
      mod.settings.i <- paste(names(mod.settings), mod.settings, collapse = "_", sep = ".")
      if(nrow(e@results.partitions) > 0) {
        nc <- ncol(e@results.partitions)
        results.partitions.na <- e@results.partitions |> dplyr::filter(tune.args == emp.mod.res$tune.args) |> dplyr::mutate(dplyr::across(3:nc, ~NA)) |> dplyr::mutate(iter = i)  
      }else{
        results.partitions.na <- e@results.partitions
      }
      
      out <- list(results = results.na, results.partitions = results.partitions.na)
    }else{
      out <- list(results = null.e.i@results, 
                  results.partitions = null.e.i@results.partitions |> dplyr::mutate(iter = i) |> dplyr::select(iter, dplyr::everything()))  
      # restore NA row if partition evaluation is missing (model was NULL)
      if(eval.type != "testing") {
        allParts <- unique(user.grp$occs.grp) %in% out$results.partitions$fold
        if(!all(allParts)) {
          inds <- which(allParts == FALSE)
          newrow <- out$results.partitions[1,]
          newrow[,4:ncol(newrow)] <- NA
          for(ind in inds) {
            out$results.partitions <- dplyr::bind_rows(out$results.partitions, newrow |> dplyr::mutate(fold = ind))  
          }
          out$results.partitions <- dplyr::arrange(out$results.partitions, fold)
        }  
      }
    }
    
    return(out)
  }
  
  if(parallel == TRUE) {
    outs <- foreach::foreach(i = 1:no.iter, .options.snow = opts, .export = c("e"), .packages = c("dplyr", "ENMeval")) %dopar% {
      null_i(i)
    }
  }else{
    outs <- list()
    for(i in 1:no.iter) {
      outs[[i]] <- null_i(i)
      if(quiet == FALSE) setTxtProgressBar(pb, i)
    }  
  }
  
  if(quiet != TRUE) close(pb)
  if(parallel == TRUE) parallel::stopCluster(cl)
  
  # assemble null evaluation statistics and take summaries
  nulls.ls <- lapply(outs, function(x) x$results)
  nulls.grp.ls <- lapply(outs, function(x) x$results.partitions)
  nulls <- dplyr::bind_rows(nulls.ls) |> dplyr::select(!dplyr::contains("AIC"))
  nulls.grp <- dplyr::bind_rows(nulls.grp.ls)
  if(eval.type %in% c("testing", "none")) {
    nulls.avgs <- nulls |> dplyr::select(dplyr::ends_with("train"), dplyr::contains(eval.stats)) |> dplyr::summarize_all(mean, na.rm = TRUE)
    nulls.sds <- nulls |> dplyr::select(dplyr::ends_with("train"), dplyr::contains(eval.stats)) |> dplyr::summarise_all(sd, na.rm = TRUE)
    # get empirical model evaluation statistics for comparison
    emp.avgs <- emp.mod.res |> dplyr::select(dplyr::ends_with("train"), dplyr::contains(eval.stats))
  }else{
    nulls.avgs <- nulls |> dplyr::select(dplyr::ends_with("train"), paste0(eval.stats, ".avg")) |> dplyr::summarize_all(mean, na.rm = TRUE)
    nulls.sds <- nulls |> dplyr::select(dplyr::ends_with("train"), paste0(eval.stats, ".sd")) |> dplyr::summarise_all(sd, na.rm = TRUE)
    # get empirical model evaluation statistics for comparison
    emp.avgs <- emp.mod.res |> dplyr::select(dplyr::ends_with("train"), paste0(eval.stats, ".avg"))
  }
  if(sum(grepl("sd", names(emp.mod.res))) > 0) {
    emp.sds <- emp.mod.res |> dplyr::select(paste0(eval.stats, ".sd"))  
  }else{
    emp.sds <- NULL
  }
  
  empNull.stats <- as.data.frame(matrix(nrow = 6, ncol = ncol(emp.avgs)+1))
  names(empNull.stats)[1] <- "statistic"
  empNull.stats$statistic <- c("emp.mean", "emp.sd", "null.mean", "null.sd", "zscore", "pvalue")
  names(empNull.stats)[-1] <- gsub(".avg", replacement = "", names(emp.avgs))
  
  # populate empirical and null means and standard deviations
  empNull.stats[1, -1] <- emp.avgs
  emp.sds.nameStrip <- gsub(".sd", "", names(emp.sds))
  if(length(emp.sds.nameStrip) > 0) empNull.stats[2, which(names(empNull.stats) %in% emp.sds.nameStrip)] <- emp.sds
  empNull.stats[3,-1] <- nulls.avgs
  empNull.stats[4,-1] <- nulls.sds
  # calculate z-scores
  empNull.stats[5,-1] <- (emp.avgs - nulls.avgs) / nulls.sds
  # find statistics that need a positive pnorm, and those that need a negative pnorm
  p.pos <- names(signs[sapply(signs, function(x) x == 1)])
  p.neg <- names(signs[sapply(signs, function(x) x == -1)])
  p.pos.inds <- which(grepl(paste(p.pos, collapse = "|"), names(empNull.stats)))
  p.neg.inds <- which(grepl(paste(p.neg, collapse = "|"), names(empNull.stats)))
  # calculate p-values
  empNull.stats[6, p.pos.inds] <- sapply(empNull.stats[5, p.pos.inds], function(x) pnorm(x, lower.tail = FALSE))
  empNull.stats[6, p.neg.inds] <- sapply(empNull.stats[5, p.neg.inds], pnorm)
  
  mod.settings.tbl <- expand.grid(mod.settings)
  mod.settings.tbl$tune.args <- apply(mod.settings.tbl, 1, function(x) paste(names(x), x, collapse = "_", sep = "."))
  mod.settings.tbl <- dplyr::mutate_all(mod.settings.tbl, as.factor)
  
  # make cbi.val column NA for jackknife partitions
  if(e@partition.method == "jackknife") {
    empNull.stats$cbi.val <- NA
    nulls$cbi.val.avg <- NA
    nulls$cbi.val.sd <- NA
    nulls.grp$cbi.val <- NA
  }
  
  # first, add clamp directions to other settings so the object can record it
  e@other.settings$clamp.directions <- e@clamp.directions
  # condense mod.args to named matrix for inserting into class slot
  e.n <- ENMnull(null.algorithm = e@algorithm,
                 null.mod.settings = mod.settings.tbl,
                 null.partition.method = e@partition.method,
                 null.partition.settings = e@partition.settings,
                 null.doClamp = e@doClamp,
                 null.other.settings = e@other.settings,
                 null.no.iter = no.iter,
                 null.results = nulls,
                 null.results.partitions = nulls.grp,
                 null.emp.results = empNull.stats,
                 emp.occs = e@occs,
                 emp.occs.grp = e@occs.grp,
                 emp.bg = e@bg,
                 emp.bg.grp = e@bg.grp)
  
  # optionally remove temp directory for maxent.jar
  if(e@algorithm == "maxent.jar" & removeMxTemp == TRUE) unlink(tmpdir, recursive = TRUE, force = TRUE)
  
  
  timed <- proc.time() - start.time
  t.min <- floor(timed[3] / 60)
  t.sec <- timed[3] - (t.min * 60)
  if(quiet == FALSE) message(paste("\nENMnulls completed in", t.min, "minutes", round(t.sec, 1), "seconds."))
  
  return(e.n)
}
