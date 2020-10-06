#' @title Generate null ecological niche models (ENMs) and compare evaluation statistics with real ENMs
#' @description \code{ENMnullSims()} builds null ecological niche models iteratively for a single set of 
#' user-specified model settings based on an input ENMevaluation object, from which all other analysis 
#' settings are copied. \code{ENMnullSims()} returns an \code{ENMnull} object with slots containing evaluation
#' summary statistics for the null models and their cross-validation results, as well as differences
#' in results between the real and null models. This comparison table includes z-scores of these differences
#' and their associated p-values (under a normal distribution).

#' 
#' @param e ENMevaluation object
#' @param mod.settings named list of one set of model settings to build null models with
#' @param no.iter numeric of number of null model iterations
#' @param envs Raster* object of environmental variables (must be in same geographic projection as occurrence data);
#' necessary when evaluating using spatial cross-validation (see details)
#' @param user.enm ENMdetails object specified by the user
#' @param userStats.signs named list of user-defined evaluation statistics attributed with
#' either 1 or -1 to designate whether the expected difference between real and null models is 
#' positive or negative; this is used to calculate the p-value of the z-score
#' @param removeMxTemp boolean (TRUE or FALSE) which if TRUE delete all temporary data generated
#' when using maxent.jar for modeling
#' @export
#'

# for split evaluation, label training occs "1" and testing evaluation occs "2" in partitions
ENMnullSims <- function(e, mod.settings, no.iter, user.enm = NULL, 
                     userStats.signs = NULL, removeMxTemp = TRUE, quiet = FALSE) {

  # assign evaluation type based on partition method
  eval.type <- switch(e@partition.method,
                      randomkfold = "knonspatial",
                      jackknife = "knonspatial",
                      block = "kspatial",
                      checkerboard1 = "kspatial",
                      checkerboard2 = "kspatial",
                      testing = "knonspatial")
  
  # checks
  if(!all(sapply(mod.settings, length) == 1)) stop("Please input a single set of model settings.")
  
  # assign directionality of sign for evaluation stats
  signs <- c(list("auc.val" = 1, "auc.train" = 1, "cbi.val" = 1, "cbi.train" = 1,
                   "auc.diff" = -1, "or.10p" = -1, "or.mtp" = -1), userStats.signs)
  
  # record start time
  start.time <- proc.time()

  # assign the number of cross validation iterations
  nk <- ifelse(e@partition.method == "testing", 1, max(as.numeric(e@occs.grp)))

  # get number of occurrence points by partition
  occs.grp.tbl <- table(e@occs.grp)
  
  # if more than one background partition exists, assume spatial CV and
  # keep existing partitions
  null.samps <- cbind(rbind(e@occs, e@bg), grp = c(e@occs.grp, e@bg.grp))
  
  # assign algorithm
  if(is.null(user.enm)) {
    enm <- lookup.enm(e@algorithm)
    if(e@algorithm == "maxent.jar") {
      # create temp directory to store maxent.jar output, for potential removal later
      tmpdir <- paste(tempdir(), runif(1,0,1), sep = "/")
      dir.create(tmpdir, showWarnings = TRUE, recursive = FALSE)
    }
  }else{
    enm <- user.enm
  }

  ############################## #
  # specify real model statistics ####
  ############################## #

  mod.tune.args  <- paste(names(mod.settings), mod.settings, collapse = "_", sep = ".")
  real.mod <- e@models[[mod.tune.args]]
  real.mod.res <- e@results %>% dplyr::filter(tune.args == mod.tune.args)

  ############################## #
  # build null models ####
  ############################## #

  # initialize list to record stats for null iteration i
  nulls.ls <- list()
  nulls.grp.ls <- list()

  if(quiet == FALSE) message(paste("Building and evaluating null ENMs with", no.iter, "iterations..."))
  if(quiet == FALSE) message("Sampling null occurrences from background values...")
  if(quiet == FALSE) pb <- txtProgressBar(0, no.iter, style = 3)

  for(i in 1:no.iter) {

    null.occs.ik <- list()
    if(eval.type == "kspatial") {
      # randomly sample the same number of training occs over each k partition
      # of envs; if kspatial evaluation, only sample over the current spatial
      # partition of envs.z
      for(k in 1:nk) {
        # sample null occurrences only from
        # the environment variable grid cells in partition group k
        null.samps.k <- null.samps %>% dplyr::filter(grp == k)
        # randomly sample n null occurrences, where n equals the number
        # of real occurrence in partition group k
        samp.k <- sample(1:nrow(null.samps.k), occs.grp.tbl[k])
        null.occs.ik[[k]] <- null.samps.k[samp.k, ]
      }
    }else if(eval.type == "knonspatial") {
      for(k in 1:nk) {
        # randomly sample n null occurrences, where n equals the number
        # of real occurrence in partition group k
        samp.k <- sample(1:nrow(null.samps), occs.grp.tbl[k])
        null.occs.ik[[k]] <- null.samps[samp.k, ]
      }
    }

    # bind rows together to make full null occurrence dataset
    null.occs.i.df <- dplyr::bind_rows(null.occs.ik)
    if(eval.type == "knonspatial") {
      if(e@partition.method == "randomkfold") null.occs.i.df$grp <- get.randomkfold(null.occs.i.df, e@bg, kfolds = e@partition.settings$kfolds)$occs.grp
      if(e@partition.method == "jackknife") null.occs.i.df$grp <- get.jackknife(null.occs.i.df, e@bg)$occs.grp
    }
    null.occs.i.z <- null.occs.i.df %>% dplyr::select(-grp)
    # assign the null occurrence partitions as user partition settings, but
    # keep the real model background partitions
    user.grp <- list(occs.grp = null.occs.i.df$grp, bg.grp = e@bg.grp)
    # shortcuts for settings
    e.s <- e@other.settings
    # assign user validation partitions to those used in the real model
    user.val.grps <- cbind(e@occs, grp = e@occs.grp)
    e.p <- e@partition.settings
    categoricals <- names(which(sapply(e@occs, is.factor)))
    if(length(categoricals) == 0) categoricals <- NULL

    null.e.i <- ENMevaluate(occs = null.occs.i.z, bg = e@bg, tune.args = mod.settings, categoricals = categoricals,
                            algorithm = e@algorithm, other.args = e.s$other.args, partitions = "user",
                            user.val.grps = user.val.grps, user.grp = user.grp, kfolds = e.p$kfolds, 
                            aggregation.factor = e.p$aggregation.factor, clamp = e.s$clamp, 
                            pred.type = e.s$pred.type, abs.auc.diff = e.s$abs.auc.diff, quiet = TRUE)
    if(quiet == FALSE) setTxtProgressBar(pb, i)

    nulls.ls[[i]] <- null.e.i@results
    nulls.grp.ls[[i]] <- null.e.i@results.partitions %>% dplyr::mutate(iter = i) %>% dplyr::select(iter, dplyr::everything())
    
    # restore NA row if partition evaluation is missing (model was NULL)
    allParts <- unique(user.grp$occs.grp) %in% nulls.grp.ls[[i]]$fold
    if(!all(allParts)) {
      inds <- which(allParts == FALSE)
      newrow <- nulls.grp.ls[[i]][1,]
      newrow[,4:ncol(newrow)] <- NA
      for(ind in inds) {
        nulls.grp.ls[[i]] <- bind_rows(nulls.grp.ls[[i]], newrow %>% mutate(fold = ind))  
      }
      nulls.grp.ls[[i]] <- arrange(nulls.grp.ls[[i]], fold)
    }
  }

  # assemble null evaluation statistics and take summaries
  nulls <- dplyr::bind_rows(nulls.ls)
  nulls.grp <- dplyr::bind_rows(nulls.grp.ls)
  nulls.avgs <- nulls %>% dplyr::select(dplyr::ends_with("train"), dplyr::ends_with("avg")) %>% dplyr::summarize_all(mean)
  nulls.sds <- nulls %>% dplyr::select(dplyr::ends_with("train"), dplyr::ends_with("avg")) %>% dplyr::summarise_all(sd)
  # get real model evaluation statistics for comparison
  real.avgs <- real.mod.res %>% dplyr::select(dplyr::ends_with("train"), dplyr::ends_with("avg"))
  real.sds <- real.mod.res %>% dplyr::select(dplyr::ends_with("sd"))
  
  realNull.stats <- as.data.frame(matrix(nrow = 6, ncol = ncol(real.avgs)+1))
  names(realNull.stats)[1] <- "statistic"
  realNull.stats$statistic <- c("real.mean", "real.sd", "null.mean", "null.sd", "zscore", "pvalue")
  names(realNull.stats)[-1] <- gsub(".avg", replacement = "", names(real.avgs))
  
  # populate real and null means and standard deviations
  realNull.stats[1, -1] <- real.avgs
  real.sds.nameStrip <- gsub(".sd", "", names(real.sds))
  realNull.stats[2, which(names(realNull.stats) %in% real.sds.nameStrip)] <- real.sds
  realNull.stats[3,-1] <- nulls.avgs
  realNull.stats[4,-1] <- nulls.sds
  # calculate z-scores
  realNull.stats[5,-1] <- (real.avgs - nulls.avgs) / nulls.sds
  # find statistics that need a positive pnorm, and those that need a negative pnorm
  p.pos <- names(signs[sapply(signs, function(x) x == 1)])
  p.neg <- names(signs[sapply(signs, function(x) x == -1)])
  p.pos.inds <- which(grepl(paste(p.pos, collapse = "|"), names(realNull.stats)))
  p.neg.inds <- which(grepl(paste(p.neg, collapse = "|"), names(realNull.stats)))
  # calculate p-values
  realNull.stats[6, p.pos.inds] <- sapply(realNull.stats[5, p.pos.inds], function(x) pnorm(x, lower.tail = FALSE))
  realNull.stats[6, p.neg.inds] <- sapply(realNull.stats[5, p.neg.inds], pnorm)
  
  mod.settings.tbl <- expand.grid(mod.settings)
  mod.settings.tbl$tune.args <- apply(mod.settings.tbl, 1, function(x) paste(names(x), x, collapse = "_", sep = "."))
  mod.settings.tbl <- dplyr::mutate_all(mod.settings.tbl, as.factor)
  
  # make cbi.val column NA for jackknife partitions
  if(e@partition.method == "jackknife") {
    nulls.grp$cbi.val <- NA
  }
  
  # condense mod.args to named matrix for inserting into class slot
  e.n <- ENMnull(null.algorithm = e@algorithm,
                 null.mod.settings = mod.settings.tbl,
                 null.partition.method = e@partition.method,
                 null.partition.settings = e@partition.settings,
                 null.other.settings = e@other.settings,
                 no.iter = no.iter,
                 null.results = nulls,
                 null.results.partitions = nulls.grp,
                 real.vs.null.results = realNull.stats,
                 real.occs = e@occs,
                 real.occs.grp = e@occs.grp,
                 real.bg = e@bg,
                 real.bg.grp = e@bg.grp)

  # optionally remove temp directory for maxent.jar
  if(e@algorithm == "maxent.jar" & removeMxTemp == TRUE) unlink(tmpdir, recursive = TRUE, force = TRUE)

  
  timed <- proc.time() - start.time
  t.min <- floor(timed[3] / 60)
  t.sec <- timed[3] - (t.min * 60)
  if(quiet == FALSE) message(paste("\nENMnullSims completed in", t.min, "minutes", round(t.sec, 1), "seconds."))
  
  return(e.n)
}
