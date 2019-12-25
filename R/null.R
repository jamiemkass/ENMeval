#' @export
#'

# for split evaluation, label training occs "1" and independent evaluation occs "2" in partitions
nullENMs <- function(e, no.iter, envs = NULL,
                     eval.type = c("split", "kfold", "kspatial"),
                     categoricals = NULL, envs.grp = NULL, abs.auc.diff = TRUE,
                     other.args = NULL, rasterBrick = TRUE, removeMxTemp = TRUE) {

  # record start time
  t1 <- proc.time()
  message("Beginning null SDM analysis...")

  # changing stack to brick to speed up analysis
  if(!is.null(envs) & class(envs) != "RasterBrick" & rasterBrick == TRUE) envs <- raster::brick(envs)

  # # if split evaluation, partition number equals 1
  # if(eval.type == "split") {
  #   nk <- 1
  #   occs.grp <- rep(0, nrow(occs))
  #   bg.grp <- rep(0, nrow(bg))
  # }else{
  #   # get total number of partitions (k)
  #   nk <- length(unique(occs.grp))
  # }

  nk <- max(as.numeric(e@occ.grp))

  # get number of occurrence points by partition
  occs.grp.tbl <- table(e@occ.grp)

  # get environmental values for occs and bg
  t2 <- proc.time()

  if(ncol(e@occ.pts > 2) & ncol(bg) > 2 & is.null(envs)) {
    message("Environmental values provided for occurrence and background records. Skipping extraction...")
    occs.vals <- occs
    bg.vals <- bg
    envs.vals <- bg
    if(!is.null(occs.indTest)) {
      if(ncol(occs.indTest) < 3) stop("Please insert fields for environmental values of occurrence test points.")
      occs.test.vals <- occs.indTest
    }
  }else{
    message("Extracting environmental values...")
    colnames(bg) <- colnames(occs)
    if(!is.null(occs.indTest)) {
      pt.vals <- as.data.frame(raster::extract(envs, rbind(occs, bg, occs.indTest)))
      occs.test.vals <- pt.vals[seq(nrow(occs) + nrow(bg) + 1, nrow(pt.vals)), ]
    }else{
      pt.vals <- as.data.frame(raster::extract(envs, rbind(occs, bg)))
    }
    occs.vals <- pt.vals[seq(1, nrow(occs)), ]
    bg.vals <- pt.vals[seq(nrow(occs)+1, nrow(occs) + nrow(bg)), ]
    envs.vals <- as.data.frame(raster::getValues(envs))
    # now remove NAs from envs.vals
    envs.vals <- na.omit(envs.vals)
    message(paste0("Environmental values extracted in ", timeCheck(t2), "."))
  }

  # convert fields for categorical data to factor class
  if(!is.null(categoricals)) {
    for(j in 1:length(categoricals)) {
      occs.vals[, categoricals[j]] <- as.factor(occs.vals[, categoricals[j]])
      bg.vals[, categoricals[j]] <- as.factor(bg.vals[, categoricals[j]])
    }
  }

  # get model function name (only Maxent functions available now)
  if(mod.name == "maxent.jar") {
    mod.fun <- dismo::maxent
    # create temp directory to store maxent.jar output, for potential removal later
    tmpdir <- paste(tempdir(), runif(1,0,1), sep = "/")
    dir.create(tmpdir, showWarnings = TRUE, recursive = FALSE)
  }else if(mod.name == "maxnet") {
    mod.fun <- maxnet::maxnet
  }else{
    message('Only Maxent functions available now. Please choose either "maxent.jar" or "maxnet".')
    return()
  }

  # initialize data frames to collect evaluation stats for each partition
  # per iteration and their averages
  all.cnames <- c("auc.train", "auc.test", "auc.diff", "or.min", "or.10")
  k.cnames <- c("auc.test", "auc.diff", "or.min", "or.10")
  null.cnames <- c("auc.train", "mean.auc.test", "sd.auc.test",
                   "mean.auc.diff", "sd.auc.diff", "mean.or.min", "sd.or.min",
                   "mean.or.10", "sd.or.10", "nparam")
  all.rnames <- c("real.mean", "real.sd", "null.mean", "null.sd", "zscore", "pvalue")
  all.stats <- data.frame(matrix(nrow = 6, ncol = length(all.cnames),
                                 dimnames = list(all.rnames, all.cnames)))
  # real.stats <- data.frame(matrix(nrow = 1, ncol = 2,
  # dimnames = list(NULL, c("auc.train", "nparam"))))
  null.stats <- data.frame(matrix(nrow = no.iter, ncol = length(null.cnames),
                                  dimnames = list(NULL, null.cnames)))

  ############################## #
  # build real model ####
  ############################## #

  t3 <- proc.time()
  message("Building and evaluating real model...")
  mod.args.real <- model.args(mod.name, mod.args, occs.vals, bg.vals, other.args)
  mod.real <- do.call(mod.fun, mod.args.real)
  # calculate training auc
  all.stats["real.mean", "auc.train"] <- dismo::evaluate(occs.vals, bg.vals, mod.real)@auc
  # real.stats$nparam <- no.params(mod.real, mod.name)
  kstats <- data.frame(matrix(nrow = nk, ncol = length(k.cnames),
                              dimnames = list(NULL, k.cnames)))

  if(eval.type == "split") {
    kstats[1,] <- evalStats(occs.vals, bg.vals, occs.test.vals, mod.real, abs.auc.diff)
  }else{
    for(k in 1:nk) {
      occs.train.real.k <- occs.vals[occs.grp != k, ]
      occs.test.real.k <- occs.vals[occs.grp == k, ]
      bg.train.k <- bg.vals[bg.grp != k, ]
      mod.args.real.k <- model.args(mod.name, mod.args, occs.train.real.k, bg.train.k, other.args)
      mod.real.k <- do.call(mod.fun, mod.args.real.k)
      kstats[k,] <- evalStats(occs.train.real.k, bg.train.k, occs.test.real.k,
                              mod.real.k, abs.auc.diff)
      message(sprintf("Completed real model evaluation for partition %i.", k))
    }
  }

  message(paste0("Real model built and evaluated in ", timeCheck(t3), "."))

  # fill in rest of real model statistics
  all.stats[1, 2:5] <- apply(kstats, 2, mean)
  all.stats[2, 2:5] <- apply(kstats, 2, sd)

  ############################## #
  # build null models ####
  ############################## #

  # initialize list to record stats for null iteration i
  null.stats.iters <- list()

  t4 <- proc.time()
  message(sprintf("Building and evaluating %i null SDMs...", no.iter))

  for(i in 1:no.iter) {
    # initialize data frame for partition statistics for current iteration
    dn <- list(NULL, c("iter", "grp", k.cnames))
    null.stats.iters[[i]] <- data.frame(matrix(nrow = nk,
                                               ncol = length(k.cnames) + 2,
                                               dimnames = dn))
    # make list to hold different null occurrence datasets
    occs.null <- list()

    # randomly sample the same number of training occs over each k partition
    # of envs; if kspatial evaluation, only sample over the current spatial
    # partition of envs.vals
    for(k in 1:nk) {
      if(eval.type == "kspatial") {
        envs.vals.k <- envs.vals[envs.grp == k,]
      }else{
        envs.vals.k <- envs.vals
      }
      samp <- sample(1:nrow(envs.vals.k), occs.grp.tbl[k])
      occs.null[[k]] <- envs.vals.k[samp, ]
    }

    # convert fields for categorical data to factor class
    if(!is.null(categoricals)) {
      for(k in 1:nk) {
        for(j in 1:length(categoricals)) {
          occs.null[[k]][, categoricals[j]] <- as.factor(occs.null[[k]][, categoricals[j]])
        }
      }
    }

    occs.null.all <- do.call(rbind, occs.null)

    mod.args.i <- model.args(mod.name, mod.args, occs.null.all, bg.vals, other.args)
    mod.i <- do.call(mod.fun, mod.args.i)
    null.stats[i,]$auc.train <- dismo::evaluate(occs.null.all, bg.vals, mod.i)@auc
    null.stats[i,]$nparam <- no.params(mod.i, mod.name)

    # get training and testing data for current iteration
    for(k in 1:nk) {
      if(eval.type == "kfold" | eval.type == "kspatial") {
        # if kfold evaluation, use non-k training and k testing partitions
        occs.null.train <- do.call(rbind, occs.null[-k])
        bg.train <- bg.vals[bg.grp != k, ]
        occs.test <- occs.vals[occs.grp == k, ]
      }else if(eval.type == "split") {
        # if split evaluation, use single random occs dataset and real
        # independent testing dataset
        occs.null.train <- occs.null[[1]]
        bg.train <- bg.vals
        occs.test <- occs.test.vals
      }

      # build custom mod.args for this iteration
      mod.args.k <- model.args(mod.name, mod.args, occs.null.train, bg.train, other.args)
      # build the model
      mod.k <- do.call(mod.fun, mod.args.k)
      # calculate evaluation statistics and put the results in null.stats.iters
      e <- evalStats(occs.null.train, bg.train, occs.test, mod.k, abs.auc.diff)
      null.stats.iters[[i]][k,] <- c(i, k, e)

      message(sprintf("Completed partition %i for null model %i.", k, i))
    }
    # average partition statistics
    # for split partition, sd will be NA because input is single fold statistic
    null.stats.means <- apply(null.stats.iters[[i]][, -c(1, 2)], 2, mean)
    null.stats.sds <- apply(null.stats.iters[[i]][, -c(1, 2)], 2, sd)
    # alternate the above values for data frame and assign to table
    null.stats[i,2:9] <- c(rbind(null.stats.means, null.stats.sds))
  }
  message(paste0("Null models built and evaluated in ", timeCheck(t4), "."))

  # condense null.stats.iters into data frame
  null.stats.iters <- do.call("rbind", null.stats.iters)
  # calculate means and standard deviations of mean k-fold values for null stats
  all.stats["null.mean",] <- apply(null.stats[c(1,2,4,6,8)], 2, mean)
  all.stats["null.sd",] <- apply(null.stats[c(1,2,4,6,8)], 2, sd)
  all.stats["zscore",] <- (all.stats["real.mean",] - all.stats["null.mean",]) / all.stats["null.sd",]
  all.stats["pvalue", 1:2] <- 1 - sapply(all.stats["zscore", 1:2], pnorm)
  all.stats["pvalue", 3:5] <- sapply(all.stats["zscore", 3:5], pnorm)

  if(is.null(model.args)) model.args <- list()
  # condense mod.args to named matrix for inserting into class slot
  out <- nullSDMResults(model = mod.name, model.args = mod.args,
                        eval.type = eval.type, occs = occs,
                        occs.grp = occs.grp, bg = bg, bg.grp = bg.grp,
                        no.iter = no.iter, all.stats = all.stats,
                        null.stats = null.stats,
                        null.stats.iters = null.stats.iters)

  # optionally remove temp directory for maxent.jar
  if(mod.name == "maxent.jar" & removeMxTemp == TRUE) unlink(tmpdir, recursive = T, force = T)

  message(paste0("Null SDM analysis completed in ", timeCheck(t1), "."))

  return(out)
}
