#' @export
#'

# for split evaluation, label training occs "1" and independent evaluation occs "2" in partitions
nullENMs <- function(e, mod.settings, no.iter, envs = NULL, user.enm = NULL, user.envs.partition = NULL,
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

  # assign the number of cross validation iterations
  nk <- ifelse(e@partition.method == "independent", 1, max(as.numeric(e@occ.grp)))
  # get names of environmental predictor variables
  envs.names <- names(e@occs)[3:ncol(e@occs)]

  # get number of occurrence points by partition
  occs.grp.tbl <- table(e@occ.grp)

  eval.type <- switch(e@partition.method,
                      random = "kfold",
                      jackknife = "kfold",
                      block = "kspatial",
                      checkerboard1 = "kspatial",
                      checkerboard2 = "kspatial",
                      independent = "kfold")

  # if envs were input, get partition groups for envs if using spatial cross validation
  if(!is.null(envs)) {
    envs.pts <- as.data.frame(na.omit(raster::rasterToPoints(envs)))
    names(envs.pts)[1:2] <- names(e@occs)[1:2]
    # convert any factor columns to factor in envs pts dataset
    envs.fact <- which(sapply(e@occs, is.factor))
    envs.pts[,envs.fact] <- factor(envs.pts[,envs.fact])
    # if using a native ENMeval spatial cross validation partitioning method, use it to assign partition groups to envs
    if(e@partition.method != "user") {
      envs.pts$grp <- switch(e@partition.method,
                         block = get.block(e@occs, envs.pts)$bg.grp,
                         checkerboard1 = get.checkerboard1(e@occs, envs.pts, aggregation.factor = strsplit(e@partition.settings, split = "=")[[1]][2]$bg.grp),
                         checkerboard2 = get.checkerboard2(e@occs, envs.pts, aggregation.factor = strsplit(e@partition.settings, split = "=")[[1]][2]$bg.grp),
                         independent = 1,
                         randomkfold = 1,
                         jackknife = 1)
    }else{
      # if user partitions, assign envs partition groups based on user input
      envs.pts$grp <- user.envs.partition
    }
  }

  # # convert fields for categorical data to factor class
  # if(!is.null(categoricals)) {
  #   for(j in 1:length(categoricals)) {
  #     occs.vals[, categoricals[j]] <- as.factor(occs.vals[, categoricals[j]])
  #     bg.vals[, categoricals[j]] <- as.factor(bg.vals[, categoricals[j]])
  #   }
  # }

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

  # # initialize data frames to collect evaluation stats for each partition
  # # per iteration and their averages
  train.stats.names <- e@results %>% dplyr::select(auc.train:ncol(e@results), -AICc, -delta.AICc, -w.AIC, -nparam) %>% names()
  train.stats <- as.data.frame(matrix(nrow = no.iter, ncol = length(train.stats.names)))
  names(train.stats) <- train.stats.names
  k.stats <- e@results.grp %>% dplyr::select(auc.test:ncol(e@results.grp)) %>% names()

  # all.cnames <- c("auc.train", "auc.test", "auc.diff", "or.min", "or.10")
  # k.cnames <- c("auc.test", "auc.diff", "or.min", "or.10")
  # null.cnames <- c("auc.train", "mean.auc.test", "sd.auc.test",
  #                  "mean.auc.diff", "sd.auc.diff", "mean.or.min", "sd.or.min",
  #                  "mean.or.10", "sd.or.10", "nparam")
  # all.rnames <- c("real.mean", "real.sd", "null.mean", "null.sd", "zscore", "pvalue")
  # all.stats <- data.frame(matrix(nrow = 6, ncol = length(all.cnames),
  #                                dimnames = list(all.rnames, all.cnames)))
  # # real.stats <- data.frame(matrix(nrow = 1, ncol = 2,
  # # dimnames = list(NULL, c("auc.train", "nparam"))))
  # null.stats <- data.frame(matrix(nrow = no.iter, ncol = length(null.cnames),
  #                                 dimnames = list(NULL, null.cnames)))

  ############################## #
  # build real model ####
  ############################## #

  mod.tune.args  <- paste(mod.settings, collapse = "_")
  real.mod <- e@models[[mod.tune.args]]
  real.mod.res <- e@results %>% dplyr::filter(tune.args == mod.tune.args)

  # t3 <- proc.time()
  # message("Building and evaluating real model...")
  # mod.args.real <- model.args(mod.name, mod.args, occs.vals, bg.vals, other.args)
  # mod.real <- do.call(mod.fun, mod.args.real)
  # # calculate training auc
  # all.stats["real.mean", "auc.train"] <- dismo::evaluate(occs.vals, bg.vals, mod.real)@auc
  # # real.stats$nparam <- no.params(mod.real, mod.name)
  # kstats <- data.frame(matrix(nrow = nk, ncol = length(k.cnames),
  #                             dimnames = list(NULL, k.cnames)))
  #
  # if(eval.type == "split") {
  #   kstats[1,] <- evalStats(occs.vals, bg.vals, occs.test.vals, mod.real, abs.auc.diff)
  # }else{
  #   for(k in 1:nk) {
  #     occs.train.real.k <- occs.vals[occs.grp != k, ]
  #     occs.test.real.k <- occs.vals[occs.grp == k, ]
  #     bg.train.k <- bg.vals[bg.grp != k, ]
  #     mod.args.real.k <- model.args(mod.name, mod.args, occs.train.real.k, bg.train.k, other.args)
  #     mod.real.k <- do.call(mod.fun, mod.args.real.k)
  #     kstats[k,] <- evalStats(occs.train.real.k, bg.train.k, occs.test.real.k,
  #                             mod.real.k, abs.auc.diff)
  #     message(sprintf("Completed real model evaluation for partition %i.", k))
  #   }
  # }
  #
  # message(paste0("Real model built and evaluated in ", timeCheck(t3), "."))

  # # fill in rest of real model statistics
  # all.stats[1, 2:5] <- apply(kstats, 2, mean)
  # all.stats[2, 2:5] <- apply(kstats, 2, sd)

  ############################## #
  # build null models ####
  ############################## #

  # initialize list to record stats for null iteration i
  nulls <- list()

  t4 <- proc.time()
  message(sprintf("Building and evaluating %i null SDMs...", no.iter))
  pb <- txtProgressBar(0, no.iter, style = 3)

  for(i in 1:no.iter) {

    null.occs.ik <- list()

    # randomly sample the same number of training occs over each k partition
    # of envs; if kspatial evaluation, only sample over the current spatial
    # partition of envs.vals
    for(k in 1:nk) {
      if(eval.type == "kspatial") {
        # if spatial cross-validation, sample null occurrences only from
        # the environment variable grid cells in partition group k
        envs.pts.k <- envs.pts %>% dplyr::filter(grp == k)
      }else{
        # if nonspatial cross-validation, sample null occurrences from
        # the entire environment variable grid for each partition group k
        envs.pts.k <- envs.pts
      }
      # randomly sample n null occurrences, where n equals the number
      # of real occurrence in partition group k
      s.k <- sample(1:nrow(envs.pts.k), occs.grp.tbl[k])
      null.occs.ik[[k]] <- envs.pts.k[s.k, ]
    }

    # bind rows together to make full null occurrence dataset
    null.occs.i.df <- dplyr::bind_rows(null.occs.ik)
    null.occs.i.vals <- null.occs.i.df %>% dplyr::select(-grp)
    # assign the null occurrence partitions as user partition settings, but
    # keep the real model background partitions
    user.grp <- list(occ.grp = null.occs.i.df$grp, bg.grp = e@bg.grp)
    # assign user test partitions to those used in the real model
    user.test.grps <- cbind(e@occs, grp = e@occ.grp)
    # shortcuts for settings
    e.s <- e@other.settings
    e.p <- e@partition.settings

    null.e.i <- ENMevaluate(null.occs.i.vals, bg = e@bg, tune.args = mod.settings, mod.name = e@algorithm, other.args = e.s$other.args, partitions = "user",
                            user.test.grps = user.test.grps, user.grp = user.grp, kfolds = e.p$kfolds, aggregation.factor = e.p$aggregation.factor,
                            doClamp = e.s$doClamp, pred.type = e.s$pred.type, abs.auc.diff = e.s$abs.auc.diff, cbi.eval = e.s$cbi.eval, quiet = TRUE)
    setTxtProgressBar(pb, i)

    nulls[[i]] <- null.e.i@results





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
      # calculate evaluation statistics and put the results in nulls.i
      e <- evalStats(occs.null.train, bg.train, occs.test, mod.k, abs.auc.diff)
      nulls.i[[i]][k,] <- c(i, k, e)

      message(sprintf("Completed partition %i for null model %i.", k, i))
    }
    # average partition statistics
    # for split partition, sd will be NA because input is single fold statistic
    null.stats.means <- apply(nulls.i[[i]][, -c(1, 2)], 2, mean)
    null.stats.sds <- apply(nulls.i[[i]][, -c(1, 2)], 2, sd)
    # alternate the above values for data frame and assign to table
    null.stats[i,2:9] <- c(rbind(null.stats.means, null.stats.sds))
  }
  message(paste0("Null models built and evaluated in ", timeCheck(t4), "."))

  # condense nulls.i into data frame
  nulls.i <- do.call("rbind", nulls.i)
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
                        nulls.i = nulls.i)

  # optionally remove temp directory for maxent.jar
  if(mod.name == "maxent.jar" & removeMxTemp == TRUE) unlink(tmpdir, recursive = T, force = T)

  message(paste0("Null SDM analysis completed in ", timeCheck(t1), "."))

  return(out)
}
