library(dplyr)
options(warn=-1)

# read in data
set.seed(48)
occs <- read.csv(file.path(system.file(package="dismo"), "/ex/bradypus.csv"))[,2:3]
envs <- raster::stack(list.files(path=paste(system.file(package='dismo'), '/ex', sep=''), 
                                 pattern='grd', full.names=TRUE))
occs.z <- cbind(occs, raster::extract(envs, occs))
occs.z$biome <- factor(occs.z$biome)
bg <- as.data.frame(dismo::randomPoints(envs, 1000))
names(bg) <- names(occs)
bg.z <- cbind(bg, raster::extract(envs, bg))
bg.z$biome <- factor(bg.z$biome)

# testing functions 
test_ENMevaluation <- function(e, alg, parts, tune.args, nparts.occs, nparts.bg, type = "") {
  tune.tbl <- expand.grid(tune.args, stringsAsFactors = FALSE)
  test_that("ENMevaluation object and slots exist", {
    expect_true(!is.null(e))
    expect_true(!is.null(e@algorithm))
    expect_true(!is.null(e@tune.settings))
    expect_true(!is.null(e@partition.method))
    expect_true(!is.null(e@results))
    expect_true(!is.null(e@results.partitions))
    expect_true(!is.null(e@models))
    expect_true(!is.null(e@predictions))
    if(type == "swd") {
      expect_true(raster::nlayers(e@predictions) == 0)  
    }else{
      expect_true(raster::nlayers(e@predictions) > 0)  
    }
    expect_true(!is.null(e@occs))
    expect_true(!is.null(e@occs.grp))
    expect_true(!is.null(e@bg))
    expect_true(!is.null(e@bg.grp))
    expect_true(!is.null(e@overlap))
  })  
  
  test_that("Data in ENMevaluation object slots have correct form", {
    # algorithm
    expect_true(e@algorithm == alg)
    # partition method 
    expect_true(e@partition.method == parts)
    # these checks relate to tune.args, which may be NULL
    if(!is.null(tune.args)) {
      # tune.settings 
      expect_true(all(as.data.frame(e@tune.settings[,1:ncol(tune.tbl)]) == tune.tbl))
      # nrow of results
      expect_true(nrow(e@results) == nrow(tune.tbl))
      # tune.args column values are concat of tuning parameters columns
      # expect_true(all(apply(e@results[names(tune.args.ls[[m]])], 1, paste, collapse = "_") == as.character(e@results$tune.args)))
      # number of models
      expect_true(length(e@models) == nrow(tune.tbl))
    }
    # number of rows for occs matches occs.grp
    expect_true(nrow(e@occs) == length(e@occs.grp))
    # number of rows for bg matches bg.grp
    expect_true(nrow(e@bg) == length(e@bg.grp))
    # no overlap is calculated for no tuning or BIOCLIM
    if(length(e@overlap) > 0) {
      # both indices exist for overlap
      expect_true(length(e@overlap) == 2)
      # number of rows of overlap D matches tune.args
      expect_true(nrow(e@overlap$D) == nrow(tune.tbl))
      # number of rows of overlap I matches tune.args
      expect_true(nrow(e@overlap$I) == nrow(tune.tbl))  
    }else{
      # no overlap matrix
      expect_true(length(e@overlap) == 0)
    }
  })
  
  test_that("Records with missing environmental values were removed", {
    expect_true(sum(is.na(e@occs)) == 0)
    expect_true(sum(is.na(e@bg)) == 0)
  })
  
  test_that("Number of partitions is correct", {
    expect_true(length(unique(e@occs.grp)) == nparts.occs)
    expect_true(length(unique(e@bg.grp)) == nparts.bg)
  })
  
  
  test_that("Results table for partitions has correct form", {
    if(parts == "none") {
      expect_true(nrow(e@results.partitions) == 0)
    }else{
      expect_true(nrow(e@results.partitions) == nparts.occs * nrow(tune.tbl))
      if(parts != "testing") {
        expect_true(max(e@results.partitions$fold) == nparts.occs)
      }else{
        expect_true(max(e@results.partitions$fold) == 0)
      }
      # jackknife has NAs for cbi.val
      if(parts == "jackknife") {
        expect_true(sum(is.na(e@results.partitions)) == nrow(e@results.partitions))
      }else{
        expect_true(sum(is.na(e@results.partitions)) == 0)
      }
    }
  })
}

test_ENMnullSims <- function(e, ns, no.iter, alg, parts, mod.settings, nparts.occs, nparts.bg, n.sims, type = "") {
  mod.settings.tbl <- expand.grid(mod.settings)
  test_that("ENMnullSims object and slots exist", {
    expect_true(!is.null(ns))
    expect_true(!is.null(ns@null.algorithm))
    expect_true(!is.null(ns@null.mod.settings))
    expect_true(!is.null(ns@null.partition.method))
    expect_true(!is.null(ns@null.partition.settings))
    expect_true(!is.null(ns@null.other.settings))
    expect_true(!is.null(ns@no.iter))
    expect_true(!is.null(ns@null.results))
    expect_true(!is.null(ns@null.results.partitions))
    expect_true(!is.null(ns@real.vs.null.results))
    expect_true(!is.null(ns@real.occs))
    expect_true(!is.null(ns@real.occs.grp))
    expect_true(!is.null(ns@real.bg))
    expect_true(!is.null(ns@real.bg.grp))
  })  
  
  test_that("Data in ENMnullSims object slots have correct form", {
    # algorithm
    expect_true(ns@null.algorithm == alg)
    # partition method 
    expect_true(ns@null.partition.method == parts)
    # mod.settings 
    if(ncol(mod.settings.tbl) > 1) {
      expect_true(all(ns@null.mod.settings[,1:ncol(mod.settings.tbl)] == mod.settings.tbl))  
    }else{
      expect_true(as.character(ns@null.mod.settings[,1]) == mod.settings.tbl)  
    }
    # no. of iterations
    expect_true(ns@no.iter == no.iter)
    # number of rows in results table
    expect_true(nrow(ns@null.results) == no.iter)
    # number of rows in results table for partitions
    if(ns@null.partition.method == "none") {
      expect_true(nrow(ns@null.results.partitions) == 0)
    }else{
      expect_true(nrow(ns@null.results.partitions) == no.iter * nparts.occs)  
    }
    
    # number of rows in real vs null results table
    expect_true(nrow(ns@real.vs.null.results) == 6)
    # there should only be two NA values for this table: read.sd for auc.train and cbi.train
    if(parts == "jackknife") {
      expect_true(sum(is.na(ns@real.vs.null.results[2,])) == 3)
      expect_true(sum(is.na(ns@real.vs.null.results[,6])) == 6) 
    }else if(parts == "testing") {
      expect_true(sum(is.na(ns@real.vs.null.results[2,])) == 7) 
    }else{
      expect_true(sum(is.na(ns@real.vs.null.results[2,])) == 2)  
    }
    # check that tables match
    expect_true(all(ns@real.occs == e@occs))
    expect_true(all(ns@real.bg == e@bg))
    expect_true(all(ns@real.occs.grp == e@occs.grp))
  })
}

algs <- list(maxnet = list(fc = c("L","LQ"), rm = 2:3),
             bioclim = list(tails = c("low", "high", "both")),
             boostedRegressionTrees = list(tc = 2:3, lr = 0.001),
             randomForest = list(ntree = 1000, mtry = 4:5))

no.iter <- 5

for(i in 1:length(algs)) {
  
  alg <- names(algs)[i]
  targs <- algs[[i]]
  mset <- lapply(targs, function(x) x[1])
  
  # block partitions
  context(paste("Testing ENMevaluate for", alg, "with block partitions..."))
  e <- ENMevaluate(occs, envs, bg, algorithm = alg, tune.args = targs, categoricals = "biome", partitions = "block", overlap = TRUE, quiet = TRUE)
  test_ENMevaluation(e, alg, "block", targs, 4, 4)
  context(paste("Testing ENMnullSims for", alg, "with block partitions..."))
  ns <- ENMnullSims(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
  test_ENMnullSims(e, ns, no.iter, alg, "block", mset, 4, 4)
  
  # checkerboard1 partitions
  context(paste("Testing ENMevaluate for", alg, "with checkerboard1 partitions..."))
  e <- ENMevaluate(occs, envs, bg, algorithm = alg, tune.args = targs, categoricals = "biome", partitions = "checkerboard1", overlap = TRUE, quiet = TRUE)
  test_ENMevaluation(e, alg, "checkerboard1", targs, 2, 2)
  context(paste("Testing ENMnullSims for", alg, "with checkerboard1 partitions..."))
  ns <- ENMnullSims(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
  test_ENMnullSims(e, ns, no.iter, alg, "checkerboard1", mset, 2, 2)
  
  # checkerboard2 partitions
  context(paste("Testing ENMevaluate for", alg, "with checkerboard2 partitions..."))
  e <- ENMevaluate(occs, envs, bg, algorithm = alg, tune.args = targs, categoricals = "biome", partitions = "checkerboard2", overlap = TRUE, quiet = TRUE)
  test_ENMevaluation(e, alg, "checkerboard2", targs, 4, 4)
  context(paste("Testing ENMnullSims for", alg, "with checkerboard2 partitions..."))
  ns <- ENMnullSims(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
  test_ENMnullSims(e, ns, no.iter, alg, "checkerboard2", mset, 4, 4)
  
  # random k-fold partitions
  context(paste("Testing ENMevaluate for", alg, "with random 4-fold partitions..."))
  e <- ENMevaluate(occs, envs, bg, algorithm = alg, tune.args = targs, categoricals = "biome", partitions = "randomkfold", overlap = TRUE, quiet = TRUE)
  test_ENMevaluation(e, alg, "randomkfold", targs, 5, 1)
  context(paste("Testing ENMnullSims for", alg, "with random 4-fold partitions..."))
  ns <- ENMnullSims(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
  test_ENMnullSims(e, ns, no.iter, alg, "randomkfold", mset, 5, 1)
  
  # jackknife partitions
  context(paste("Testing ENMevaluate for", alg, "with jackknife partitions..."))
  e <- ENMevaluate(occs[1:10,], envs, bg, algorithm = alg, tune.args = targs, partitions = "jackknife", overlap = TRUE, quiet = TRUE)
  test_ENMevaluation(e, alg, "jackknife", targs, nrow(e@occs), 1)
  context(paste("Testing ENMnullSims for", alg, "with jackknife partitions..."))
  ns <- ENMnullSims(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
  test_ENMnullSims(e, ns, no.iter, alg, "jackknife", mset, nrow(e@occs), 1)
  
  # testing partition
  context(paste("Testing ENMevaluate for", alg, "with testing partitions..."))
  e <- ENMevaluate(occs[1:100,], envs, bg, algorithm = alg, tune.args = targs, categoricals = "biome", partitions = "testing", occs.testing = occs[101:nrow(occs),], overlap = TRUE, quiet = TRUE)
  test_ENMevaluation(e, alg, "testing", targs, 1, 1)
  context(paste("Testing ENMnullSims for", alg, "with testing partitions..."))
  ns <- ENMnullSims(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
  test_ENMnullSims(e, ns, no.iter, alg, "testing", mset, 1, 1)
  
  # no partitions
  context(paste("Testing ENMevaluate for", alg, "with no partitions..."))
  e <- ENMevaluate(occs, envs, bg, algorithm = alg, tune.args = targs, categoricals = "biome", partitions = "none", overlap = TRUE, quiet = TRUE)
  test_ENMevaluation(e, alg, "none", targs, 1, 1)
  context(paste("Testing ENMnullSims for", alg, "with no partitions..."))
  ns <- ENMnullSims(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
  test_ENMnullSims(e, ns, no.iter, alg, "none", mset, 1, 1)
  
  # user partitions
  context(paste("Testing ENMevaluate for", alg, "with user partitions..."))
  user.grp <- list(occs.grp = round(runif(nrow(occs), 1, 4)), bg.grp = round(runif(nrow(bg), 1, 4)))
  e <- ENMevaluate(occs, envs, bg, algorithm = alg, tune.args = targs, categoricals = "biome", user.grp = user.grp, partitions = "user", overlap = TRUE, quiet = TRUE)
  test_ENMevaluation(e, alg, "user", targs, 4, 4)
  context(paste("Testing ENMnullSims for", alg, "with user partitions..."))
  ns <- ENMnullSims(e, mod.settings = mset, no.iter = no.iter, user.eval.type = "kspatial", quiet = TRUE)
  test_ENMnullSims(e, ns, no.iter, alg, "user", mset, 4, 4)
  
  # no envs (SWD)
  context(paste("Testing ENMevaluate for", alg, "with random 4-fold partitions and no raster environmental variables..."))
  e <- ENMevaluate(occs.z, bg = bg.z, algorithm = alg, tune.args = targs, categoricals = "biome", partitions = "randomkfold", quiet = TRUE)
  test_ENMevaluation(e, alg, "randomkfold", targs, 5, 1, type = "swd")
  context(paste("Testing ENMnullSims for", alg, "with random 4-fold partitions and no raster environmental variables..."))
  ns <- ENMnullSims(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
  test_ENMnullSims(e, ns, no.iter, alg, "randomkfold", mset, 5, 1)
  
  # no bg
  context(paste("Testing ENMevaluate for", alg, "with random 4-fold partitions and no input background data..."))
  e <- ENMevaluate(occs, envs, algorithm = alg, n.bg = 1000, tune.args = targs, categoricals = "biome", partitions = "randomkfold", overlap = TRUE, quiet = TRUE)
  test_ENMevaluation(e, alg, "randomkfold", targs, 5, 1) 
  context(paste("Testing ENMnullSims for", alg, "with random 4-fold partitions and no input background data..."))
  ns <- ENMnullSims(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
  test_ENMnullSims(e, ns, no.iter, alg, "randomkfold", mset, 5, 1)
}


