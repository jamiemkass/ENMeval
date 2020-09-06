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


# functions 
test_all <- function(e, alg, parts, tune.args, nparts.occs, nparts.bg, type = "") {
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
      # these checks relate to tune.args, which is NULL for BIOCLIM
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
      expect_true(max(e@results.partitions$fold) == nparts.occs)
      # jackknife has NAs for cbi.val
      if(parts == "jackknife") {
        expect_true(sum(is.na(e@results.partitions)) == nrow(e@results.partitions))
      }else{
        expect_true(sum(is.na(e@results.partitions)) == 0)
      }
    }
  })
}


algs <- list(maxnet = list(fc = c("L","LQ"), rm = 2:3),
             bioclim = list(tails = c("low", "high", "both")),
             boostedRegressionTrees = list(tc = 1:2, lr = 0.01),
             randomForest = list(ntree = 1000, mtry = 4:5))

for(i in 1:length(algs)) {
  
  alg <- names(algs)[i]
  targs <- algs[[i]]
  
  # block partitions
  context(paste("Testing", alg, "with block partitions..."))
  e <- ENMevaluate(occs, envs, bg, algorithm = alg, tune.args = targs, categoricals = "biome", partitions = "block", overlap = TRUE, quiet = TRUE)
  # ns <- ENMnullSims(e, mod.settings = targs, no.iter = 10)
  test_all(e, alg, "block", targs, 4, 4)
  
  # checkerboard1 partitions
  context(paste("Testing", alg, "with checkerboard1 partitions..."))
  e <- ENMevaluate(occs, envs, bg, algorithm = alg, tune.args = targs, categoricals = "biome", partitions = "checkerboard1", overlap = TRUE, quiet = TRUE)
  test_all(e, alg, "checkerboard1", targs, 2, 2)
  
  # checkerboard2 partitions
  context(paste("Testing", alg, "with checkerboard2 partitions..."))
  e <- ENMevaluate(occs, envs, bg, algorithm = alg, tune.args = targs, categoricals = "biome", partitions = "checkerboard2", overlap = TRUE, quiet = TRUE)
  test_all(e, alg, "checkerboard2", targs, 4, 4)
  
  # random k-fold partitions
  context(paste("Testing", alg, "with random 4-fold partitions..."))
  e <- ENMevaluate(occs, envs, bg, algorithm = alg, tune.args = targs, categoricals = "biome", partitions = "randomkfold", kfolds = 4, overlap = TRUE, quiet = TRUE)
  test_all(e, alg, "randomkfold", targs, 4, 1)
  
  # jackknife partitions
  context(paste("Testing", alg, "with jackknife partitions..."))
  e <- ENMevaluate(occs[1:10,], envs, bg, algorithm = alg, tune.args = targs, categoricals = "biome", partitions = "jackknife", overlap = TRUE, quiet = TRUE)
  test_all(e, alg, "jackknife", targs, nrow(e@occs), 1)
  
  # testing partition
  context(paste("Testing", alg, "with testing partitions..."))
  e <- ENMevaluate(occs[1:100,], envs, bg, algorithm = alg, tune.args = targs, categoricals = "biome", partitions = "testing", occs.testing = occs[101:nrow(occs),], overlap = TRUE, quiet = TRUE)
  test_all(e, alg, "testing", targs, 1, 1)
  
  # no partitions
  context(paste("Testing", alg, "with no partitions..."))
  e <- ENMevaluate(occs, envs, bg, algorithm = alg, tune.args = targs, categoricals = "biome", partitions = "none", overlap = TRUE, quiet = TRUE)
  test_all(e, alg, "none", targs, 1, 1)
  
  # user partitions
  context(paste("Testing", alg, "with user partitions..."))
  user.grp <- list(occs.grp = round(runif(nrow(occs), 1, 4)), bg.grp = round(runif(nrow(bg), 1, 4)))
  e <- ENMevaluate(occs, envs, bg, algorithm = alg, tune.args = targs, categoricals = "biome", user.grp = user.grp, partitions = "user", overlap = TRUE, quiet = TRUE)
  test_all(e, alg, "user", targs, 4, 4)
  
  # no envs (SWD)
  context(paste("Testing", alg, "with random 4-fold partitions and no raster environmental variables..."))
  e <- ENMevaluate(occs.z, bg = bg.z, algorithm = alg, tune.args = targs, categoricals = "biome", partitions = "randomkfold", kfolds = 4, quiet = TRUE)
  test_all(e, alg, "randomkfold", targs, 4, 1, type = "swd")
  
  # no bg
  context(paste("Testing", alg, "with random 4-fold partitions and no input background data..."))
  e <- ENMevaluate(occs, envs, algorithm = alg, n.bg = 1000, tune.args = targs, categoricals = "biome", partitions = "randomkfold", kfolds = 4, overlap = TRUE, quiet = TRUE)
  test_all(e, alg, "randomkfold", targs, 4, 1) 
}


