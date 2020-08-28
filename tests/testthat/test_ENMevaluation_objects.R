context("Check ENMevaluation objects")

# read in data
set.seed(48)
occs <- read.csv(file.path(system.file(package="dismo"), "/ex/bradypus.csv"))[,2:3]
envs <- raster::stack(list.files(path=paste(system.file(package='dismo'), '/ex', sep=''), 
                                 pattern='grd', full.names=TRUE))
occs.z <- cbind(occs, raster::extract(envs, occs))
bg <- as.data.frame(dismo::randomPoints(envs, 1000))
names(bg) <- names(occs)
bg.z <- cbind(bg, raster::extract(envs, bg))
# tune.args <- list(fc = c("L","LQ","H"), rm = 1:5)
tune.args.ls <- list("maxnet" = list(fc = c("L","LQ"), rm = 2:3),
                  "boostedRegressionTrees" = list(ntree = 1000, tc = 1:2, lr = 0.01),
                  "randomForest" = list(ntree = 1000, mtry = 4:5))
tune.args.tbl.ls <- lapply(tune.args.ls, expand.grid, stringsAsFactors = FALSE)
parts <- c("block", "checkerboard1", "checkerboard2", "randomkfold", "jackknife", "testing", "none", "user", rep("randomkfold", 5))
kfolds.n <- 4
user.grp <- list(occs.grp = round(runif(nrow(occs), 1, 4)), bg.grp = round(runif(nrow(bg), 1, 4)))
# random sample of occs for runs that use subsets
i <- sample(1:nrow(occs))

e.ls <- list()
e.ls.alg <- c(rep("maxnet", 10), "bioclim", "boostedRegressionTrees", "randomForest")
# maxnet run with tuning parameters, categorical variable, block partitions, and niche overlap
e.ls$block <- ENMevaluate(occs, envs, bg, algorithm = "maxnet", tune.args = tune.args.ls$maxnet, categoricals = "biome", 
                      partitions = "block", overlap = TRUE)
# checkerboard1 partitions
e.ls$cb1 <- ENMevaluate(occs, envs, bg, algorithm = "maxnet", tune.args = tune.args.ls$maxnet, categoricals = "biome", 
                        partitions = "checkerboard1", overlap = TRUE)
# checkerboard2 partitions
e.ls$cb2 <- ENMevaluate(occs, envs, bg, algorithm = "maxnet", tune.args = tune.args.ls$maxnet, categoricals = "biome", 
                        partitions = "checkerboard2", overlap = TRUE)
# random k-fold partitions
e.ls$rand <- ENMevaluate(occs, envs, bg, algorithm = "maxnet", tune.args = tune.args.ls$maxnet, categoricals = "biome", 
                      partitions = "randomkfold", kfolds = kfolds.n, overlap = TRUE)
# jackknife partitions
e.ls$jack <- ENMevaluate(occs[i[1:10],], envs, bg, algorithm = "maxnet", tune.args = tune.args.ls$maxnet, categoricals = "biome", 
                      partitions = "jackknife", overlap = TRUE)
# testing partition
e.ls$ind <- ENMevaluate(occs[i[1:100],], envs, bg, algorithm = "maxnet", tune.args = tune.args.ls$maxnet, categoricals = "biome", 
                      partitions = "testing", occs.testing = occs[i[101:nrow(occs)],], overlap = TRUE)
# no partitions
e.ls$nopart <- ENMevaluate(occs, envs, bg, algorithm = "maxnet", tune.args = tune.args.ls$maxnet, categoricals = "biome", 
                      partitions = "none", overlap = TRUE)
# user partitions
e.ls$user <- ENMevaluate(occs, envs, bg, algorithm = "maxnet", tune.args = tune.args.ls$maxnet, categoricals = "biome", 
                      user.grp = user.grp, partitions = "user", overlap = TRUE)
# no envs (SWD)
e.ls$swd <- ENMevaluate(occs.z, bg = bg.z, algorithm = "maxnet", tune.args = tune.args.ls$maxnet, categoricals = "biome", 
                         partitions = "randomkfold", kfolds = kfolds.n)
# no bg
e.ls$nobg <- ENMevaluate(occs, envs, algorithm = "maxnet", tune.args = tune.args.ls$maxnet, categoricals = "biome", 
                        partitions = "randomkfold", kfolds = kfolds.n, overlap = TRUE)
# bioclim
e.ls$bioclim <- ENMevaluate(occs, envs, bg, algorithm = "bioclim", categoricals = "biome", 
                        partitions = "randomkfold", kfolds = kfolds.n, overlap = TRUE)
# boostedRegressionTrees
e.ls$boostedRegressionTrees <- ENMevaluate(occs, envs, bg, algorithm = "boostedRegressionTrees", tune.args = tune.args.ls$boostedRegressionTrees, categoricals = "biome", 
                            partitions = "randomkfold", kfolds = kfolds.n, overlap = TRUE)
# randomForest
e.ls$randomForest <- ENMevaluate(occs, envs, bg, algorithm = "randomForest", tune.args = tune.args.ls$randomForest, categoricals = "biome", 
                        partitions = "randomkfold", kfolds = kfolds.n, overlap = TRUE)

test_that("ENMevaluation object and slots exist", {
  for(x in 1:length(e.ls)) {
    e <- e.ls[[x]]
    expect_true(!is.null(e))
    expect_true(!is.null(e@algorithm))
    expect_true(!is.null(e@tune.settings))
    expect_true(!is.null(e@partition.method))
    expect_true(!is.null(e@results))
    expect_true(!is.null(e@results.partitions))
    expect_true(!is.null(e@models))
    expect_true(!is.null(e@predictions))
    if(names(e.ls)[x] != "swd") {
      expect_true(raster::nlayers(e@predictions) > 0)  
    }else{
      expect_true(raster::nlayers(e@predictions) == 0)  
    }
    expect_true(!is.null(e@occs))
    expect_true(!is.null(e@occs.grp))
    expect_true(!is.null(e@bg))
    expect_true(!is.null(e@bg.grp))
    expect_true(!is.null(e@overlap))
  }
})

test_that("Data in ENMevaluation object slots have correct form", {
  for(x in 1:length(e.ls)) {
    e <- e.ls[[x]]
    m <- ifelse(names(e.ls)[x] == "boostedRegressionTrees", "boostedRegressionTrees", ifelse(names(e.ls)[x] == "randomForest", "randomForest", "maxnet"))
    # algorithm
    expect_true(e@algorithm == e.ls.alg[x])
    # partition method 
    expect_true(e@partition.method == parts[x])
    # these checks relate to tune.args, which is NULL for BIOCLIM
    if(x != 11) {
      # tune.settings 
      expect_true(all(e@tune.settings[,1:ncol(tune.args.tbl.ls[[m]])] == tune.args.tbl.ls[[m]]))
      # nrow of results
      expect_true(nrow(e@results) == nrow(tune.args.tbl.ls[[m]]))
      # tune.args column values are concat of tuning parameters columns
      # expect_true(all(apply(e@results[names(tune.args.ls[[m]])], 1, paste, collapse = "_") == as.character(e@results$tune.args)))
      # number of models
      expect_true(length(e@models) == nrow(tune.args.tbl.ls[[m]]))
    }
    # number of rows for occs matches occs.grp
    expect_true(nrow(e@occs) == length(e@occs.grp))
    # number of rows for bg matches bg.grp
    expect_true(nrow(e@bg) == length(e@bg.grp))
    # no overlap is calculated for no tuning or BIOCLIM
    if(!(x %in% c(9,11))) {
      # both indicies exist for overlap
      expect_true(length(e@overlap) == 2)
      # number of rows of overlap D matches tune.args
      expect_true(nrow(e@overlap$D) == nrow(tune.args.tbl.ls[[m]]))
      # number of rows of overlap I matches tune.args
      expect_true(nrow(e@overlap$I) == nrow(tune.args.tbl.ls[[m]]))  
    }else{
      # no overlap matrix
      expect_true(length(e@overlap) == 0)
    }
  }
})

# check NA env records

test_that("Records with missing environmental values were removed", {
  expect_true(sum(is.na(e.ls$block@occs)) == 0)
  expect_true(sum(is.na(e.ls$block@bg)) == 0)
})

# check partition numbers

test_that("Spatial block has correct number of partitions", {
  expect_true(length(unique(e.ls$block@occs.grp)) == 4)
  expect_true(length(unique(e.ls$block@bg.grp)) == 4)
})

test_that("Checkerboard 1 has correct number of partitions", {
  expect_true(length(unique(e.ls$cb1@occs.grp)) == 2)
  expect_true(length(unique(e.ls$cb1@bg.grp)) == 2)
})

test_that("Checkerboard 2 has correct number of partitions", {
  expect_true(length(unique(e.ls$cb2@occs.grp)) == 4)
  expect_true(length(unique(e.ls$cb2@bg.grp)) == 4)
})

test_that("Random k-fold has correct number of partitions", {
  expect_true(length(unique(e.ls$rand@occs.grp)) == kfolds.n)
  expect_true(length(unique(e.ls$rand@bg.grp)) == 1)
})

test_that("Jackknife has correct number of partitions", {
  expect_true(length(unique(e.ls$jack@occs.grp)) == nrow(e.ls$jack@occs))
  expect_true(length(unique(e.ls$jack@bg.grp)) == 1)
})

test_that("Testing data has correct number of partitions", {
  expect_true(length(unique(e.ls$ind@occs.grp)) == 1)
  expect_true(length(unique(e.ls$ind@bg.grp)) == 1)
})

test_that("No partitions has no partitions", {
  expect_true(length(unique(e.ls$nopart@occs.grp)) == 1)
  expect_true(length(unique(e.ls$nopart@bg.grp)) == 1)
})

test_that("User has correct number of partitions", {
  expect_true(length(unique(e.ls$user@occs.grp)) == length(unique(user.grp$occs.grp)))
  expect_true(length(unique(e.ls$user@bg.grp)) == length(unique(user.grp$bg.grp)))
})

# check results of specific parameterizations

test_that("Block test results table has correct form", {
  block.res <- e.ls$block@results.partitions
  expect_true(nrow(block.res) == 4 * nrow(tune.args.tbl.ls$maxnet))
  expect_true(max(block.res$fold) == 4)
  expect_true(sum(is.na(block.res)) == 0)
})

test_that("Checkerboard1 test results table has correct form", {
  cb1.res <- e.ls$cb1@results.partitions
  expect_true(nrow(cb1.res) == 2 * nrow(tune.args.tbl.ls$maxnet))
  expect_true(max(cb1.res$fold) == 2)
  expect_true(sum(is.na(cb1.res)) == 0)
})

test_that("Checkerboard2 test results table has correct form", {
  cb2.res <- e.ls$cb2@results.partitions
  expect_true(nrow(cb2.res) == 4 * nrow(tune.args.tbl.ls$maxnet))
  expect_true(max(cb2.res$fold) == 4)
  expect_true(sum(is.na(cb2.res)) == 0)
})

test_that("Random test results table has correct form", {
  rand.res <- e.ls$rand@results.partitions
  expect_true(nrow(rand.res) == kfolds.n * nrow(tune.args.tbl.ls$maxnet))
  expect_true(max(rand.res$fold) == kfolds.n)
  expect_true("cbi.val" %in% names(rand.res))
  expect_true(sum(is.na(rand.res)) == 0)
})

test_that("Jackknife test results table has correct form", {
  jack.res <- e.ls$jack@results.partitions
  noccs <- nrow(e.ls$jack@occs)
  expect_true(nrow(jack.res) == noccs * nrow(tune.args.tbl.ls$maxnet))
  expect_true(max(jack.res$fold) == noccs)
  expect_true(sum(is.na(jack.res$cbi.val)) == 36)
  expect_true(sum(is.na(jack.res)) == 36)
})

test_that("Testing data test results table has correct form", {
  ind.res <- e.ls$ind@results.partitions
  expect_true(nrow(ind.res) == nrow(tune.args.tbl.ls$maxnet))
  expect_true(max(ind.res$fold) == 1)
  expect_true("cbi.val" %in% names(ind.res))
  expect_true(sum(is.na(ind.res)) == 0)
})

test_that("No partition has no test results table", {
  no.res <- e.ls$nopart@results.partitions
  expect_true(nrow(no.res) == 0)
})

