context("Check ENMevaluation objects")

# read in data
set.seed(48)
occs <- readRDS("data/bvariegatus.rds")
envs <- raster::stack(list.files(path=paste(system.file(package='dismo'), '/ex', sep=''), 
                                 pattern='grd', full.names=TRUE))
bg <- as.data.frame(dismo::randomPoints(envs, 1000))
names(bg) <- names(occs)
# tune.args <- list(fc = c("L","LQ","H"), rm = 1:5)
tune.args <- list(fc = "L", rm = 2:3)
tune.args.tbl <- expand.grid(tune.args, stringsAsFactors = FALSE)
parts <- c("block", "checkerboard1", "checkerboard2", "randomkfold", "jackknife", "independent", "none", "user")
kfolds.n <- 4
user.grp <- list(occ.grp = round(runif(nrow(occs), 1, 4)), bg.grp = round(runif(nrow(bg), 1, 4)))
# random sample of occs for runs that use subsets
i <- sample(1:nrow(occs))

e.ls <- list()
# maxnet run with tuning parameters, categorical variable, block partitions, and niche overlap
e.ls$block <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", 
                      partitions = parts[1], overlap = TRUE)
# checkerboard1 partitions
e.ls$cb1 <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", 
                        partitions = parts[2], overlap = TRUE)
# checkerboard2 partitions
e.ls$cb2 <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", 
                        partitions = parts[3], overlap = TRUE)
# random k-fold partitions
e.ls$rand <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", 
                      partitions = parts[4], kfolds = kfolds.n, overlap = TRUE)
# jackknife partitions
e.ls$jack <- ENMevaluate(occs[i[1:10],], envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", 
                      partitions = parts[5], overlap = TRUE)
# independent partition
e.ls$ind <- ENMevaluate(occs[i[1:450],], envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", 
                      partitions = parts[6], occs.ind = occs[i[451:nrow(occs)],], overlap = TRUE)
# no partitions
e.ls$no <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", 
                      partitions = parts[7], overlap = TRUE)
# user partitions
e.ls$user <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", 
                      user.grp = user.grp, partitions = parts[8], overlap = TRUE)

test_that("ENMevaluation object and slots exist", {
  for(x in 1:length(e.ls)) {
    e <- e.ls[[x]]
    expect_true(!is.null(e))
    expect_true(!is.null(e@algorithm))
    expect_true(!is.null(e@tune.settings))
    expect_true(!is.null(e@partition.method))
    expect_true(!is.null(e@results))
    expect_true(!is.null(e@results.grp))
    expect_true(!is.null(e@models))
    expect_true(!is.null(e@predictions))
    expect_true(!is.null(e@occ.pts))
    expect_true(!is.null(e@occ.grp))
    expect_true(!is.null(e@bg.pts))
    expect_true(!is.null(e@bg.grp))
    expect_true(!is.null(e@overlap))
  }
})

test_that("Data in ENMevaluation object slots have correct form", {
  for(x in 1:length(e.ls)) {
    e <- e.ls[[x]]
    # algorithm
    expect_true(e@algorithm == "maxnet")
    # tune.settings 
    expect_true(all(e@tune.settings[,1:2] == tune.args.tbl))
    # partition method 
    expect_true(e@partition.method == parts[x])
    # nrow of results
    expect_true(nrow(e@results) == nrow(tune.args.tbl))
    # tune.args column values are concat of tuning parameters columns
    expect_true(all(apply(e@results[names(tune.args)], 1, paste, collapse = "_") == as.character(e@results$tune.args)))
    # number of models
    expect_true(length(e@models) == nrow(tune.args.tbl))
    # number of rows for occs.pts matches occ.grp
    expect_true(nrow(e@occ.pts) == length(e@occ.grp))
    # number of rows for bg.pts matches bg.grp
    expect_true(nrow(e@bg.pts) == length(e@bg.grp))
    # both indicies exist for overlap
    expect_true(length(e@overlap) == 2)
    # number of rows of overlap D matches tune.args
    expect_true(nrow(e@overlap$D) == nrow(tune.args.tbl))
    # number of rows of overlap I matches tune.args
    expect_true(nrow(e@overlap$I) == nrow(tune.args.tbl))
  }
})

# check NA env records

test_that("Records with missing environmental values were removed", {
  x.occs <- raster::extract(envs, e.ls$block@occ.pts)
  x.bg <- raster::extract(envs, e.ls$block@bg.pts)
  expect_true(sum(is.na(x.occs)) == 0)
  expect_true(sum(is.na(x.bg)) == 0)
})

# check partition numbers

test_that("Spatial block has correct number of partitions", {
  expect_true(length(unique(e.ls$block@occ.grp)) == 4)
  expect_true(length(unique(e.ls$block@bg.grp)) == 4)
})

test_that("Checkerboard 1 has correct number of partitions", {
  expect_true(length(unique(e.ls$cb1@occ.grp)) == 2)
  expect_true(length(unique(e.ls$cb1@bg.grp)) == 2)
})

test_that("Checkerboard 2 has correct number of partitions", {
  expect_true(length(unique(e.ls$cb2@occ.grp)) == 4)
  expect_true(length(unique(e.ls$cb2@bg.grp)) == 4)
})

test_that("Random k-fold has correct number of partitions", {
  expect_true(length(unique(e.ls$rand@occ.grp)) == kfolds.n)
  expect_true(length(unique(e.ls$rand@bg.grp)) == 1)
})

test_that("Jackknife has correct number of partitions", {
  expect_true(length(unique(e.ls$jack@occ.grp)) == 10)
  expect_true(length(unique(e.ls$jack@bg.grp)) == 1)
})

test_that("Independent has correct number of partitions", {
  expect_true(length(unique(e.ls$ind@occ.grp)) == 2)
  expect_true(length(unique(e.ls$ind@bg.grp)) == 1)
})

test_that("No partitions has no partitions", {
  expect_true(length(unique(e.ls$no@occ.grp)) == 1)
  expect_true(length(unique(e.ls$no@bg.grp)) == 1)
})

test_that("User has correct number of partitions", {
  expect_true(length(unique(e.ls$user@occ.grp)) == length(unique(user.grp$occ.grp)))
  expect_true(length(unique(e.ls$user@bg.grp)) == length(unique(user.grp$bg.grp)))
})

# check results of specific parameterizations

test_that("Independent test results table has correct form", {
  ind.res <- e.ls$ind@results.grp
  expect_true(nrow(ind.res) == nrow(tune.args.tbl))
  expect_true("cbi.test" %in% names(ind.res))
})

test_that("No partition has no test results table", {
  no.res <- e.ls$no@results.grp
  expect_true(nrow(no.res) == 0)
})
