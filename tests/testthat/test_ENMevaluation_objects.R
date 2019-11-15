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
parts <- c("block", "checkerboard1", "checkerboard2", "randomkfold", "jackknife", "independent", "none")
# random sample of occs for runs that use subsets
i <- sample(1:nrow(occs))

e.ls <- list()
# maxnet run with tuning parameters, categorical variable, block partitions, and niche overlap
e.ls$b <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", 
                      partitions = parts[1], overlap = TRUE)
# checkerboard1 partitions
e.ls$cb1 <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", 
                        partitions = parts[2], overlap = TRUE)
# checkerboard2 partitions
e.ls$cb2 <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", 
                        partitions = parts[3], overlap = TRUE)
# random k-fold partitions
e.ls$r <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", 
                      partitions = parts[4], kfolds = 4, overlap = TRUE)
# jackknife partitions
e.ls$j <- ENMevaluate(occs[i[1:10],], envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", 
                      partitions = parts[5], overlap = TRUE)
# independent partition
e.ls$i <- ENMevaluate(occs[i[1:450],], envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", 
                      partitions = parts[6], occs.ind = occs[i[451:nrow(occs)],], overlap = TRUE)
# no partitions
e.ls$n <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", 
                      partitions = parts[7], overlap = TRUE)


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
