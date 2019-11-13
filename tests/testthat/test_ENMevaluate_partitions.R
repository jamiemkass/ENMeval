context("ENMevaluate partitions")

# read in data
set.seed(48)
occs <- readRDS("data/bvariegatus.rds")
envs <- raster::stack(list.files(path=paste(system.file(package='dismo'), '/ex', sep=''), 
                                 pattern='grd', full.names=TRUE))
bg <- as.data.frame(dismo::randomPoints(envs, 1000))
names(bg) <- names(occs)
tune.args <- list(fc = c("L","LQ","H"), rm = 1:5)
tune.args <- list(fc = "L", rm = 2:3)
# random sample of occs for runs that use subsets
i <- sample(1:nrow(occs))

# maxnet run with tuning parameters, categorical variable, block partitions, and niche overlap
e.block <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", 
                  partitions = "block", overlap = TRUE)
# checkerboard1 partitions
e.cb1 <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", 
                  partitions = "checkerboard1", overlap = TRUE)
# checkerboard2 partitions
e.cb2 <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", 
                  partitions = "checkerboard2", overlap = TRUE)
# random k-fold partitions
e <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", 
                  partitions = "randomkfold", kfolds = 4, overlap = TRUE)
# jackknife partitions
e5 <- ENMevaluate(occs[i[1:10],], envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", 
                  partitions = "jackknife", overlap = TRUE)
# independent partition
e6 <- ENMevaluate(occs[i[1:450],], envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", 
                  partitions = "independent", occs.ind = occs[i[451:nrow(occs)],], overlap = TRUE)
# no partitions
e7 <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", 
                  partitions = "none", overlap = TRUE)


# test_that("")

