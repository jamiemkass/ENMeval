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

algorithm <- "randomForest"
tune.args <- list(num.trees = 1000, mtry = 4:6)
mset <- lapply(tune.args, function(x) x[1])
no.iter <- 5

# block partitions
context(paste("Testing ENMevaluate for", algorithm, "with block partitions..."))
e <- ENMevaluate(occs, envs, bg, tune.args = tune.args, partitions = "block", algorithm = algorithm, categoricals = "biome", overlap = TRUE, quiet = TRUE)
test_ENMevaluation(e, algorithm, "block", tune.args, 4, 4)
context(paste("Testing ENMnullSims for", algorithm, "with block partitions..."))
ns <- ENMnullSims(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
test_ENMnullSims(e, ns, no.iter, algorithm, "block", mset, 4, 4)

# checkerboard1 partitions
context(paste("Testing ENMevaluate for", algorithm, "with checkerboard1 partitions..."))
e <- ENMevaluate(occs, envs, bg, tune.args = tune.args, partitions = "checkerboard1", algorithm = algorithm, categoricals = "biome", overlap = TRUE, quiet = TRUE)
test_ENMevaluation(e, algorithm, "checkerboard1", tune.args, 2, 2)
context(paste("Testing ENMnullSims for", algorithm, "with checkerboard1 partitions..."))
ns <- ENMnullSims(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
test_ENMnullSims(e, ns, no.iter, algorithm, "checkerboard1", mset, 2, 2)

# checkerboard2 partitions
context(paste("Testing ENMevaluate for", algorithm, "with checkerboard2 partitions..."))
e <- ENMevaluate(occs, envs, bg, tune.args = tune.args, partitions = "checkerboard2", algorithm = algorithm, categoricals = "biome", overlap = TRUE, quiet = TRUE)
test_ENMevaluation(e, algorithm, "checkerboard2", tune.args, 4, 4)
context(paste("Testing ENMnullSims for", algorithm, "with checkerboard2 partitions..."))
ns <- ENMnullSims(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
test_ENMnullSims(e, ns, no.iter, algorithm, "checkerboard2", mset, 4, 4)

# random k-fold partitions
context(paste("Testing ENMevaluate for", algorithm, "with random 4-fold partitions..."))
e <- ENMevaluate(occs, envs, bg, tune.args = tune.args, partitions = "randomkfold", algorithm = algorithm, categoricals = "biome", overlap = TRUE, quiet = TRUE)
test_ENMevaluation(e, algorithm, "randomkfold", tune.args, 5, 1)
context(paste("Testing ENMnullSims for", algorithm, "with random 4-fold partitions..."))
ns <- ENMnullSims(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
test_ENMnullSims(e, ns, no.iter, algorithm, "randomkfold", mset, 5, 1)

# jackknife partitions
context(paste("Testing ENMevaluate for", algorithm, "with jackknife partitions..."))
e <- ENMevaluate(occs[1:10,], envs, bg, tune.args = tune.args, partitions = "jackknife", algorithm = algorithm, overlap = TRUE, quiet = TRUE)
test_ENMevaluation(e, algorithm, "jackknife", tune.args, nrow(e@occs), 1)
context(paste("Testing ENMnullSims for", algorithm, "with jackknife partitions..."))
ns <- ENMnullSims(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
test_ENMnullSims(e, ns, no.iter, algorithm, "jackknife", mset, nrow(e@occs), 1)

# testing partition
context(paste("Testing ENMevaluate for", algorithm, "with testing partitions..."))
e <- ENMevaluate(occs[1:100,], envs, bg, tune.args = tune.args, partitions = "testing", algorithm = algorithm, categoricals = "biome",  occs.testing = occs[101:nrow(occs),], overlap = TRUE, quiet = TRUE)
test_ENMevaluation(e, algorithm, "testing", tune.args, 1, 1)
context(paste("Testing ENMnullSims for", algorithm, "with testing partitions..."))
ns <- ENMnullSims(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
test_ENMnullSims(e, ns, no.iter, algorithm, "testing", mset, 1, 1)

# no partitions
context(paste("Testing ENMevaluate for", algorithm, "with no partitions..."))
e <- ENMevaluate(occs, envs, bg, tune.args = tune.args, partitions = "none", algorithm = algorithm, categoricals = "biome", overlap = TRUE, quiet = TRUE)
test_ENMevaluation(e, algorithm, "none", tune.args, 1, 1)
context(paste("Testing ENMnullSims for", algorithm, "with no partitions..."))
ns <- ENMnullSims(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
test_ENMnullSims(e, ns, no.iter, algorithm, "none", mset, 1, 1)

# user partitions
context(paste("Testing ENMevaluate for", algorithm, "with user partitions..."))
user.grp <- list(occs.grp = round(runif(nrow(occs), 1, 4)), bg.grp = round(runif(nrow(bg), 1, 4)))
e <- ENMevaluate(occs, envs, bg, tune.args = tune.args, partitions = "user", algorithm = algorithm, categoricals = "biome", user.grp = user.grp, overlap = TRUE, quiet = TRUE)
test_ENMevaluation(e, algorithm, "user", tune.args, 4, 4)
context(paste("Testing ENMnullSims for", algorithm, "with user partitions..."))
ns <- ENMnullSims(e, mod.settings = mset, no.iter = no.iter, user.eval.type = "kspatial", quiet = TRUE)
test_ENMnullSims(e, ns, no.iter, algorithm, "user", mset, 4, 4)

# no envs (SWD)
context(paste("Testing ENMevaluate for", algorithm, "with random 4-fold partitions and no raster environmental variables..."))
e <- ENMevaluate(occs.z, bg = bg.z, tune.args = tune.args, partitions = "randomkfold", algorithm = algorithm, categoricals = "biome", quiet = TRUE)
test_ENMevaluation(e, algorithm, "randomkfold", tune.args, 5, 1, type = "swd")
context(paste("Testing ENMnullSims for", algorithm, "with random 4-fold partitions and no raster environmental variables..."))
ns <- ENMnullSims(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
test_ENMnullSims(e, ns, no.iter, algorithm, "randomkfold", mset, 5, 1)

# no bg
context(paste("Testing ENMevaluate for", algorithm, "with random 4-fold partitions and no input background data..."))
e <- ENMevaluate(occs, envs, tune.args = tune.args, partitions = "randomkfold", algorithm = algorithm, n.bg = 1000, categoricals = "biome", overlap = TRUE, quiet = TRUE)
test_ENMevaluation(e, algorithm, "randomkfold", tune.args, 5, 1) 
context(paste("Testing ENMnullSims for", algorithm, "with random 4-fold partitions and no input background data..."))
ns <- ENMnullSims(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
test_ENMnullSims(e, ns, no.iter, algorithm, "randomkfold", mset, 5, 1)



