library(dplyr)
options(warn=-1)

# read in data
set.seed(48)
occs <- read.csv(file.path(system.file(package="dismo"), "/ex/bradypus.csv"))[,2:3]
envs <- raster::stack(list.files(path=paste(system.file(package='dismo'), '/ex', sep=''), 
                                 pattern='grd', full.names=TRUE))
# remove categorical variable for BIOCLIM to work
envs <- envs[[-9]]
occs.z <- cbind(occs, raster::extract(envs, occs))
bg <- as.data.frame(dismo::randomPoints(envs, 1000))
names(bg) <- names(occs)
bg.z <- cbind(bg, raster::extract(envs, bg))

algorithm <- "bioclim"
tune.args <- list(tails = c("low", "high", "both"))
mset <- lapply(tune.args, function(x) x[1])
no.iter <- 5

# block partitions
context(paste("Testing ENMevaluate for", algorithm, "with block partitions..."))
e <- ENMevaluate(occs, envs, bg, tune.args = tune.args, partitions = "block", algorithm = algorithm, overlap = TRUE, quiet = TRUE)
test_ENMevaluation(e, algorithm, "block", tune.args, 4, 4)

context(paste("Testing evalplot.stats for", algorithm, "with block partitions..."))
test_evalplot.stats(e)
grps <- get.block(occs, bg)
context(paste("Testing evalplot.envSim.hist for", algorithm, "with block partitions..."))
test_evalplot.envSim.hist(e, occs.z, bg.z, grps$occs.grp, grps$bg.grp)
context(paste("Testing evalplot.envSim.map for", algorithm, "with block partitions..."))
test_evalplot.envSim.map(e, envs, occs.z, bg.z, grps$occs.grp, grps$bg.grp)

context(paste("Testing ENMnulls for", algorithm, "with block partitions..."))
ns <- ENMnulls(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
test_ENMnulls(e, ns, no.iter, algorithm, "block", mset, 4, 4)

context(paste("Testing evalplot.nulls for", algorithm, "with block partitions..."))
test_evalplot.nulls(ns)


# checkerboard1 partitions
context(paste("Testing ENMevaluate for", algorithm, "with checkerboard1 partitions..."))
e <- ENMevaluate(occs, envs, bg, tune.args = tune.args, partitions = "checkerboard1", algorithm = algorithm, overlap = TRUE, quiet = TRUE)
test_ENMevaluation(e, algorithm, "checkerboard1", tune.args, 2, 2)

context(paste("Testing evalplot.stats for", algorithm, "with checkerboard1 partitions..."))
test_evalplot.stats(e)
grps <- get.checkerboard1(occs, envs, bg, aggregation.factor = 2)
context(paste("Testing evalplot.envSim.hist for", algorithm, "with checkerboard1 partitions..."))
test_evalplot.envSim.hist(e, occs.z, bg.z, grps$occs.grp, grps$bg.grp)
context(paste("Testing evalplot.envSim.map for", algorithm, "with checkerboard1 partitions..."))
test_evalplot.envSim.map(e, envs, occs.z, bg.z, grps$occs.grp, grps$bg.grp)

context(paste("Testing ENMnulls for", algorithm, "with checkerboard1 partitions..."))
ns <- ENMnulls(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
test_ENMnulls(e, ns, no.iter, algorithm, "checkerboard1", mset, 2, 2)

context(paste("Testing ENMnulls plotting function for", algorithm, "with checkerboard1 partitions..."))
test_evalplot.nulls(ns)

# checkerboard2 partitions
context(paste("Testing ENMevaluate for", algorithm, "with checkerboard2 partitions..."))
e <- ENMevaluate(occs, envs, bg, tune.args = tune.args, partitions = "checkerboard2", algorithm = algorithm, overlap = TRUE, quiet = TRUE)
test_ENMevaluation(e, algorithm, "checkerboard2", tune.args, 4, 4)

context(paste("Testing evalplot.stats for", algorithm, "with checkerboard2 partitions..."))
test_evalplot.stats(e)
grps <- get.checkerboard2(occs, envs, bg, aggregation.factor = 2)
context(paste("Testing evalplot.envSim.hist for", algorithm, "with checkerboard2 partitions..."))
test_evalplot.envSim.hist(e, occs.z, bg.z, grps$occs.grp, grps$bg.grp)
context(paste("Testing evalplot.envSim.map for", algorithm, "with checkerboard2 partitions..."))
test_evalplot.envSim.map(e, envs, occs.z, bg.z, grps$occs.grp, grps$bg.grp)

context(paste("Testing ENMnulls for", algorithm, "with checkerboard2 partitions..."))
ns <- ENMnulls(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
test_ENMnulls(e, ns, no.iter, algorithm, "checkerboard2", mset, 4, 4)

context(paste("Testing ENMnulls plotting function for", algorithm, "with checkerboard2 partitions..."))
test_evalplot.nulls(ns)

# random k-fold partitions
context(paste("Testing ENMevaluate for", algorithm, "with random 5-fold partitions..."))
e <- ENMevaluate(occs, envs, bg, tune.args = tune.args, partitions = "randomkfold", algorithm = algorithm, overlap = TRUE, quiet = TRUE)
test_ENMevaluation(e, algorithm, "randomkfold", tune.args, 5, 1)

context(paste("Testing evalplot.stats for", algorithm, "with random 5-fold partitions..."))
test_evalplot.stats(e)
grps <- get.randomkfold(occs, bg, kfolds = 5)
context(paste("Testing evalplot.envSim.hist for", algorithm, "with random 5-fold partitions..."))
test_evalplot.envSim.hist(e, occs.z, bg.z, grps$occs.grp, grps$bg.grp, bg.sel = 0)
context(paste("Testing evalplot.envSim.map for", algorithm, "with random 5-fold partitions..."))
test_evalplot.envSim.map(e, envs, occs.z, bg.z, grps$occs.grp, grps$bg.grp, bg.sel = 0)

context(paste("Testing ENMnulls for", algorithm, "with random 5-fold partitions..."))
ns <- ENMnulls(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
test_ENMnulls(e, ns, no.iter, algorithm, "randomkfold", mset, 5, 1)

context(paste("Testing ENMnulls plotting function for", algorithm, "with random 5-fold partitions..."))
test_evalplot.nulls(ns)

# jackknife partitions
context(paste("Testing ENMevaluate for", algorithm, "with jackknife partitions..."))
e <- ENMevaluate(occs[1:10,], envs, bg, tune.args = tune.args, partitions = "jackknife", algorithm = algorithm, overlap = TRUE, quiet = TRUE)
test_ENMevaluation(e, algorithm, "jackknife", tune.args, nrow(e@occs), 1)

context(paste("Testing evalplot.stats for", algorithm, "with testing partition..."))
test_evalplot.stats(e)

context(paste("Testing ENMnulls for", algorithm, "with jackknife partitions..."))
ns <- ENMnulls(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
test_ENMnulls(e, ns, no.iter, algorithm, "jackknife", mset, nrow(e@occs), 1)

context(paste("Testing ENMnulls plotting function for", algorithm, "with jackknife partitions..."))
test_evalplot.nulls(ns)

# testing partition
context(paste("Testing ENMevaluate for", algorithm, "with testing partition..."))
e <- ENMevaluate(occs[1:100,], envs, bg, tune.args = tune.args, partitions = "testing", algorithm = algorithm, occs.testing = occs[101:nrow(occs),], overlap = TRUE, quiet = TRUE)
test_ENMevaluation(e, algorithm, "testing", tune.args, 1, 1)

context(paste("Testing evalplot.stats for", algorithm, "with testing partition..."))
test_evalplot.stats(e)
grps <- list(occs.grp = rep(0, nrow(occs)), bg.grp = rep(0, nrow(bg)))
context(paste("Testing evalplot.envSim.hist for", algorithm, "with testing partition..."))
test_evalplot.envSim.hist(e, occs.z, bg.z, grps$occs.grp, grps$bg.grp, bg.sel = 0, occs.testing.z = e@occs.testing)
context(paste("Testing evalplot.envSim.map for", algorithm, "with testing partition..."))
test_evalplot.envSim.map(e, envs, occs.z, bg.z, grps$occs.grp, grps$bg.grp, bg.sel = 0, occs.testing.z = e@occs.testing)

context(paste("Testing ENMnulls for", algorithm, "with testing partitions..."))
ns <- ENMnulls(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
test_ENMnulls(e, ns, no.iter, algorithm, "testing", mset, 1, 1)

context(paste("Testing ENMnulls plotting function for", algorithm, "with testing partition..."))
test_evalplot.nulls(ns)

# no partitions
context(paste("Testing ENMevaluate for", algorithm, "with no partitions..."))
e <- ENMevaluate(occs, envs, bg, tune.args = tune.args, partitions = "none", algorithm = algorithm, overlap = TRUE, quiet = TRUE)
test_ENMevaluation(e, algorithm, "none", tune.args, 1, 1)

context(paste("Testing ENMnulls for", algorithm, "with no partitions..."))
ns <- ENMnulls(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
test_ENMnulls(e, ns, no.iter, algorithm, "none", mset, 1, 1)

context(paste("Testing ENMnulls plotting function for", algorithm, "with no partitions..."))
test_evalplot.nulls(ns)

# user partitions
context(paste("Testing ENMevaluate for", algorithm, "with user partitions..."))
user.grp <- list(occs.grp = round(runif(nrow(occs), 1, 4)), bg.grp = round(runif(nrow(bg), 1, 4)))
e <- ENMevaluate(occs, envs, bg, tune.args = tune.args, partitions = "user", algorithm = algorithm, user.grp = user.grp, overlap = TRUE, quiet = TRUE)
test_ENMevaluation(e, algorithm, "user", tune.args, 4, 4)

context(paste("Testing evalplot.stats for", algorithm, "with user partitions..."))
test_evalplot.stats(e)
context(paste("Testing evalplot.envSim.hist for", algorithm, "with user partitions..."))
test_evalplot.envSim.hist(e, occs.z, bg.z, user.grp$occs.grp, user.grp$bg.grp)
context(paste("Testing evalplot.envSim.map for", algorithm, "with user partitions..."))
test_evalplot.envSim.map(e, envs, occs.z, bg.z, user.grp$occs.grp, user.grp$bg.grp)

context(paste("Testing ENMnulls for", algorithm, "with user partitions..."))
ns <- ENMnulls(e, mod.settings = mset, no.iter = no.iter, user.eval.type = "kspatial", quiet = TRUE)
test_ENMnulls(e, ns, no.iter, algorithm, "user", mset, 4, 4)

context(paste("Testing ENMnulls plotting function for", algorithm, "with user partitions..."))
test_evalplot.nulls(ns)

# no envs (SWD)
context(paste("Testing ENMevaluate for", algorithm, "with random 5-fold partitions and no raster environmental variables..."))
e <- ENMevaluate(occs.z, bg = bg.z, tune.args = tune.args, partitions = "randomkfold", algorithm = algorithm, quiet = TRUE)
test_ENMevaluation(e, algorithm, "randomkfold", tune.args, 5, 1, type = "swd")

context(paste("Testing evalplot.stats for", algorithm, "with random 5-fold partitions and no raster environmental variables..."))
test_evalplot.stats(e)
grps <- get.randomkfold(occs, bg, kfolds = 5)
context(paste("Testing evalplot.envSim.hist for", algorithm, "with random 5-fold partitions and no raster environmental variables..."))
test_evalplot.envSim.hist(e, occs.z, bg.z, grps$occs.grp, grps$bg.grp, bg.sel = 0)
context(paste("Testing evalplot.envSim.map for", algorithm, "with random 5-fold partitions and no raster environmental variables..."))
test_evalplot.envSim.map(e, envs, occs.z, bg.z, grps$occs.grp, grps$bg.grp, bg.sel = 0)

context(paste("Testing ENMnulls for", algorithm, "with random 5-fold partitions and no raster environmental variables..."))
ns <- ENMnulls(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
test_ENMnulls(e, ns, no.iter, algorithm, "randomkfold", mset, 5, 1)

context(paste("Testing ENMnulls plotting function for", algorithm, "with random 5-fold partitions and no raster environmental variables..."))
test_evalplot.nulls(ns)

# no bg
context(paste("Testing ENMevaluate for", algorithm, "with random 5-fold partitions and no input background data..."))
e <- ENMevaluate(occs, envs, tune.args = tune.args, partitions = "randomkfold", algorithm = algorithm, n.bg = 1000, overlap = TRUE, quiet = TRUE)
test_ENMevaluation(e, algorithm, "randomkfold", tune.args, 5, 1) 

context(paste("Testing evalplot.stats for", algorithm, "with random 5-fold partitions and no input background data..."))
test_evalplot.stats(e)
grps <- get.randomkfold(occs, bg, kfolds = 5)
context(paste("Testing evalplot.envSim.hist for", algorithm, "with random 5-fold partitions and no input background data..."))
test_evalplot.envSim.hist(e, occs.z, bg.z, grps$occs.grp, grps$bg.grp, bg.sel = 0)
context(paste("Testing evalplot.envSim.map for", algorithm, "with random 5-fold partitions and no input background data..."))
test_evalplot.envSim.map(e, envs, occs.z, bg.z, grps$occs.grp, grps$bg.grp, bg.sel = 0)

context(paste("Testing ENMnulls for", algorithm, "with random 5-fold partitions and no input background data..."))
ns <- ENMnulls(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
test_ENMnulls(e, ns, no.iter, algorithm, "randomkfold", mset, 5, 1)

context(paste("Testing ENMnulls plotting function for", algorithm, "with random 5-fold partitions and no input background data..."))
test_evalplot.nulls(ns)