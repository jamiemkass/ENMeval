# set to FALSE to run a comprehensive set of tests
# when TRUE, only some essential tests are run to avoid lagging when 
# submitting to CRAN
skip_tests_for_cran <- TRUE

library(dplyr)
options(warn=-1)

# read in data
set.seed(48)
occs <- read.csv(file.path(system.file(package="predicts"), 
                           "/ex/bradypus.csv"))[,2:3]
envs.orig <- terra::rast(list.files(path = paste(system.file(package='predicts'), 
                                                 '/ex', sep=''), 
                                    pattern='tif$', full.names=TRUE))
bg <- predicts::backgroundSample(envs.orig, n = 1000) |> as.data.frame()
names(bg) <- names(occs)

# define SWD tables for testing
occs.z <- cbind(occs, terra::extract(envs.orig, occs, ID = FALSE))
occs.z$biome <- factor(occs.z$biome)
bg.z <- cbind(bg, terra::extract(envs.orig, bg, ID = FALSE))
bg.z$biome <- factor(bg.z$biome)

alg <- "bioclim"
no.iter <- 5

# special settings for BIOCLIM
envs <- envs.orig[[-10]]
occs.z$biome <- NULL
bg.z$biome <- NULL
cats1 <- NULL

# define tune args
tune.args <- list(tails = c("low", "high", "both"))
mset <- lapply(tune.args, function(x) x[1])

if(skip_tests_for_cran == FALSE) {
  # block partitions
  context(paste("Testing ENMevaluate for", alg, "with block partitions..."))
  e <- ENMevaluate(occs, envs, bg, tune.args = tune.args, partitions = "block", algorithm = alg, categoricals = cats1, overlap = TRUE, quiet = TRUE)
  test_ENMevaluation(e, alg, "block", tune.args, 4, 4)
  
  context(paste("Testing evalplot.stats for", alg, "with block partitions..."))
  test_evalplot.stats(e)
  context(paste("Testing evalplot.envSim.hist for", alg, "with block partitions..."))
  test_evalplot.envSim.hist(e, e@occs, e@bg, e@occs.grp, e@bg.grp)
  context(paste("Testing evalplot.envSim.map for", alg, "with block partitions..."))
  test_evalplot.envSim.map(e, envs, e@occs, e@bg, e@occs.grp, e@bg.grp)
  
  context(paste("Testing ENMnulls for", alg, "with block partitions..."))
  ns <- ENMnulls(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
  test_ENMnulls(e, ns, no.iter, alg, "block", mset, 4, 4)
  
  context(paste("Testing evalplot.nulls for", alg, "with block partitions..."))
  test_evalplot.nulls(ns)
}

# block partitions with parallel processing
if(skip_tests_for_cran == FALSE) {
  context(paste("Testing ENMevaluate for", alg, "with block partitions in parallel..."))
  e <- ENMevaluate(occs, envs, bg, tune.args = tune.args, partitions = "block",
                   algorithm = alg, categoricals = cats1, overlap = TRUE, quiet = TRUE,
                   parallel = TRUE)
  test_ENMevaluation(e, alg, "block", tune.args, 4, 4)
  
  context(paste("Testing evalplot.stats for", alg, "with block partitions in parallel..."))
  test_evalplot.stats(e)
  context(paste("Testing evalplot.envSim.hist for", alg, "with block partitions in parallel..."))
  test_evalplot.envSim.hist(e, e@occs, e@bg, e@occs.grp, e@bg.grp)
  context(paste("Testing evalplot.envSim.map for", alg, "with block partitions in parallel..."))
  test_evalplot.envSim.map(e, envs, e@occs, e@bg, e@occs.grp, e@bg.grp)
  
  context(paste("Testing ENMnulls for", alg, "with block partitions in parallel..."))
  ns <- ENMnulls(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE, parallel = TRUE)
  test_ENMnulls(e, ns, no.iter, alg, "block", mset, 4, 4)
  
  context(paste("Testing evalplot.nulls for", alg, "with block partitions in parallel..."))
  test_evalplot.nulls(ns)
}

# checkerboard1 partitions
if(skip_tests_for_cran == FALSE) {
  context(paste("Testing ENMevaluate for", alg, "with checkerboard1 partitions..."))
  e <- ENMevaluate(occs, envs, bg, tune.args = tune.args, partitions = "checkerboard", algorithm = alg, categoricals = cats1, overlap = TRUE, quiet = TRUE)
  test_ENMevaluation(e, alg, "checkerboard", tune.args, 2, 2)
  
  context(paste("Testing evalplot.stats for", alg, "with checkerboard1 partitions..."))
  test_evalplot.stats(e)
  context(paste("Testing evalplot.envSim.hist for", alg, "with checkerboard1 partitions..."))
  test_evalplot.envSim.hist(e, e@occs, e@bg, e@occs.grp, e@bg.grp)
  context(paste("Testing evalplot.envSim.map for", alg, "with checkerboard1 partitions..."))
  test_evalplot.envSim.map(e, envs, e@occs, e@bg, e@occs.grp, e@bg.grp)
  
  context(paste("Testing ENMnulls for", alg, "with checkerboard1 partitions..."))
  ns <- ENMnulls(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
  test_ENMnulls(e, ns, no.iter, alg, "checkerboard", mset, 2, 2)
  
  context(paste("Testing ENMnulls plotting function for", alg, "with checkerboard1 partitions..."))
  test_evalplot.nulls(ns)  
}


# checkerboard2 partitions
if(skip_tests_for_cran == FALSE) {
  context(paste("Testing ENMevaluate for", alg, "with checkerboard2 partitions..."))
  e <- ENMevaluate(occs, envs, bg, tune.args = tune.args, partitions = "checkerboard", partition.settings = list(aggregation.factor = c(2,2)), algorithm = alg, categoricals = cats1, overlap = TRUE, quiet = TRUE)
  test_ENMevaluation(e, alg, "checkerboard", tune.args, 4, 4)
  
  context(paste("Testing evalplot.stats for", alg, "with checkerboard2 partitions..."))
  test_evalplot.stats(e)
  context(paste("Testing evalplot.envSim.hist for", alg, "with checkerboard2 partitions..."))
  test_evalplot.envSim.hist(e, e@occs, e@bg, e@occs.grp, e@bg.grp)
  context(paste("Testing evalplot.envSim.map for", alg, "with checkerboard2 partitions..."))
  test_evalplot.envSim.map(e, envs, e@occs, e@bg, e@occs.grp, e@bg.grp)
  
  context(paste("Testing ENMnulls for", alg, "with checkerboard2 partitions..."))
  ns <- ENMnulls(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
  test_ENMnulls(e, ns, no.iter, alg, "checkerboard", mset, 4, 4)
  
  context(paste("Testing ENMnulls plotting function for", alg, "with checkerboard2 partitions..."))
  test_evalplot.nulls(ns)
}

# random k-fold partitions
if(skip_tests_for_cran == FALSE) {
  context(paste("Testing ENMevaluate for", alg, "with random 5-fold partitions..."))
  e <- ENMevaluate(occs, envs, bg, tune.args = tune.args, partitions = "randomkfold", algorithm = alg, categoricals = cats1, overlap = TRUE, quiet = TRUE)
  test_ENMevaluation(e, alg, "randomkfold", tune.args, 5, 1)
  
  context(paste("Testing evalplot.stats for", alg, "with random 5-fold partitions..."))
  test_evalplot.stats(e)
  context(paste("Testing evalplot.envSim.hist for", alg, "with random 5-fold partitions..."))
  test_evalplot.envSim.hist(e, e@occs, e@bg, e@occs.grp, e@bg.grp, bg.sel = 0)
  context(paste("Testing evalplot.envSim.map for", alg, "with random 5-fold partitions..."))
  test_evalplot.envSim.map(e, envs, e@occs, e@bg, e@occs.grp, e@bg.grp, bg.sel = 0)
  
  context(paste("Testing ENMnulls for", alg, "with random 5-fold partitions..."))
  ns <- ENMnulls(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
  test_ENMnulls(e, ns, no.iter, alg, "randomkfold", mset, 5, 1)
  
  context(paste("Testing ENMnulls plotting function for", alg, "with random 5-fold partitions..."))
  test_evalplot.nulls(ns)
}

# jackknife partitions
if(skip_tests_for_cran == FALSE) {
  context(paste("Testing ENMevaluate for", alg, "with jackknife partitions..."))
  e <- ENMevaluate(occs[1:10,], envs, bg, tune.args = tune.args, partitions = "jackknife", algorithm = alg, overlap = TRUE, quiet = TRUE)
  test_ENMevaluation(e, alg, "jackknife", tune.args, nrow(e@occs), 1)
  
  context(paste("Testing evalplot.stats for", alg, "with testing partition..."))
  test_evalplot.stats(e)
  
  context(paste("Testing ENMnulls for", alg, "with jackknife partitions..."))
  ns <- ENMnulls(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
  test_ENMnulls(e, ns, no.iter, alg, "jackknife", mset, nrow(e@occs), 1)
  
  context(paste("Testing ENMnulls plotting function for", alg, "with jackknife partitions..."))
  test_evalplot.nulls(ns)
}

# testing partition
if(skip_tests_for_cran == FALSE) {
  context(paste("Testing ENMevaluate for", alg, "with testing partition..."))
  e <- ENMevaluate(occs[1:100,], envs, bg, tune.args = tune.args, partitions = "testing", algorithm = alg, categoricals = cats1,  occs.testing = occs[101:nrow(occs),], overlap = TRUE, quiet = TRUE)
  test_ENMevaluation(e, alg, "testing", tune.args, 1, 1)
  
  context(paste("Testing evalplot.stats for", alg, "with testing partition..."))
  test_evalplot.stats(e)
  context(paste("Testing evalplot.envSim.hist for", alg, "with testing partition..."))
  test_evalplot.envSim.hist(e, e@occs, e@bg, e@occs.grp, e@bg.grp, bg.sel = 0, occs.testing.z = e@occs.testing)
  context(paste("Testing evalplot.envSim.map for", alg, "with testing partition..."))
  test_evalplot.envSim.map(e, envs, e@occs, e@bg, e@occs.grp, e@bg.grp, bg.sel = 0, occs.testing.z = e@occs.testing)
  
  context(paste("Testing ENMnulls for", alg, "with testing partitions..."))
  ns <- ENMnulls(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
  test_ENMnulls(e, ns, no.iter, alg, "testing", mset, 1, 1)
  
  context(paste("Testing ENMnulls plotting function for", alg, "with testing partition..."))
  test_evalplot.nulls(ns)
}

# no partitions
if(skip_tests_for_cran == FALSE) {
  context(paste("Testing ENMevaluate for", alg, "with no partitions..."))
  e <- ENMevaluate(occs, envs, bg, tune.args = tune.args, partitions = "none", algorithm = alg, categoricals = cats1, overlap = TRUE, quiet = TRUE)
  test_ENMevaluation(e, alg, "none", tune.args, 1, 1)
  
  context(paste("Testing ENMnulls for", alg, "with no partitions..."))
  ns <- ENMnulls(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
  test_ENMnulls(e, ns, no.iter, alg, "none", mset, 1, 1)
  
  context(paste("Testing ENMnulls plotting function for", alg, "with no partitions..."))
  test_evalplot.nulls(ns)
}

# user partitions
if(skip_tests_for_cran == FALSE) {
  context(paste("Testing ENMevaluate for", alg, "with user partitions..."))
  user.grp <- list(occs.grp = round(runif(nrow(occs), 1, 4)), bg.grp = round(runif(nrow(bg), 1, 4)))
  e <- ENMevaluate(occs, envs, bg, tune.args = tune.args, partitions = "user", algorithm = alg, categoricals = cats1, user.grp = user.grp, overlap = TRUE, quiet = TRUE)
  test_ENMevaluation(e, alg, "user", tune.args, 4, 4)
  
  context(paste("Testing evalplot.stats for", alg, "with user partitions..."))
  test_evalplot.stats(e)
  context(paste("Testing evalplot.envSim.hist for", alg, "with user partitions..."))
  test_evalplot.envSim.hist(e, e@occs, e@bg, e@occs.grp, e@bg.grp)
  context(paste("Testing evalplot.envSim.map for", alg, "with user partitions..."))
  test_evalplot.envSim.map(e, envs, e@occs, e@bg, e@occs.grp, e@bg.grp)
  
  context(paste("Testing ENMnulls for", alg, "with user partitions..."))
  ns <- ENMnulls(e, mod.settings = mset, no.iter = no.iter, user.eval.type = "kspatial", quiet = TRUE)
  test_ENMnulls(e, ns, no.iter, alg, "user", mset, 4, 4)
  
  context(paste("Testing ENMnulls plotting function for", alg, "with user partitions..."))
  test_evalplot.nulls(ns)
}

# no envs (SWD)
if(skip_tests_for_cran == FALSE) {
  context(paste("Testing ENMevaluate for", alg, "with random 5-fold partitions and no raster environmental variables..."))
  
  e <- ENMevaluate(occs.z, bg = bg.z, tune.args = tune.args, partitions = "randomkfold", algorithm = alg, categoricals = cats1, quiet = TRUE)
  test_ENMevaluation(e, alg, "randomkfold", tune.args, 5, 1, type = "swd")
  
  context(paste("Testing evalplot.stats for", alg, "with random 5-fold partitions and no raster environmental variables..."))
  test_evalplot.stats(e)
  context(paste("Testing evalplot.envSim.hist for", alg, "with random 5-fold partitions and no raster environmental variables..."))
  test_evalplot.envSim.hist(e, e@occs, e@bg, e@occs.grp, e@bg.grp, bg.sel = 0)
  context(paste("Testing evalplot.envSim.map for", alg, "with random 5-fold partitions and no raster environmental variables..."))
  test_evalplot.envSim.map(e, envs, e@occs, e@bg, e@occs.grp, e@bg.grp, bg.sel = 0)
  
  context(paste("Testing ENMnulls for", alg, "with random 5-fold partitions and no raster environmental variables..."))
  ns <- ENMnulls(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
  test_ENMnulls(e, ns, no.iter, alg, "randomkfold", mset, 5, 1)
  
  context(paste("Testing ENMnulls plotting function for", alg, "with random 5-fold partitions and no raster environmental variables..."))
  test_evalplot.nulls(ns)
}

# no bg
if(skip_tests_for_cran == FALSE) {
  context(paste("Testing ENMevaluate for", alg, "with random 5-fold partitions and no input background data..."))
  e <- ENMevaluate(occs, envs, tune.args = tune.args, partitions = "randomkfold", algorithm = alg, n.bg = 1000, categoricals = cats1, overlap = TRUE, quiet = TRUE)
  test_ENMevaluation(e, alg, "randomkfold", tune.args, 5, 1) 
  
  context(paste("Testing evalplot.stats for", alg, "with random 5-fold partitions and no input background data..."))
  test_evalplot.stats(e)
  context(paste("Testing evalplot.envSim.hist for", alg, "with random 5-fold partitions and no input background data..."))
  test_evalplot.envSim.hist(e, e@occs, e@bg, e@occs.grp, e@bg.grp, bg.sel = 0)
  context(paste("Testing evalplot.envSim.map for", alg, "with random 5-fold partitions and no input background data..."))
  test_evalplot.envSim.map(e, envs, e@occs, e@bg, e@occs.grp, e@bg.grp, bg.sel = 0)
  
  context(paste("Testing ENMnulls for", alg, "with random 5-fold partitions and no input background data..."))
  ns <- ENMnulls(e, mod.settings = mset, no.iter = no.iter, quiet = TRUE)
  test_ENMnulls(e, ns, no.iter, alg, "randomkfold", mset, 5, 1)
  
  context(paste("Testing ENMnulls plotting function for", alg, "with random 5-fold partitions and no input background data..."))
  test_evalplot.nulls(ns)
}

if(skip_tests_for_cran == FALSE) {
  # clamping
  context(paste("Testing clamping function for", alg, "with..."))
  test_clamp(envs, e@occs, e@bg, categoricals = cats1)
  context(paste("Testing clamping function for", alg, "with two categorical variables..."))
  if(skip_tests_for_cran == FALSE & alg != "bioclim") test_clamp(envs.2cat, e.2cat@occs, e.2cat@bg, categoricals = c("biome.1", "biome.2"))
}
