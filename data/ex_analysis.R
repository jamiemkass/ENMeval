set.seed(48)
bv <- spocc::occ('Bradypus variegatus', 'gbif', limit=500, has_coords=TRUE)
occs <- as.data.frame(bv$gbif$data$Bradypus_variegatus[,2:3])
occs <- occs[!duplicated(occs),]
envs <- raster::stack(list.files(path=paste(system.file(package='dismo'), '/ex', sep=''), pattern='grd', full.names=TRUE))
which(rowSums(is.na(raster::extract(envs, occs))) > 0)
bg <- as.data.frame(dismo::randomPoints(envs, 1000))
names(bg) <- names(occs)
tune.args <- list(fc = c("L", "LQ"), rm = seq(2,3,0.5))
# tune.args <- list(fc = "L", rm = 1)
partitions <- "randomkfold"
kfolds <- 2
categoricals <- "biome"
skipRasters <- FALSE
other.args <- NULL
updateProgress <- NULL
doClamp <- TRUE
abs.auc.diff <- TRUE
pred.type <- "cloglog"
user.grp = NULL
occs.ind = NULL
cbi.eval = NULL
quiet = FALSE

# user groups
user.grp <- list(occ.grp = rep(1,nrow(occs)), bg.grp = rep(0, nrow(bg)))

# independent partitions
occs.ind <- occs[1:50,]
occs <- occs[51:nrow(occs),]

# null models
mod.settings=list(fc="L",rm=2)
userStats.exp.sign=list(maxKappa = 1, maxTSS = 1)

# SWD

# 
# # divide all grid cells in study extent into same partition groups
# # as the real occurrence data
# # envs.xy <- rasterToPoints(envs[[1]], spatial = TRUE)
# # envs.grp <- ENMeval::get.block(occ=occs, bg.coords=envs.xy@coords)$bg.grp
# 
# 
# ## old ENMeval
# e <- ENMevaluate(occs, envs, bg, alg = "maxnet", fc = c("L", "LQ"), RMvalues = 1:4, categoricals = "biome", method = "block")

## regular run
# e <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", partitions = "block", overlap = TRUE)
## regular run with user partitions
# user.grp <- list(occ.grp = round(runif(nrow(occs), 1, 2)), bg.grp = round(runif(nrow(bg), 1, 2)))
# e <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", partitions = "user", user.grp = user.grp, overlap = TRUE)
## regular run with NA in occurrences
# occs[18,] <- c(0,0)
# e <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", partitions = "block", overlap = TRUE, parallel = TRUE)
## run with SWD
# e <- ENMevaluate(cbind(occs, raster::extract(envs, occs)), bg = cbind(bg, raster::extract(envs, bg)), mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", partitions = "block")
## run with independent testing data
# e <- ENMevaluate(occs[1:250,], envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", partitions = "independent", occs.ind = occs[251:nrow(occs),])
## run with just AIC
# e <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", partitions = "none")
## run with maxent.jar
# e <- ENMevaluate(occs, envs, bg, mod.name = "maxent.jar", tune.args = tune.args, categoricals = "biome", partitions = "block")
## run with BRT
# tune.args <- list(tree.complexity = 1:2, learning.rate = 0.1, bag.fraction = 0.5)
# e <- ENMevaluate(occs, envs, bg, mod.name = "brt", tune.args = tune.args, categoricals = "biome", partitions = "block")
# run with BIOCLIM
# e <- ENMevaluate(occs, envs, bg, mod.name = "bioclim", categoricals = "biome", partitions = "block")
# run parallel
# e <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", partitions = "block", overlap = TRUE, parallel = TRUE)
# run with legacy parameters
# e <- ENMevaluate(occ = occs, env = envs, bg.coords = bg, algorithm = "maxnet", fc = c("L", "LQ"), RMvalues = 1:3, categoricals = "biome", method = "block", overlap = TRUE)
