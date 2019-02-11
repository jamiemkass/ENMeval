# set.seed(48)
# bv <- spocc::occ('Bradypus variegatus', 'gbif', limit=300, has_coords=TRUE)
# occs <- as.data.frame(bv$gbif$data$Bradypus_variegatus[,2:3])
# occs <- occs[!duplicated(occs),]
# envs <- raster::stack(list.files(path=paste(system.file(package='dismo'), '/ex', sep=''), pattern='grd', full.names=TRUE))
# which(rowSums(is.na(raster::extract(envs, occs))) > 0)
# bg <- dismo::randomPoints(envs, 1000)
# tune.args <- list(rm = 1:4, fc = c("L", "LQ"))
# partitions <- "block"
# categoricals <- "biome"
# 
# # SWD
# colnames(bg) <- c("longitude", "latitude")
# occs.vals <- raster::extract(envs, occs)
# bg.vals <- raster::extract(envs, bg)
# 
# # divide all grid cells in study extent into same partition groups
# # as the real occurrence data
# # envs.xy <- rasterToPoints(envs[[1]], spatial = TRUE)
# # envs.folds <- ENMeval::get.block(occ=occs, bg.coords=envs.xy@coords)$bg.grp
# 
# 
# ## old ENMeval
# e <- ENMevaluate(occs, envs, bg, alg = "maxnet", fc = c("L", "LQ"), RMvalues = 1:4, categoricals = "biome", method = "block")
# 
# ## regular run
# e <- ENMevaluate(occs, envs, bg, mod.fun = maxnet::maxnet, tune.args = tune.args, categoricals = "biome", partitions = "block")
# ## run with SWD
# e <- ENMevaluate(occs, bg = bg, occs.vals = occs.vals, bg.vals = bg.vals, mod.fun = maxnet::maxnet, tune.args = tune.args, categoricals = "biome", partitions = "block")
# ## run with independent testing data
# e <- ENMevaluate(occs[1:250,], envs, bg, mod.fun = maxnet::maxnet, tune.args = tune.args, categoricals = "biome", partitions = "independent", occs.ind = occs[251:nrow(occs),])
# ## run with just AIC
# e <- ENMevaluate(occs, envs, bg, mod.fun = maxnet::maxnet, tune.args = tune.args, categoricals = "biome", partitions = "none")
# ## run with maxent.jar
# e <- ENMevaluate(occs, envs, bg, mod.fun = dismo::maxent, tune.args = tune.args, categoricals = "biome", partitions = "block")
# ## run with BRT
# tune.args <- list(tree.complexity = 1, learning.rate = 0.1, bag.fraction = 0.5)
# e <- ENMevaluate(occs, envs, bg, mod.fun = dismo::gbm.step, tune.args = tune.args, categoricals = "biome", partitions = "block")
# ## run with BIOCLIM
# e <- ENMevaluate(occs, envs, bg = bg, mod.fun = dismo::bioclim, tune.args = NULL, categoricals = "biome", partitions = "block")
