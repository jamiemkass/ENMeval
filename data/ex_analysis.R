# set.seed(48)
# bv <- spocc::occ('Bradypus variegatus', 'gbif', limit=100, has_coords=TRUE)
# occs <- as.data.frame(bv$gbif$data$Bradypus_variegatus[,2:3])
# occs <- occs[!duplicated(occs),]
# envs <- raster::stack(list.files(path=paste(system.file(package='dismo'), '/ex', sep=''), pattern='grd', full.names=TRUE))
# # which(rowSums(is.na(raster::extract(envs, occs))) > 0)
# bg <- as.data.frame(dismo::randomPoints(envs, 10000))
# names(bg) <- names(occs)
# tune.args <- list(fc = c("L", "LQ"), rm = seq(1,2,0.5))
# # tune.args <- list(tails = c("low", "high", "both"))
# # tune.args <- list(fc = "L", rm = 1)
# partitions <- "randomkfold"
# kfolds <- 4
# categoricals <- "biome"
# skipRasters <- FALSE
# other.args <- NULL
# updateProgress <- NULL
# clamp <- TRUE
# abs.auc.diff <- TRUE
# pred.type <- "cloglog"
# user.grp = NULL
# occs.testing = NULL
# quiet = FALSE
# 
# # user groups
# user.grp <- list(occs.grp = rep(1,nrow(occs)), bg.grp = rep(0, nrow(bg)))
# 
# # testing partitions
# occs.testing <- occs[1:50,]
# occs <- occs[51:nrow(occs),]
# 
# # null models
# mod.settings=list(fc="L",rm=2)
# userStats.exp.sign=list(maxKappa = 1, maxTSS = 1)
# 
# # SWD
# 
# #
# # # divide all grid cells in study extent into same partition groups
# # # as the real occurrence data
# # # envs.xy <- rasterToPoints(envs[[1]], spatial = TRUE)
# # # envs.grp <- ENMeval::get.block(occ=occs, bg.coords=envs.xy@coords)$bg.grp
# #
# #
# # ## old ENMeval
# # e <- ENMevaluate(occs, envs, bg, alg = "maxnet", fc = c("L", "LQ"), RMvalues = 1:4, categoricals = "biome", method = "block")
# 
# ## regular run
# # e <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", partitions = "block", overlap = TRUE)
# ## regular run with user partitions
# # user.grp <- list(occs.grp = round(runif(nrow(occs), 1, 2)), bg.grp = round(runif(nrow(bg), 1, 2)))
# # e <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", partitions = "user", user.grp = user.grp, overlap = TRUE)
# ## regular run with NA in occurrences
# # occs[18,] <- c(0,0)
# # e <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", partitions = "block", overlap = TRUE, parallel = TRUE)
# ## run with SWD
# # e <- ENMevaluate(cbind(occs, raster::extract(envs, occs)), bg = cbind(bg, raster::extract(envs, bg)), mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", partitions = "block")
# ## run with testing validation data
# # e <- ENMevaluate(occs[1:50,], envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", partitions = "testing", occs.testing = occs[51:nrow(occs),])
# ## run with just AIC
# # e <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", partitions = "none")
# ## run with maxent.jar
# # e <- ENMevaluate(occs, envs, bg, mod.name = "maxent.jar", tune.args = tune.args, categoricals = "biome", partitions = "block")
# ## run with BRT
# # tune.args <- list(tree.complexity = 1:2, learning.rate = 0.1, bag.fraction = 0.5)
# # e <- ENMevaluate(occs, envs, bg, mod.name = "brt", tune.args = tune.args, categoricals = "biome", partitions = "block")
# # run with BIOCLIM
# # e <- ENMevaluate(occs, envs, bg, mod.name = "bioclim", categoricals = "biome", partitions = "block")
# # run parallel
# # e <- ENMevaluate(occs, envs, bg, mod.name = "maxnet", tune.args = tune.args, categoricals = "biome", partitions = "block", overlap = TRUE, parallel = TRUE)
# # null models
# ns <- ENMnullSims(e, mod.settings = list(fc = "L", rm = 2), no.iter = 10)
