# bv <- spocc::occ('Bradypus variegatus', 'gbif', limit=300, has_coords=TRUE)
# occs <- as.data.frame(bv$gbif$data$Bradypus_variegatus[,2:3])
# occs <- occs[!duplicated(occs),]
# envs <- raster::stack(list.files(path=paste(system.file(package='dismo'), '/ex', sep=''), pattern='grd', full.names=TRUE))
# which(rowSums(is.na(raster::extract(envs, occs))) > 0)
# bg <- dismo::randomPoints(envs, 1000)
# mod.settings <- list("rm" = 1:4, "fc" = c("L", "LQ"))
# partitions <- "block"
# categoricals <- "biome"
# 
# # divide all grid cells in study extent into same partition groups
# # as the real occurrence data
# envs.xy <- rasterToPoints(envs[[1]], spatial = TRUE)
# envs.folds <- ENMeval::get.block(occ=occs, bg.coords=envs.xy@coords)$bg.grp
# 
