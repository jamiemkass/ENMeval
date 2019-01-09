library(spocc)
library(raster)
library(dismo)

bv <- occ('Bradypus variegatus', 'gbif', limit=300, has_coords=TRUE)
occs <- as.data.frame(bv$gbif$data$Bradypus_variegatus[,2:3])
occs <- occs[!duplicated(occs),]
envs <- stack(list.files(path=paste(system.file(package='dismo'), '/ex', sep=''), pattern='grd', full.names=TRUE))
which(rowSums(is.na(extract(envs, occs))) > 0)
bg <- randomPoints(envs, 10000)
mod.settings <- list("rm" = 1:4, "fc" = c("L", "LQ"))
