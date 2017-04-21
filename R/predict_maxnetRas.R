#################################################
#########	PREDICT MAXNET RASTER #############
#################################################

# function to make a raster prediction from a maxnet object
predict_maxnetRas <- function(mod, envs, type, clamp) {
  envs.n <- nlayers(envs)
  envs.pts <- rasterToPoints(envs)
  mxnet.p <- predict(mod, envs.pts, type=type, clamp=clamp)
  envs.pts <- cbind(envs.pts, as.numeric(mxnet.p))
  mxnet.p <- rasterFromXYZ(envs.pts[,c(1,2,envs.n+1)], res=res(envs))
  return(mxnet.p)
}
