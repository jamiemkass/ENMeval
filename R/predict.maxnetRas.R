#################################################
#########	PREDICT MAXNET RASTER #############
#################################################

# function to make a raster prediction from a maxnet object
predict.maxnetRas <- function(mod, env, type, clamp) {
  env.n <- nlayers(env)
  env.pts <- rasterToPoints(env)
  origNrow <- nrow(env.pts)
  env.pts <- na.omit(env.pts)
  naOmitNrow <- nrow(env.pts)
  rowDiff <- origNrow - naOmitNrow
  if (rowDiff > 0) {
    message(paste('\n', rowDiff, "grid cells found with at least one NA value: these cells were not processed."))
  }
  mxnet.p <- predict(mod, env.pts, type=type, clamp=clamp)
  env.pts <- cbind(env.pts, as.numeric(mxnet.p))
  mxnet.p <- rasterFromXYZ(env.pts[,c(1, 2, env.n+3)], res=res(env))
  return(mxnet.p)
}
