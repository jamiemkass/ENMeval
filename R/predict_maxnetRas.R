#################################################
#########	PREDICT MAXNET RASTER #############
#################################################

# function to make a raster prediction from a maxnet object
predict_maxnetRas <- function(mod, env, type, clamp) {
  env.n <- nlayers(env)
  env.pts <- rasterToPoints(env)
  mxnet.p <- predict(mod, env.pts, type=type, clamp=clamp)
  env.pts <- cbind(env.pts, as.numeric(mxnet.p))
  mxnet.p <- rasterFromXYZ(env.pts[,c(1, 2, env.n+3)], res=res(env))
  return(mxnet.p)
}
