# Predict raster for maxnet

As maxnet does not have native functionality for making prediction
rasters, this function does it. It is a wrapper for the internal
enm.maxnet@predict function, which is not really intended for use
outside the package.

## Usage

``` r
maxnet.predictRaster(mod, envs, pred.type = "cloglog", doClamp = TRUE, ...)
```

## Arguments

- mod:

  maxnet object

- envs:

  SpatRaster

- pred.type:

  character: the type of maxnet prediction to make; the default is
  "cloglog"

- doClamp:

  Boolean: whether to clamp predictions or not

- ...:

  any additional parameters
