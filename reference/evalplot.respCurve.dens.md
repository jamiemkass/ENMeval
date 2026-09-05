# Response curve and density plots for one variable

A wrapper function to plot response curves and density plot.

## Usage

``` r
evalplot.respCurve.dens(
  mod,
  data,
  envs = NULL,
  var,
  fun = mean,
  type = c(1, 2),
  exp.curve = 0.025,
  nr.curve = 100,
  clamp.tails = TRUE,
  bw.envs = 10
)
```

## Arguments

- mod:

  A maxent.jar or maxnet model object.

- data:

  Data frame of training data (occurrences + background).

- envs:

  Raster data (SpatRaster) of environmental variables for model
  projection. If \`NULL\` (default), only the training-data response
  curve and density are plotted, with no transfer-environment
  comparison.

- var:

  A character string specifying the variable name for the response
  curve.

- fun:

  A function to compute constant values for other variables (default is
  \`median\`).

- type:

  Number (1 or 2) to specify type of response curve to plot. See details
  for explanation.

- exp.curve:

  Numeric value indicating the range expansion for plotting (default is
  0.025).

- nr.curve:

  Integer specifying the number of points for the response curve
  (default is 100).

- clamp.tails:

  Logical; if \`TRUE\`, clamping tails in plot (default is \`TRUE\`).

- bw.envs:

  The smoothing bandwidth to be used in the environmental variables

## Value

A combined patchwork plot of all response curves with a shared y-axis
label.

## References

Pinilla-Buitrago, G.E., Kass, J.M., & Anderson, R.P. (2026).
Extrapolation strategy matters when transferring ecological niche
models: new visualization tools for informed decisions. Ecography,
e08590. https://doi.org/10.1002/ecog.08590

## Author

Gonzalo E. Pinilla-Buitrago

## Examples

``` r
if (FALSE) { # \dontrun{
occs <- read.csv(file.path(system.file(package="predicts"), "/ex/bradypus.csv"))[,2:3]
envs <- rast(list.files(path=paste(system.file(package="predicts"), "/ex", sep=""),
                        pattern="tif$", full.names=TRUE))
# No biome
envs <- envs[[!(names(envs) %in% "biome")]]
occs.z <- cbind(occs, terra::extract(envs, occs, ID = FALSE))
bg <- as.data.frame(predicts::backgroundSample(envs, n = 10000))
names(bg) <- names(occs)
bg.z <- cbind(bg, terra::extract(envs, bg, ID = FALSE))
os <- list(abs.auc.diff = FALSE, pred.type = "cloglog", validation.bg = "partition")
ps <- list(orientation = "lat_lat")
# Transfer envs
tr_envs <- envs * 1.5
# Plot
mod <- e@models[[1]]
# Define data as combined training values with coordinates removed
data <- rbind(e@occs, e@bg)[,3:11]
evalplot.respCurve.dens(mod, data, envs = tr_envs, var = "bio1", fun = median)
# Plot training data only, without a transfer environment
evalplot.respCurve.dens(mod, data, var = "bio1", fun = median)
} # }
```
