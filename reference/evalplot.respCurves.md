# Plot Response Curves for All Variables with Shared Y-Axis

A wrapper function to plot response curves for all contributing
variables and combine them using patchwork. The plots share a common
y-axis label.

## Usage

``` r
evalplot.respCurves(
  mod,
  data,
  envs = NULL,
  fun = mean,
  type = c(1, 2),
  exp.curve = 0.025,
  nr.curve = 100,
  clamp.tails = TRUE
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
  curves are plotted, with no transfer-environment comparison.

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

## Value

A combined patchwork plot of all response curves with a shared y-axis
label.

## References

Pinilla-Buitrago, G.E., Kass, J.M., & Anderson, R.P. (2026).
Extrapolation strategy matters when transferring ecological niche
models: new visualization tools for informed decisions. Ecography,
e08590. https://doi.org/10.1002/ecog.08590

## Author

Gonzalo E. Pinilla- Buitrago

## Examples

``` r
if (FALSE) { # \dontrun{
library(ENMeval)
occs <- read.csv(file.path(system.file(package="predicts"), "/ex/bradypus.csv"))[,2:3]
envs <- terra::rast(list.files(path=paste(system.file(package="predicts"), "/ex", sep=""),
                        pattern="tif$", full.names=TRUE))
# No biome
envs <- envs[[!(names(envs) %in% "biome")]]
occs.z <- cbind(occs, terra::extract(envs, occs, ID = FALSE))
bg <- as.data.frame(predicts::backgroundSample(envs, n = 10000))
names(bg) <- names(occs)
bg.z <- cbind(bg, terra::extract(envs, bg, ID = FALSE))
os <- list(abs.auc.diff = FALSE, pred.type = "cloglog", validation.bg = "partition")
ps <- list(orientation = "lat_lat")
e <- ENMevaluate(occs, envs, bg, tune.args = list(fc = "LQ", rm = 1),
                 partitions = "block", other.settings = os,
                 partition.settings = ps, algorithm = "maxnet", overlap = TRUE)
# Transfer envs
tr_envs <- envs * 1.5
mod <- e@models[[1]]
# Define data as combined training values with coordinates removed
data <- rbind(e@occs, e@bg)[,3:11]
# Plot
evalplot.respCurves(mod, data, envs = tr_envs)
# Plot training data only, without a transfer environment
evalplot.respCurves(mod, data)
} # }
```
