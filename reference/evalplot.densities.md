# Density plots for All Variables with Shared Y-Axis

A wrapper function to plot response curves for all contributing
variables and combine them using patchwork. The plots share a common
y-axis label.

## Usage

``` r
evalplot.densities(data, envs = NULL, vars = NULL, bw.envs = 10)
```

## Arguments

- data:

  Data frame of training data (occurrences + background).

- envs:

  Raster data (SpatRaster) of environmental variables for model
  projection. If \`NULL\` (default), only the training-data densities
  are plotted, with no transfer-environment comparison.

- vars:

  Vector specifying the variable names for the response curve. Default
  is all variables (from \`envs\` if supplied, otherwise from \`data\`).

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
e <- ENMevaluate(occs, envs, bg,
                 tune.args = list(fc = "LQ", rm = 1),
                 partitions = "block", other.settings = os, partition.settings = ps,
                 algorithm = "maxnet", overlap = TRUE)
# Transfer envs
tr_envs <- envs * 1.5
# Define data as combined training values with coordinates removed
data <- rbind(e@occs, e@bg)[,3:11]
# Plot
evalplot.densities(data, envs = tr_envs)
# Plot training data only, without a transfer environment
evalplot.densities(data)
} # }
```
