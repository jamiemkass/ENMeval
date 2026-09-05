# Plot Response Curve for Maxent Models

This function plots a response curve for a given environmental variable
based on a maxent.jar or maxnet model. It allows plotting clamping on or
off and supports multiple variables via a wrapper that combines plots
using the patchwork package.

## Usage

``` r
evalplot.respCurve(
  mod,
  data,
  envs = NULL,
  var,
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
  curve is plotted, with no transfer-environment comparison.

- var:

  A character string specifying the variable name for the response
  curve.

- fun:

  If maxent.jar a function to compute constant values for other
  variables (default is \`mean\`). Maxnet models always use mean.

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

A ggplot object of the response curve.

## Details

The type 1 option (default) sets the focal variable to values along a
range "r" from its minimum to maximum (buffered by exp.curve) while
setting all other variables to static values defined by fun (which
defaults to their means), then makes a model prediction for this table.
The type 2 option sets the focal variable to one static value along r
while keeping all other variables at their original values, makes a
model predicton for this table, then repeats this process for all values
along the range, resulting in 100 model predictions for an r of length
100. The final curve for type 2 plots the means of these prediction
tables.

The original maxent.jar software and the dismo package implemented type
1 response curves, but the pdp package and the predicts package
implement type 2, so the user can choose which to visualize in order to
directly compare to one of these outputs.

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
e <- ENMevaluate(occs, envs, bg,
                 tune.args = list(fc = "LQ", rm = 1),
                 partitions = "block", other.settings = os, partition.settings = ps,
                 algorithm = "maxnet", overlap = TRUE)
# Transfer envs
tr_envs <- envs * 1.5
# Plot
# Plot with clamp tails
mod <- e@models[[1]]
# Define data as combined training values with coordinates removed
data <- rbind(e@occs, e@bg)[,3:11]
# Plot
evalplot.respCurve(mod, data, envs = tr_envs, var = "bio1")
# Without tails
evalplot.respCurve(mod, data, envs = tr_envs, var = "bio1", clamp.tails = FALSE)
} # }
```
