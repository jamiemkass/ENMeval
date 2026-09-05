# Partition group plots

Plot occurrence partition groups over an environmental predictor raster.

## Usage

``` r
evalplot.grps(
  e = NULL,
  envs,
  pts = NULL,
  pts.grp = NULL,
  ref.data = "occs",
  pts.size = 1.5,
  return.tbl = FALSE
)
```

## Arguments

- e:

  ENMevaluation object

- envs:

  SpatRaster: environmental predictor variable used to build the models
  in "e"

- pts:

  matrix / data frame: coordinates for occurrence or background data

- pts.grp:

  numeric vector: partition groups corresponding to data in "pts"

- ref.data:

  character: plot occurrences ("occs") or background ("bg"), with
  default "occs"

- pts.size:

  numeric: custom point size for ggplot

- return.tbl:

  boolean: if TRUE, return the data frames used to make the ggplot
  instead of the plot itself

## Details

This function serves as a quick way to visualize occurrence or
background partitions over the extent of an environmental predictor
raster. It can be run with an existing ENMevaluation object, or
alternatively with occurrence or background coordinates and the
corresponding partitions.

## Examples

``` r
if (FALSE) { # \dontrun{
library(ENMeval)
occs <- read.csv(file.path(system.file(package="predicts"), 
                           "/ex/bradypus.csv"))[,2:3]
envs <- terra::rast(list.files(path=paste(system.file(package="predicts"), 
                                   "/ex", sep=""), pattern="tif$", full.names=TRUE))
bg <- as.data.frame(predicts::backgroundSample(envs, n = 10000))
names(bg) <- names(occs)

parts <- get.block(occs, bg, orientation = "lat_lon")

# now, plot the partition groups for occurrence and background points
evalplot.grps(envs = envs, pts = occs, pts.grp = parts$occs.grp)
evalplot.grps(envs = envs, pts = bg, pts.grp = parts$bg.grp)

# you can also plot with an ENMevaluation object
ps <- list(orientation = "lat_lon")
e <- ENMevaluate(occs, envs, bg, 
                 tune.args = list(fc = c("L","LQ"), rm = 1:3), 
                 partitions = "block", partition.settings = ps, 
                 algorithm = "maxnet", categoricals = "biome", 
                 parallel = TRUE)

evalplot.grps(e = e, envs = envs, ref.data = "occs")
evalplot.grps(e = e, envs = envs, ref.data = "bg")
} # }
```
