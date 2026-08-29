# ENMevaluation statistics plot

Plot evaluation statistics over tuning parameter ranges to visualize
differences in performance

## Usage

``` r
evalplot.stats(
  e,
  stats,
  x.var,
  color.var,
  dodge = NULL,
  error.bars = TRUE,
  facet.labels = NULL,
  metric.levels = NULL,
  return.tbl = FALSE
)
```

## Arguments

- e:

  ENMevaluation object

- stats:

  character vector: names of statistics from results table to be
  plotted; if more than one statistic is specified, the plot will be
  faceted

- x.var:

  character: variable to be plotted on x-axis

- color.var:

  character: variable used to assign symbology colors

- dodge:

  numeric: dodge width for points and lines; this improves visibility
  when there is high overlap (optional)

- error.bars:

  boolean: if TRUE, plot error bars

- facet.labels:

  character vector: custom names for the metric facets

- metric.levels:

  character vector: custom factor levels for metrics; this controls the
  order that metric statistics are plotted

- return.tbl:

  boolean: if TRUE, return the data frames of results used to make the
  ggplot instead of the plot itself

## Value

A ggplot of evaluation statistics.

## Details

In this plot, the x-axis represents a tuning parameter range, while the
y-axis represents the average of a statistic over all partitions.
Different colors represent the categories or values of another tuning
parameter. Error bars represent the standard deviation of a statistic
around the mean. Currently, this function can only plot two tuning
parameters at a time.

## Examples

``` r
if (FALSE) { # \dontrun{
# first, let's tune some models
occs <- read.csv(file.path(system.file(package="predicts"), 
"/ex/bradypus.csv"))[,2:3]
envs <- rast(list.files(path=paste(system.file(package="predicts"), 
"/ex", sep=""), pattern="tif$", full.names=TRUE))
bg <- as.data.frame(predicts::backgroundSample(envs, n = 10000))
names(bg) <- names(occs)
 
ps <- list(orientation = "lat_lat")
e <- ENMevaluate(occs, envs, bg, 
               tune.args = list(fc = c("L","LQ","LQH"), rm = 1:5), 
               partitions = "block", partition.settings = ps, 
               algorithm = "maxnet", categoricals = "biome", 
               parallel = TRUE)
               
evalplot.stats(e, c("cbi.val", "or.mtp"), x.var = "rm", color.var = "fc", 
             error.bars = FALSE)
} # }
```
