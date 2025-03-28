% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{evalplot.envSim.map}
\alias{evalplot.envSim.map}
\title{Similarity maps for partition groups}
\usage{
evalplot.envSim.map(
  e = NULL,
  envs,
  occs.z = NULL,
  bg.z = NULL,
  occs.grp = NULL,
  bg.grp = NULL,
  ref.data = "occs",
  envs.vars = NULL,
  bb.buf = 0,
  occs.testing.z = NULL,
  plot.bg.pts = FALSE,
  sim.palette = NULL,
  pts.size = 1.5,
  gradient.colors = c("red", "white", "blue"),
  na.color = "gray",
  return.tbl = FALSE,
  return.ras = FALSE,
  quiet = FALSE
)
}
\arguments{
\item{e}{ENMevaluation object (optional)}

\item{envs}{SpatRaster: environmental predictor variables used to build the models in "e"; categorical variables will be 
removed before internally as they cannot be used to calculate MESS}

\item{occs.z}{data frame: longitude, latitude, and environmental predictor variable values for occurrence records, in that order (optional);
the first two columns must be named "longitude" and "latitude"}

\item{bg.z}{data frame: longitude, latitude, and environmental predictor variable values for background records, in that order (optional);
the first two columns must be named "longitude" and "latitude"}

\item{occs.grp}{numeric vector: partition groups for occurrence records (optional)}

\item{bg.grp}{numeric vector: partition groups for background records (optional)}

\item{ref.data}{character: the reference to calculate MESS based on occurrences ("occs") or background ("bg"), with default "occs"}

\item{envs.vars}{character vector: names of a predictor variable to plot similarities for; if left NULL, calculations are done
with respect to all variables (optional)}

\item{bb.buf}{numeric: distance used to buffer (extend) the mapping extent in map units; for latitude/longitude, this is in degrees (optional)}

\item{occs.testing.z}{data frame: longitude, latitude, and environmental predictor variable values for fully withheld testing records, 
in that order; this is for use only with the "testing" partition option when an ENMevaluation object is not input (optional)}

\item{plot.bg.pts}{boolean: if TRUE, plot background points when using ref.data = "bg"}

\item{sim.palette}{character: RColorBrewer palette name to use for plotting discrete variables; if NULL, default is "Set1"}

\item{pts.size}{numeric: custom point size for ggplot}

\item{gradient.colors}{character vector: colors used for ggplot2::scale_fill_gradient2}

\item{na.color}{character: color used for NA values}

\item{return.tbl}{boolean: if TRUE, return the data frames of similarity values used to make the ggplot instead of the plot itself}

\item{return.ras}{boolean: if TRUE, return the SpatRaster of similarity values used to make the ggplot instead of the plot itself}

\item{quiet}{boolean: if TRUE, silence all function messages (but not errors)}
}
\value{
A ggplot of environmental similarities between the occurrence or background data 
for each partition and all predictor variable values in the extent.
}
\description{
Maps environmental similarity of reference partitions (occurrences or 
background) to all cells with values in the predictor variable rasters. This function uses 
raster data, and thus cannot map similarity values using only tables of environmental values f
or occurrences or background. Further, this function does not calculate similarity for 
categorical variables.
}
\details{
When fully withheld testing groups are used, make sure to input either an ENMevaluation 
object or the argument occs.testing.z. In the resulting plot, partition 1 refers to the training data,
while partition 2 refers to the fully withheld testing group.

Rasters are plotted showing the environmental similarity estimates for each 
partition group in relation to the others. The similarity between environmental values associated with the 
validation occurrence or background records per partition group and those associated with 
the entire study extent (specified by the extent of the input SpatRaster "envs") are 
calculated with the MESS algorithm, and the minimum similarity per grid cell is returned. Higher 
negative values indicate greater environmental difference between the validation occurrences 
and the study extent, and higher positive values indicate greater similarity. This function 
uses the `mess()` function from the package `predicts` to calculate the similarities. Please see the below 
reference for details on MESS.
}
\examples{
\dontrun{
library(terra)
library(ENMeval)

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

# now, plot the environmental similarity of each partition to the others               
evalplot.envSim.map(e, envs)
}

}
\references{
Baumgartner J, Wilson P (2021). _rmaxent: Tools for working with Maxent in R_. R package version 0.8.5.9000, <URL: https://github.com/johnbaums/rmaxent>.
Elith, J., Kearney, M., and Phillips, S. (2010) The art of modelling range-shifting species. \emph{Methods in Ecology and Evolution}, \bold{1}: 330-342. \doi{doi:10.1111/j.2041-210X.2010.00036.x}
}
