\name{enmeval_results }
\alias{enmeval_results}
\docType{data}
\title{ An object of class "ENMevaluation" }
\description{ An example results file based on a call of ENMevaluate (see example). }
\usage{ data(enmeval_results) }
\format{
An object of class 'ENMevaluation' with nine slots:

\code{@ results}         : data.frame of evaluation metrics

\code{@ predictions}     : RasterStack of model predictions

\code{@ models}: list of \code{MaxEnt} model objects (see \code{MaxEnt} documentation for details)

\code{@ partition.method}: character giving method of data partitioning

\code{@ occ.pts}         : data.frame of latitude and longitude of occurrence localities

\code{@ occ.grp}         : data.frame of bins for occurrence localities

\code{@ bg.pts}          : data.frame of latitude and longitude of background localities

\code{@ bg.grp}          : data.frame of bins for background localities

\code{@ overlap}         : matrix of pairwise niche overlap

}

\details{
The dataset is based on the simulated dataset and call of \code{\link{ENMevaluate}} shown in the example section below.
}

\examples{
require(raster)

### Simulated data environmental covariates
set.seed(1)
r1 <- raster(matrix(nrow=50, ncol=50, data=runif(10000, 0, 25)))
r2 <- raster(matrix(nrow=50, ncol=50, data=rep(1:100, each=100), byrow=TRUE))
r3 <- raster(matrix(nrow=50, ncol=50, data=rep(1:100, each=100)))
r4 <- raster(matrix(nrow=50, ncol=50, data=c(rep(1,1000),rep(2,500)),byrow=TRUE))
values(r4) <- as.factor(values(r4))
env <- stack(r1,r2,r3,r4)

### Simulate occurrence localities
nocc <- 50
x <- (rpois(nocc, 2) + abs(rnorm(nocc)))/11
y <- runif(nocc, 0, .99)
occ <- cbind(x,y)

\dontrun{
### This gives the results that are loaded below:
enmeval_results <- ENMevaluate(occ, env, method="block", n.bg=500, 
  categoricals=4, algorithm='maxent.jar')
}

data(enmeval_results)
enmeval_results

### See table of evaluation metrics
enmeval_results@results

### Plot prediction with lowest AICc
plot(enmeval_results@predictions[[which (enmeval_results@results$delta.AICc == 0) ]])
points(enmeval_results@occ.pts, pch=21, bg= enmeval_results@occ.grp)

### Niche overlap statistics between model predictions
enmeval_results@overlap
}
