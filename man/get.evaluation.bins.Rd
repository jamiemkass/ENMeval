\name{get.evaluation.bins }

\alias{get.evaluation.bins}
\alias{get.block}
\alias{get.checkerboard1}
\alias{get.checkerboard2}
\alias{get.randomkfold}
\alias{get.user}
\alias{get.jackknife}

\title{ Methods to partition data for evaluation }

\description{ \code{ENMeval} provides six methods to partition occurrence and background localities into bins for training and testing (or, evaluation and calibration).  Users should carefully consider the objectives of their study and the influence of spatial bias when deciding on a method of data partitioning. }

\usage{
get.block (occ, bg.coords)
get.checkerboard1(occ, env, bg.coords, aggregation.factor)
get.checkerboard2(occ, env, bg.coords, aggregation.factor)
get.jackknife(occ, bg.coords)
get.randomkfold(occ, bg.coords, kfolds)
get.user(occ.grp, bg.grp)
}

\arguments{
\item{occ}{ Two-column matrix or data.frame of longitude and latitude (in that order) of occurrence localities. }

\item{bg.coords}{ Two-column matrix or data.frame of longitude and latitude (in that order) of background localities. }

\item{env}{ RasterStack of environmental predictor variables. }

\item{aggregation.factor}{ A vector or list of 1 or 2 numbers giving the scale for aggregation used for the \code{get.checkerboard1} and \code{get.checkerboard2} methods.  If a single number is given and \code{get.checkerboard2} partitioning method is used, the single value is used for both scales of aggregation.}

\item{kfolds}{ Number of random \emph{k}-folds for \code{get.randomkfold} method. }

\item{occ.grp}{ Vector of user-defined bins for occurrence localities for \code{get.user} method. }

\item{bg.grp}{ Vector of user-defined bins for background localities for \code{get.user} method. }
}

\details{
These functions are used internally to partition data during a call of \code{\link{ENMevaluate}}.

The \code{get.block} method partitions occurrence localities by finding the latitude and longitude that divide the occurrence localities into four groups of (insofar as possible) equal numbers.  Background localities are assigned to each of the four groups based on their position with respect to these lines.  While the \code{get.block} method results in (approximately) equal division of occurrence localities among four groups, the number of background localities (and, consequently, environmental and geographic space) in each group depends on the distribution of occurrence localities across the study area.

The \code{get.checkerboard1} and \code{get.checkerboard2} methods are variants of a checkerboard approach to partition occurrence localities.  These methods use the \code{gridSample} function of the \pkg{dismo} package (Hijmans \emph{et al.} 2011) to partition records according to checkerboard grids across the study extent.  The spatial grain of these grids is determined by resampling (or aggregating) the original environmental input grids based on the user-defined \code{aggregation factor} (e.g., an aggregation factor of 2 results in a checkerboard with grid cells four times as large in area as the original input grids).  The \code{get.checkerboard1} method partitions data into two groups according to a single checkerboard pattern, and the \code{get.checkerboard2} method partitions data into four groups according to two nested checkerboard grids.  In contrast to the \code{get.block} method, both the \code{get.checkerboard1} and \code{get.checkerboard2} methods subdivide geographic space equally but do not ensure a balanced number of occurrence localities in each group.  The two \code{get.checkerboard} methods give warnings (and potentially errors) if zero points (occurrence or background) fall in any of the expected bins.

The \code{get.jackknife} method is a special case of \emph{k}-fold cross validation where the number of bins (\emph{k}) is equal to the number of occurrence localities (\emph{n}) in the dataset.  It is suggested for datasets of relatively small sample size (generally < 25 localities) (Pearson \emph{et al.} 2007; Shcheglovitova and Anderson 2013).

The \code{get.randomkfold} method partitions occurrence localities randomly into a user-specified number of (\emph{k}) bins.   This is equivalent to the method of \emph{k}-fold cross valiation currently provided by Maxent. 

The \code{get.user} method is flexible and enables users to define evaluation bins \emph{a priori}.  With this method, occurrence and background localities, as well as evaluation bin designation for each locality, are supplied by the user.

}

\value{
A named list of two items:
\item{$occ.grp}{ A vector of bin designation for occurrence localities in the same order they were provided.}
\item{$bg.grp}{ A vector of bin designation for background localities in the same order they were provided.}
}

\note{
The \code{checkerboard1} and \code{checkerboard2} methods are designed to partition occurrence localities into two and four evaluation bins, respectively.  They may give fewer bins, however, depending on where the occurrence localities fall with respect to the grid cells (e.g., all records happen to fall in the "black" squares).  A warning is given if the number of bins is < 4 for the \code{checkerboard2} method, and an error is given if all localities fall into a single evaluation bin.
}

\references{ 
Hijmans, R. J., Phillips, S., Leathwick, J. and Elith, J. 2011. dismo package for R. Available online at: \url{http://cran.r-project.org/web/packages/dismo/index.html}.

Pearson, R. G., Raxworthy, C. J., Nakamura, M. and Peterson, A. T. 2007. Predicting species distributions from small numbers of occurrence records: a test case using cryptic geckos in Madagascar. \emph{Journal of Biogeography}, \bold{34}: 102-117.

Shcheglovitova, M. and Anderson, R. P. (2013) Estimating optimal complexity for ecological niche models: a jackknife approach for species with small sample sizes. \emph{Ecological Modelling}, \bold{269}: 9-17.
}

\author{ Robert Muscarella <bob.muscarella@gmail.com> and Jamie M. Kass <jkass@gc.cuny.edu> }

\examples{

set.seed(1)

### Create environmental extent (raster) 
env <- raster(matrix(nrow=25, ncol=25))

### Create presence localities
set.seed(1)
nocc <- 25
xocc <- rnorm(nocc, sd=0.25) + 0.5
yocc <- runif(nocc, 0, 1)
occ.pts <- as.data.frame(cbind(xocc, yocc))

### Create background points
nbg <- 500
xbg <- runif(nbg, 0, 1)
ybg <- runif(nbg, 0, 1)
bg.pts <- as.data.frame(cbind(xbg, ybg))

### Show points
plot(env)
points(bg.pts)
points(occ.pts, pch=21, bg=2)

### Block partitioning method
blk.pts <- get.block(occ.pts, bg.pts)
plot(env)
points(occ.pts, pch=23, bg=blk.pts$occ.grp)
plot(env)
points(bg.pts, pch=21, bg=blk.pts$bg.grp)

### Checkerboard1 partitioning method
chk1.pts <- get.checkerboard1(occ.pts, env, bg.pts, 4)
plot(env)
points(occ.pts, pch=23, bg=chk1.pts$occ.grp)
plot(env)
points(bg.pts, pch=21, bg=chk1.pts$bg.grp)

### Checkerboard2 partitioning method
chk2.pts <- get.checkerboard2(occ.pts, env, bg.pts, c(2,2))
plot(env)
points(occ.pts, pch=23, bg=chk2.pts$occ.grp)
plot(env)
points(bg.pts, pch=21, bg=chk2.pts$bg.grp)

### Random k-fold partitions
# Note that k random does not partition the background
krandom.pts <- get.randomkfold(occ.pts, bg.pts, 4)
plot(env)
points(occ.pts, pch=23, bg=krandom.pts$occ.grp)
plot(env)
points(bg.pts, pch=21, bg=krandom.pts$bg.grp)

### k-1 jackknife partitions
# Note background is not partitioned
jack.pts <- get.jackknife(occ.pts, bg.pts)
plot(env)
points(occ.pts, pch=23, bg=rainbow(length(jack.pts$occ.grp)))
plot(env)
points(bg.pts, pch=21, bg=jack.pts$bg.grp)

### User-defined partitions
# Note background is not partitioned
occ.grp <- c(rep(1, 10), rep(2, 5), rep(3, 10))
bg.grp <- c(rep(1, 200), rep(2, 100), rep(3, 200))
user.pts <- get.user(occ.grp, bg.grp)
plot(env)
points(occ.pts, pch=23, bg=user.pts$occ.grp)
plot(env)
points(bg.pts, pch=21, bg=user.pts$bg.grp)

}
