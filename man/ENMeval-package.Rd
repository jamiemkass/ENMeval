\name{ENMeval-package}
\alias{ENMeval-package}
\alias{ENMeval}
\docType{package}
\title{ Automated runs and evaluations of ecological niche models }
\description{Automatically partitions data into bins for model training and testing, executes ecological niche models (ENMs) across a range of user-defined settings, and calculates evaluation metrics to help achieve a balance between goodness-of-fit and model complexity.}
\details{
\tabular{ll}{
Package: \tab ENMeval\cr
Type: \tab Package\cr
Version: \tab 0.1.2\cr
Date: \tab 2014-08-25\cr
License: \tab GNU 3.0\cr
}

The \pkg{ENMeval} package (Muscarella \emph{et al.} 2014) (1) automatically partitions data into training and testing bins using one of six methods (including several options for spatially independent partitions as well as user-defined bins), (2) executes a series of ENMs using Maxent (Phillips \emph{et al.} 2006) with a variety of user-defined settings (i.e., feature classes and regularization multipliers), conducting \emph{k}-fold cross validation, and (3) calculates multiple evaluation metrics to aid in selecting model settings that balance model goodness-of-fit and complexity (i.e., "model tuning" or "smoothing").

\code{\link{ENMevaluate}} is the primary function of the \pkg{ENMeval} package, and multiple other functions highlighted below are called when it is run.  The six options for partitioning occurrence data into training and testing (i.e., calibration and evaluation) bins are: \emph{n}-1 jackknife, random \emph{k}-fold, user-specified bins, and three explicit methods of masked geographically structured \emph{k}-fold partitioning (see: \code{\link{get.evaluation.bins}}).  After model training, these bins are used to calculate five metrics of model performance for each combination of settings: model discrimination (AUC of test localities), the difference between training and testing AUC, two different threshold-based omission rates, and the small sample-size corrected version of the Akaike information criterion (AICc), the latter using the unpartitioned dataset.  A model prediction (as a raster layer) using the full (unpartitioned) dataset is generated for each combination of feature class and regularization multiplier settings.  Similarity of these models in geographic space (i.e., "niche overlap") can be calculated to better understand how model settings change predictions (see \code{\link{calc.niche.overlap}}).  The results of \code{ENMevaluate} are returned as an object of class \code{\link{ENMevaluation-class}}.  A basic plotting function (\code{\link{eval.plot}}) can be used to visualize how evaluation metrics depend on model settings.

As of version 0.1.2, \code{\link{ENMevaluate}} includes an option for parallel computing.  Setting \code{parallel = TRUE} can significantly speed up processing time, particularly for large analyses.  For very small analyses, it may actually take longer than running with \code{parallel = FALSE}.

}

\note{
Currently, \pkg{ENMeval} only implements the Maxent algorithm, but we eventually plan to expand it to work with other algorithms.  All calculations are based on the raw Maxent output (i.e., \emph{not} logistic or cumulative transformations) and users can choose whether to use 'clamping' (see Maxent documentation for details on this option).  Additionally, Maxent models are run with the arguments: \code{noaddsamplestobackground} and \code{noremoveDuplicates}.  Users should consult Maxent documentation (Phillips \emph{et al.} 2006) and other references (e.g., Phillips and Dudik 2008) for more information on these options.  We note that interested users can edit the source code of \code{ENMeval} (in particular, the \code{\link{make.args}} and \code{\link{tuning}} functions) if they desire to change these or other options.

\code{ENMevaluate} directly uses several functions from the \pkg{dismo} package (Hijmans \emph{et al.} 2011), the most important of which is the \code{maxent} function that runs the Maxent algorithm (Phillips \emph{et al.} 2006) in Java.  Before running this command, the user must first download Maxent from \href{http://www.cs.princeton.edu/~schapire/maxent/}{this website}.  Then, place the file 'maxent.jar' in the 'java' folder of the \pkg{dismo} package. The user can locate that folder by typing: \code{system.file("java", package="dismo")}.  For additional details, users should consult the documentation of the \pkg{dismo} package.

}

\author{
Robert Muscarella, Peter J. Galante, Mariano Soley-Guardia, Robert A. Boria, Jamie M. Kass, Maria Uriarte and Robert P. Anderson

Maintainer: Robert Muscarella <bob.muscarella@gmail.com>
}
\references{ 
Hijmans, R. J., Phillips, S., Leathwick, J. and Elith, J. 2011. dismo package for R. Available online at: \url{http://cran.r-project.org/web/packages/dismo/index.html}.

Muscarella, R., Galante, P. J., Soley-Guardia, M., Boria, R. A., Kass, J. M., Uriarte, M., and Anderson, R. P. 2014. ENMeval: An R package for conducting spatially independent evaluations and estimating optimal model complexity for Maxent ecological niche models. \emph{Methods in Ecology and Evolution}, \bold{5}: 1198-1205.

Phillips, S. J., Anderson, R. P., and Schapire, R. E. 2006. Maximum entropy modeling of species geographic distributions. \emph{Ecological Modelling}, \bold{190}: 231-259.

Phillips, S. J., and Dudik, M. 2008. Modeling of species distributions with Maxent: new extensions and a comprehensive evaluation. \emph{Ecography}, \bold{31}: 161-175.

}
\keyword{ niche }
\keyword{ ENM }
\keyword{ SDM }
\seealso{
\code{maxent} in the \pkg{dismo} package
}
