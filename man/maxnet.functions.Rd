\name{ maxnet.functions  }

\alias{maxentJARversion.R}
\alias{modelTune.maxentJar.R}
\alias{modelTune.maxnet.R}
\alias{predict.maxnetRas.R}

\title{ Functions for compatability with maxnet package }

\description{ These functions are used internally for compatability with the \pkg{maxnet} package (Phillips \emph{et al.} 2017, Phillips 2017). }

\usage{
  maxentJARversion.R ()
  modelTune.maxentJar.R (pres, bg, env, nk, group.data, args.i, userArgs, rasterPreds, clamp)
  modelTune.maxnet.R (pres, bg, env, nk, group.data, args.i, rasterPreds, clamp)
  predict.maxnetRas.R (mod, env, type, clamp)
}

\arguments{
  \item{pres}{ Occurrence points. }
  
  \item{bg}{ Background points. }
  
  \item{env}{ Environmental predictor variables. }
  
  \item{nk}{ Number of k-fold partitions. }
  
  \item{group.data}{ Input data grouped for k-fold evaluations (output of the \code{get.evaluation.bins} functions). }
  
  \item{args.i}{ Internal arguments. }
  
  \item{userArgs}{ User arguments.}
  
  \item{rasterPreds}{ Raster(s) of model predictions. }
  
  \item{clamp}{ Logical. }
  
  \item{mod}{ maxnet model object. }
  
  \item{type}{ see \pkg{maxnet} package documentation. }
}

\details{
  These functions are used internally for compatability with the \pkg{maxnet} package.}

\value{
  Depends on which function is used.
}


\references{ 
  Phillips, S. J. 2017. maxnet package for R. Available online at: \url{https://CRAN.R-project.org/package=maxnet}.
  
  Phillips, S. J., Anderson, R. P., Dudík, M., Schapire, R. E. and Blair, M. E. 2017. Opening the black box: an open-source release of Maxent. \emph{Ecography}, \bold{40}: 887–893.
  }

\author{ Robert Muscarella <bob.muscarella@gmail.com> and Jamie M. Kass <jkass@gc.cuny.edu> }
