% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enm.maxnet.R
\docType{data}
\name{enm.maxnet}
\alias{enm.maxnet}
\title{ENMdetails maxnet}
\format{
An object of class \code{ENMdetails} of length 1.
}
\usage{
enm.maxnet
}
\description{
This is the ENMdetails implementation for maxnet, the R version of
the Maxent algorithm. The configuration for running the model now includes addsamplestobackground = TRUE,
which explicitly adds presences to the background for model training, though as the current 
version of maxnet has this set to TRUE as default, behavior between ENMeval versions should not differ.
}
\keyword{datasets}
