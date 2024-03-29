% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enm.maxent.jar.R
\docType{data}
\name{enm.maxent.jar}
\alias{enm.maxent.jar}
\title{ENMdetails maxent.jar}
\format{
An object of class \code{ENMdetails} of length 1.
}
\usage{
enm.maxent.jar
}
\description{
This is the ENMdetails implementation for maxent.jar, the Java version of
the Maxent algorithm. The configuration for running the model differs slightly from that
in previous versions of ENMeval (0.3.0 and before) in that this version (>=2.0.0) uses the
default of adding presences to the background for model training, while previous versions
had turned this off. Specifically, previous versions ran maxent() with "noaddsamplestobackground"
in the "args" vector argument, while this version does not.
}
\keyword{datasets}
