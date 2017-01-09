\name{eval2}
\alias{eval2}
\docType{data}
\title{ An object of class "ENMevaluation" }
\description{ An example results file based on a call of ENMevaluate for use in the ENMeval vignette. }
\usage{ data(eval2) }
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
The dataset is used for the ENMeval vignette.
}
