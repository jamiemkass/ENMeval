\name{corrected.var }
\alias{corrected.var}
\title{ Calculate variance corrected for non-independence of \emph{k}-fold iterations }
\description{This function calculates variance corrected for non-independence of \emph{k}-fold iterations.  See Appendix of Shcheglovitova & Anderson (2013) and other references (Miller 1974; Parr 1985; Shao and Wu 1989) for additional details. }
\usage{ corrected.var(x, nk) }
\arguments{

  \item{x}{ A numeric vector. }

  \item{nk}{ Number of \emph{k}-fold iterations. }

}
\value{ A numeric value of the corrected variance. }
\details{ This function calculates variance that is corrected for the non-independence of \emph{k} cross-validation iterations.  Following Shao and Wu (1989): 

	\deqn{Sum Of Squares * ((n-1)/n)} 

where \emph{n} = the number of \emph{k}-fold iterations.
}
\references{
Miller, R. G. (1974) The jackknife - a review. \emph{Biometrika}, \bold{61}: 1-15.

Parr, W. C. (1985) Jackknifing differentiable statistical functionals. \emph{Journal of the Royal Statistics Society, Series B}, \bold{47}: 56-66.

Shao J. and Wu, C. F. J. (1989) A general theory for jackknife variance estimation. \emph{Annals of Statistics}, \bold{17}: 1176-1197.

Shcheglovitova, M. and Anderson, R. P. (2013) Estimating optimal complexity for ecological niche models: a jackknife approach for species with small sample sizes. \emph{Ecological Modelling}, \bold{269}: 9-17.

}
\author{ Robert Muscarella <bob.muscarella@gmail.com> }
