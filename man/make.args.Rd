\name{make.args }
\alias{make.args}
\title{ Generate arguments for Maxent }
\description{ This function generates a list of arguments to pass to \code{Maxent} or to use as convenient labels for plotting. }
\usage{
make.args(RMvalues = seq(0.5, 4, 0.5), 
	fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), 
	labels = FALSE)
}
\arguments{

  \item{RMvalues}{ Vector of (non-negative) values to use for the regularization multiplier. }

  \item{fc}{ Character vector of feature class combinations to be included in analysis. }

  \item{labels}{ logical; If \code{FALSE} (default), provides arguments to pass directly to Maxent; if  \code{TRUE}, provides more intuitive labels to use, for example, in plotting.}
}

\details{
When \code{labels = FALSE}, the following additional arguments are added: 

\code{noaddsamplestobackground, noremoveDuplicates, noautofeature}.

For details on these arguments, see Phillips \emph{et al.} (2006) and the help documentation and tutorial of the Maxent software and the tutorial that can be downloaded from \href{https://www.cs.princeton.edu/~schapire/maxent/}{this website}.
}

\value{
If \code{labels = FALSE}, a \code{list} the length of the total number of unique combinations of feature class(es) and regularization multiplier(s).

If \code{labels = TRUE}, a \code{list} of two items:
\item{$}{character vector of feature class combinations in the same order they were provided.}
\item{$}{numeric vector of regularization multiplier values in the same order they were provided.}

}

\author{ Robert Muscarella <bob.muscarella@gmail.com> and Jamie M. Kass <jkass@gc.cuny.edu> }

\references{ 

Phillips, S. J., Anderson, R. P. and Schapire, R. E. 2006. Maximum entropy modeling of species geographic distributions. \emph{Ecological Modelling}, \bold{190}: 231-259.

}

\seealso{
\code{maxent} in the \pkg{dismo} package.
}

\examples{
make.args(RMvalues=c(1:3), fc=c("L","LQ"))

make.args(RMvalues=c(1:3), fc=c("L","LQ"), labels=TRUE)
}
