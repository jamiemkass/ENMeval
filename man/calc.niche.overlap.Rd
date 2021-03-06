% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{calc.niche.overlap}
\alias{calc.niche.overlap}
\title{Calculate Similarity of ENMs in Geographic Space}
\usage{
calc.niche.overlap(predictors, overlapStat, quiet = FALSE)
}
\arguments{
\item{predictors}{RasterStack: at least 2 Maxent raster predictions}

\item{overlapStat}{character: either "D" or "I", the statistic calculated by the \code{nicheOverlap} function of the \pkg{dismo} package (default: "D")}

\item{quiet}{boolean: if TRUE, silence all function messages (but not errors)}
}
\value{
A matrix with the lower triangle giving values of pairwise "niche overlap" in geographic space.  Row and column names correspond to the results table output by \code{\link{ENMevaluate}()}.
}
\description{
Compute pairwise "niche overlap" in geographic space for Maxent predictions. The value ranges from 0 (no overlap) to 1 (identical predictions).  The function uses the \code{nicheOverlap} function of the \pkg{dismo} package (Hijmans \emph{et al.} 2011).
}
\details{
"D" refers to Schoeners \emph{D} (Schoener 1968), while "I" refers to the \emph{I} similarity statistic from Warren \emph{et al.} (2008).
}
\references{
Hijmans, R. J., Phillips, S., Leathwick, J. & Elith, J. (2011) dismo package for R. Available online at: \url{https://cran.r-project.org/package=dismo}.

Schoener, T. W. (1968) The \emph{Anolis} lizards of Bimini: resource partitioning in a complex fauna. \emph{Ecology}, \bold{49}: 704-726. \doi{10.2307/1935534}

Warren, D. L., Glor, R. E., Turelli, M. & Funk, D. (2008) Environmental niche equivalency versus conservatism: quantitative approaches to niche evolution. \emph{Evolution}, \bold{62}: 2868-2883. \doi{10.1111/j.1558-5646.2008.00482.x}
}
\seealso{
`nicheOverlap` in the \pkg{dismo} package
}
\author{
Based on \pkg{dismo}::\code{nicheOverlap}, which is based on \pkg{SDMTools}::\code{Istat}
Robert Muscarella <bob.muscarella@gmail.com>
}
