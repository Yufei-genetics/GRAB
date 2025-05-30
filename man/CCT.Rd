% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CCT.R
\name{CCT}
\alias{CCT}
\title{An analytical p-value combination method using the Cauchy distribution}
\usage{
CCT(pvals, weights = NULL)
}
\arguments{
\item{pvals}{a numeric vector of p-values, where each of the element is
between 0 to 1, to be combined. If a p-value equals to 1, we set it as 0.999. If a p-value equals to 0, an error action is executed.}

\item{weights}{a numeric vector of non-negative weights. If \code{NULL}, the
equal weights are assumed.}
}
\value{
the aggregated p-value combining p-values from the vector \code{pvals}.
}
\description{
The \code{CCT} function takes in a numeric vector of p-values, a numeric
vector of non-negative weights, and return the aggregated p-value using Cauchy method.
}
\examples{
pvalues <- c(2e-02,4e-04,0.2,0.1,0.8)
CCT(pvals=pvalues)
}
\references{
Liu, Y., & Xie, J. (2020). Cauchy combination test: a powerful test
with analytic p-value calculation under arbitrary dependency structures.
\emph{Journal of the American Statistical Association 115}(529), 393-402.
(\href{https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1554485}{pub})
}
