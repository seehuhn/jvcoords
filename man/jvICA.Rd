\name{jvICA}
\alias{jvICA}
\title{Perform Independent Component Analysis (ICA)}
\description{Performs an independent components analysis on the given data
  matrix and returns the results as an object of class \code{coords}.}
\usage{
  jvICA(x, n.comp, maxit = 5000)
}
\arguments{
  \item{x}{A data matrix, rows are samples, columns are variables.}
  \item{n.comp}{How many independent components to compute.}
  \item{maxit}{Number of iterations to use for minimizing the entropy.}
}
\details{
  This function performs Independent Component Analysis (PCA) on the
  data.
}
\value{
  An object of class \code{coords}.
}
\author{
  Jochen Voss <voss@seehuhn.de>
}
\seealso{
  \code{\link{coords}}
}
