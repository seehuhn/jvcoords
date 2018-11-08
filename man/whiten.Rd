\name{whiten}
\alias{whiten}
\title{Whiten data}
\description{Performs whitening of data and returns the results as an object of
  class \code{coords}.}
\usage{
  whiten(x, compute.scores = TRUE)
}
\arguments{
  \item{x}{A data matrix, rows are observations, columns are variables.}
  \item{compute.scores}{Whether to compute the scores (i.e. \code{x} in the new
    basis).}
}
\details{
  This function whitens the data by finding an affine transformation such
  that the transformed data has mean 0 and identity covariance matrix.

  If \code{compute.scores} is set to \code{FALSE}, only the
  information required for the \code{toCoords()} and
  \code{fromCoords()} to work is stored in the returned \code{coords}
  object; otherwise the scores (transformed data) will be stored in
  the \code{$y} field of the \code{coords} object.
}
\value{
  An object of class \code{coords}, with the following
  additional components added:
  \item{loadings}{the loadings, each column is one of the new basis vectors}
  \item{y}{if \code{compute.scores==TRUE}, this is \code{x} expressed in the new basis}
}
\author{
  Jochen Voss <voss@seehuhn.de>
}
\seealso{
  \code{\link{coords}}
}
\examples{
  w <- whiten(iris[,1:4])
  colMeans(w$y)
  round(cov(w$y), 3)
}
