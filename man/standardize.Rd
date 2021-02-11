\name{standardize}
\alias{standardize}
\title{Standardize data}
\description{Standardize each column of a data matrix and return the
  results as an object of class \code{coords}.}
\usage{
  standardize(x, compute.scores = TRUE)
}
\arguments{
  \item{x}{A data matrix, rows are observations, columns are variables.}
  \item{compute.scores}{Whether to compute the scores (i.e. \code{x}
    in the new basis).}
}
\details{
  This function standardizes the columns of \code{x} by subtracting the
  mean of each column and then dividing by the standard deviation.  The
  transformed data is stored in the \code{$y} field of the returned
  \code{coords} object.

  If \code{compute.scores} is set to \code{FALSE}, only the information
  required for the \code{toCoords()} and \code{fromCoords()} to work is
  stored in the returned \code{coords} object; otherwise the scores
  (transformed data) will be stored in the \code{$y} field of the
  \code{coords} object.
}
\value{
  An object of class \code{coords}, with the following additional
  components added:
  \item{y}{if \code{compute.scores==TRUE}, this is \code{x} expressed in
    the new basis}
}
\author{
  Jochen Voss <voss@seehuhn.de>
}
\seealso{
  \code{\link{coords}};
  alternative implementation \code{\link{scale}}
}
\examples{
  w <- standardize(iris[, 1:4])
  colMeans(w$y)
  apply(w$y, 2, sd)
}
