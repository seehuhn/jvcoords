\name{standardize}
\alias{standardize}
\title{function to standardize data}
\description{Standardizes each column of the data and returns the
  results as an object of class \code{coords}.}
\usage{
  standardize(x)
}
\arguments{
  \item{x}{A data matrix, rows are observations, columns are variables.}
}
\details{
  This function standardizes the columns of \code{x} by subtracting
  the mean of each column and then dividing by the standard devition.
  The transformed data is stored in the \code{$y} field of the
  returned \code{coords} object.
}
\value{
  An object of class \code{coords}.
}
\author{
  Jochen Voss <voss@seehuhn.de>
}
\seealso{
  \code{\link{coords}};
  alternative implementation \code{\link{scale}}
}
\examples{
  w <- standardize(iris[,1:4])
  colMeans(w$y)
  apply(w$y, 2, sd)
}
