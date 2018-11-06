\name{coords}
\alias{coords}
\alias{AppendTrfm}
\alias{FromCoords}
\alias{ToCoords}
\title{Affine Coordinate Transforms}
\description{Perform affine coordinate transformations.}
\usage{
  coords(p, name = NULL)
  AppendTrfm(trfm, op = c("add", "mult", "orth"), val)
  ToCoords(trfm, x)
  FromCoords(trfm, y)
}
\arguments{
  \item{p}{The dimension of the original data}
  \item{name}{A short name for the coordinate transformation (optional)}
  \item{trfm}{an object of class \code{coords}}
  \item{x}{Data matrix, rows are observations, columns are variables.}
  \item{y}{Transformed data matrix, rows are observations, columns are variables.}
  \item{loadings}{The loadings, as a matrix with orthogonal columns.}
}
\details{
  The function \code{coords()} creates an object representing an affine
  coordinate transformation.
}
\value{
  An object of class \code{coords}, as a list with the following
  components:
  \item{p}{the number of variables in the original data set}
  \item{q}{the number of variables in the transformed data set}
  \item{name}{the name of the transformation}
  \item{loadings}{the loadings: each column is one of the new basis vectors}
}
\author{
  Jochen Voss <voss@seehuhn.de>
}
\seealso{
  \code{\link{standardize}},
  \code{\link{whiten}},
  \code{\link{PCA}}
}
\examples{
  pc <- PCA(iris[,1:4], n.comp=3)
  ToCoords(pc, c(5, 3, 4, 1))
  FromCoords(pc, c(1, 1, 0))
}
