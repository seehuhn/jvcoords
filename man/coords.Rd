\name{coords}
\alias{coords}
\alias{fromCoords}
\alias{toCoords}
\title{Affine Coordinate Transforms}
\description{Perform affine coordinate transformations.}
\usage{
  coords(loadings=NULL, pre.sub=0, pre.div=1, post.mul=1)
  toCoords(trans, x)
  fromCoords(trans, y)
}
\arguments{
  \item{loadings}{The loadings, as a matrix with orthogonal columns.}
  \item{pre.sub}{pre-transformation centring}
  \item{pre.div}{pre-transformation scaling}
  \item{post.mul}{post-transformation scaling}
  \item{trans}{an object of class \code{coords}}
  \item{x}{Data matrix, rows are observations, columns are variables.}
  \item{y}{Transformed data matrix, rows are observations, columns are variables.}
}
\details{
  The function \code{coords()} creates an object representing a coordinate
  transformation.
}
\value{
  An object of class \code{coords}, as a list with the following
  components:
  \item{p}{the number of variables in the original data set}
  \item{q}{the number of variables in the transformed data set}
  \item{loadings}{the loadings: each column is one of the new basis vectors}
  \item{pre.sub}{pre-transformation centring}
  \item{pre.div}{pre-transformation scaling}
  \item{post.mul}{post-transformation scaling}
}
\author{
  Jochen Voss <voss@seehuhn.de>
}
\seealso{
  \code{\link{standardize}},
  \code{\link{whiten}},
  \code{\link{jvPCA}}
}
\examples{
  pc <- jvPCA(iris[,1:4], n.comp=3)
  toCoords(pc, c(5, 3, 4, 1))
  fromCoords(pc, c(1, 1, 0))
}
