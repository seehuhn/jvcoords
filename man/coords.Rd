\name{coords}
\alias{coords}
\alias{fromCoords}
\alias{toCoords}
\title{Affine Coordinate Transforms}
\description{Perform affine coordinate transformations.}
\usage{
  coords(loadings, pre.sub=0, pre.div=1, post.mul=1)
  coords.print(trans, ...)
  toCoords(trans, x)
  fromCoords(trans, y)
}
\arguments{
  \item{loadings}{The loadings, as a matrix with orthogonal columns.}
  \item{pre.sub}{pre-transformation centring}
  \item{pre.div}{pre-transformation scaling}
  \item{post.mul}{post-transformation scaling}
  \item{trans}{an object of class \code{coords}}
  \item{x}{Data matrix, rows are samples, columns are variables.}
  \item{y}{Transformed data matrix, rows are samples, columns are variables.}
}
\details{
  The function \code{coords()} creates an object representing a coordinate
  transformation.
}
\value{
  An object of class \code{coords}, as a list with the following
  components:
  \item{loadings}{the loadings: each column is one of the new basis vectors}
  \item{pre.sub}{pre-transformation centring}
  \item{pre.div}{pre-transformation scaling}
  \item{post.mul}{post-transformation scaling}
}
\author{
  Jochen Voss <voss@seehuhn.de>
}
\seealso{
  \code{\link{jvPCA}},
  \code{\link{whiten}}
}
\examples{
  pc <- jvPCA(iris[,1:4], scale=TRUE)
  pc
  head(pc$y)
}
