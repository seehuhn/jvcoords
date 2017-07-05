\name{jvmulti-package}
\alias{jvmulti-package}
\docType{package}
\title{
  \packageTitle{jvmulti}
}
\description{
  This package implements an alternative to the \code{prcomp()} function.
  The main advantage over \code{prcomp()} is that the jvPCA package provides
  functions to convert between the original basis and the principal
  component basis.

  See the documentation for \code{jvPCA()} for information on how
  to use this package.
}
\author{
  Jochen Voss <voss@seehuhn.de>
}
\keyword{
  package
  multivariate
}
\seealso{
  \code{\link{jvPCA}}; alternative implementation are
  \code{\link{prcomp}},
  \code{\link{princomp}}
}
\examples{
  pc <- jvPCA(iris[,1:4], scale=TRUE)
  plot(pc$scores[,1], pc$scores[,2], col=iris$Species)
}
