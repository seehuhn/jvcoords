\name{jvPCA-package}
\alias{jvPCA-package}
\docType{package}
\title{
  \packageTitle{jvPCA}
}
\description{
  \packageDescription{jvPCA}
}
\details{
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
  plot(pc$scores[,1], pc$scores[,2])
}
