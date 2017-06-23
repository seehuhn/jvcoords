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
\packageAuthor{jvPCA}

Maintainer: \packageMaintainer{jvPCA}
}
\keyword{
  package
  multivariate
}
\seealso{
  \code{\link{prcomp}}
  \code{\link{toPC}}
  \code{\link{fromPC}}
}
\examples{
  pc <- jvPCA(iris[,1:4], scale=TRUE)
  plot(pc$scores[,1], pc$scores[,2])
}
