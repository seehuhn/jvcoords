\name{fromPC}
\alias{fromPC}
\title{Convert data from the PCA space back to the original space}
\description{Convert PCA scores back into the space of the original
  data, by inverting the PCA, inverting the scaling and the adding the
  mean.}
\usage{
fromPC(pca, y)
}
\arguments{
  \item{pca}{A \code{jvPCA} object.}
  \item{y}{Either a vector where the length is the number of principal
    components, or a matrix where the same number of columns equals the
    number of principal components.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{The input \code{y} expressed in the basis of the original data.}
\author{Jochen Voss}
\examples{
  pc <- jvPCA(iris[, 1:4], scale=TRUE)
  head(pc$scores)
  x <- fromPC(pc, pc$scores)
  all.equal(x, as.matrix(iris[, 1:4]), check.attributes = FALSE)
}
