\name{toPC}
\alias{toPC}
\title{Convert data from the original space to the PCA space}
\description{Convert data into the principal component basis: this
  subtracts the original mean for each variable, applies the scaling (if
  any) to each variable, and then computes the PCA scores.}
\usage{
  toPC(trans, x)
}
\arguments{
  \item{trans}{A jvPCA object.}
  \item{x}{Either a vector where the length is the number of the
    original variables, or a matrix with the same number of columns as
    the original data.}
}
\details{
  This function implements to formula to compute the PCA scores,
  so that the scores for additional data can be computed.  If all
  components are included in \code{pca} (i.e., if no dimension reduction
  is performed), this operation can be reverted using \code{fromPC()}.

  When applied to the original data \code{x}, the output of
  \code{toPC(pca, x)} equals \code{pca$scores}.
}
\value{The input x expressed in the principal component basis, i.e. the
  PCA scores.}
\author{
  Jochen Voss <voss@seehuhn.de>
}
\seealso{
  \code{\link{jvPCA}},
  \code{\link{fromPC}}
}
\examples{
  pc <- jvPCA(iris[, 1:4], scale=TRUE)
  y <- toPC(pc, iris[, 1:4])
  all.equal(pc$scores, y)
}
