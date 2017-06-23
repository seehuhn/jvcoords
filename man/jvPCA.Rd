\name{jvPCA}
\alias{jvPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{function to perform Principal Component Analysis (PCA)}
\description{Performs a principal components analysis on the given data
  matrix and returns the results as an object of class \code{jvPCA}.}
\usage{
  jvPCA(x, n.comp, scale = FALSE, compute.scores = TRUE)
}
\arguments{
  \item{x}{A data matrix, rows correspond to samples, columns to variables.}
  \item{n.comp}{How many principal components to compute.}
  \item{scale}{Whether to standardise the columns before doing PCA.}
  \item{compute.scores}{Whether to compute the scores (i.e. x in the new
    basis).}
}
\details{
  This function performs Principal Component Analysis (PCA) on the
  data.  The mean of each variable is always subtracted before
  the PCA is performed and, if \code{scale} is set, the data
  will also be rescaled to unit variance.
}
\value{
  An object of class \code{jvPCA}, as a list with the following
  components:
  \item{loadings}{the loadings: each column is one of the new basis vectors}
  \item{scores}{if compute.scores==TRUE, this is x expressed in the new basis}
  \item{var}{the variance of the data along each of the new basis vectors}
  \item{total.var}{the total variance of the data}
  \item{mean}{the column means of the original data}
  \item{scale}{the pre-scaling used for each column}
}
\author{Jochen Voss}
\seealso{
  \code{\link{toPC}}
  \code{\link{fromPC}}
}
\examples{
  pc <- jvPCA(iris[,1:4])
}
