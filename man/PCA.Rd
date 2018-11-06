\name{PCA}
\alias{PCA}
\title{perform Principal Component Analysis (PCA)}
\description{Performs a principal components analysis on the given data
  matrix and returns the results as an object of class \code{coords}.}
\usage{
  PCA(x, n.comp, scale = FALSE, compute.scores = TRUE)
}
\arguments{
  \item{x}{A data matrix, rows are observations, columns are variables.}
  \item{n.comp}{How many principal components to compute.}
  \item{scale}{Whether to standardize the columns before doing PCA.}
  \item{compute.scores}{Whether to compute the scores (i.e. x in the new
    basis).}
}
\details{
  This function performs Principal Component Analysis (PCA) on the
  data.  Variables are always centred before
  the PCA is performed and, if \code{scale} is set, the variables
  will also be rescaled to unit variance.

  If \code{compute.scores} is set to \code{FALSE}, only the information
  required for the \code{toPC()} and \code{fromPC()} to work is stored
  in the returned \code{coords} object; otherwise the scores will
  be stored in the \code{$y} field of the \code{coords} object.

  The \code{PCA()} function is an alternative to
  the \code{prcomp()} command from the standard library.
  The main advantage of \code{PCA()} is that the \code{coords}
  class provides
  functions to convert between the original basis and the principal
  component basis.
}
\value{
  An object of class \code{coords}, as a list with the following
  components:
  \item{loadings}{the loadings: each column is one of the new basis vectors}
  \item{y}{if compute.scores==TRUE, this is x expressed in the new basis}
  \item{var}{the variance of the data along each of the new basis vectors}
  \item{total.var}{the total variance of the data}
  \item{pre.sub}{the column means of the original data}
  \item{pre.div}{the standard deviations of the original variables}
}
\author{
  Jochen Voss <voss@seehuhn.de>
}
\seealso{
  \code{\link{coords}};
  alternative implementations: \code{\link{prcomp}}, \code{\link{princomp}}
}
\examples{
  pc <- PCA(iris[,1:4], scale=TRUE)
  pc
  plot(pc$y[,1], pc$y[,2], col=iris$Species)
}
