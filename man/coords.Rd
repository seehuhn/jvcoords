\encoding{UTF-8}
\name{coords}
\alias{coords}
\alias{appendTrfm}
\alias{fromCoords}
\alias{toCoords}
\title{An S3 class to represent affine coordinate transforms}
\description{Perform affine coordinate transformations.}
\usage{
  coords(p, name = NULL, shift = 0)
  appendTrfm(trfm, op = c("diag", "orth"), val)
  toCoords(trfm, x)
  fromCoords(trfm, y, apply.shift = TRUE)
}
\arguments{
  \item{p}{The dimension of the original data}
  \item{name}{A short name for the coordinate transformation (optional)}
  \item{shift}{A value subtract from the data as the first step of the coordinate transformation}
  \item{trfm}{an object of class \code{coords}}
  \item{op}{the type of transformation to append}
  \item{val}{data for the transformation to append}
  \item{x}{Data matrix, rows are observations, columns are variables.}
  \item{y}{Transformed data matrix, rows are observations, columns are variables.}
  \item{apply.shift}{Whether to apply the final shift of coordinates.}
}
\details{
  The function \code{coords()} creates a new object representing an affine
  coordinate transformation.  Initially, the object represents a shift of
  vectors vectors of length \code{p} by the amount \code{shift}.  The function
  \code{appendTrfm()} can be used to modify the transformation as needed.
  The optional argument \code{name}, if set, is used when printing
  objects of class \code{coords}.

  The function \code{appendTrfm()} modifies the transformation \code{trfm}
  by concatenating it with a new affine transformation and returning the
  result.  The arguments \code{op} and \code{val} specify which affine
  transformation to append.  There are three choices for \code{op}:
  \itemize{
  \item \code{diag} denotes multiplication by a diagonal matrix: an input
    vector \code{x} is mapped to the output \code{x * val}.  The scaling factor
    \code{val} can either be a vector of length \code{trfm$q} (for element-wise
    scaling), or a number.
  \item \code{orth} denotes multiplication with an orthogonal matrix.
    \code{val} must be a matrix with orthogonal columns (not necessarily
    square) and \code{trfm$q} rows.  An input vector \code{x} is mapped
    to the output \code{x \%*\% orth}.
  }
  The new transformation is applied after any other transformations
  already associated with \code{trfm}.

  The function \code{toCoords()} applies the affine transformation \code{trfm}
  to the data \code{x}.  The data \code{x} must either be a vector of
  length \code{trfm$p}, in which case the result is a vector of
  length \code{trfm$q}, or a matrix with \code{trfm$p} columns, in which case
  the transformation is applied to each row of the matrix separately.

  The function \code{fromCoords()} implements the inverse transform
  to \code{toCoords()}.  The output always satisfies
  \code{toCoords(trfm, fromCoords(trfm, y)) == y}.  If
  \code{trfm$p == trfm$q}, \emph{i.e.} if the transformation is bijective,
  the \code{fromCoords(trfm, toCoords(trfm, x)) == x} also holds.
  The argument \code{apply.shift} can be set to false to apply only
  the linear part of the (inverse) transformation, leaving out the
  final shift.
}
\value{
  An object of class \code{coords}, as a list with the following
  components:
  \item{p}{the number of variables in the original data set}
  \item{q}{the number of variables in the transformed data set}
  \item{shift}{the affine part of the transformation}
  \item{name}{the name of the transformation}
  \item{cmds}{a representation of the transformation (internal use only)}
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
  toCoords(pc, c(5, 3, 4, 1))
  fromCoords(pc, c(1, 0, 0))
}
