\name{toPC}
\alias{toPC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert data from the original space to the PCA space}
\description{Convert data into the principal component basis: this
  subtracts the original mean for each variable, applies the scaling (if
  any) to each variable and then computes the PCA scores.}
\usage{
toPC(pca, x)
}
\arguments{
  \item{pca}{A jvPCA object.}
  \item{x}{Either a vector where the length is the number of the
    original variables, or a matrix with the same number of columns as
    the original data.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{The input x expressed in the principal component basis, i.e. the
  PCA scores.}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (pca, x)
{
    if (is.null(dim(x))) {
	stopifnot(length(x) == length(pca$mean))
	(t(pca$loadings) \%*\% ((x - pca$mean)/pca$scale))[, 1]
    }
    else {
	t((t(x) - pca$mean)/pca$scale) \%*\% pca$loadings
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
