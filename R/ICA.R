# Estimate the entropy of a continuous distribution from data.
# Based on the estimate (13) in "Nonparametric entropy estimation:
# an overview" by J. Beirlant et al.
# TODO: is this the best way to estimate entropy for our purpose?
entropy <- function(x, m=3) {
  x <- sort(x)
  n <- length(x)
  d <- x[(m+1):n] - x[1:(n-m)]
  # If data is rounded/integer-valued, the spacing between two samples
  # may be 0 and the estimate for the entropy is negative infinity.
  # While this is the mathematically correct answer, it confuses
  # the optimisation algorithm, so we add 1e-6 before taking the log.
  sum(log(n * d / m + 1e-6)) / n - digamma(m) + log(m)
}

jvICA <- function(X, n.comp, maxit=5000) {
  trans <- whiten(X)
  z <- trans$w
  q <- ncol(z)

  if (missing(n.comp)) {
    n.comp <- q
  } else {
    stopifnot(n.comp > 0 && n.comp <= q)
  }

  # IC is an orthogonal matrix.  The following loop will modify IC to produce
  # the independent components in the columns of IC one by one, while keeping
  # the matrix orthogonal.
  IC <- diag(q)
  # If the remaining search space is one-dimensional (i.e. if k=q), no search
  # is required and the corresponding column of IC already contains a normed
  # vector in this space.  Thus we need to find at most q-1 directions.
  for (k in 1:min(n.comp, q-1)) {
  	# Find the direction in the space spanned by IC[,k], ..., IC[,q]
  	# which minimises entropy.
  	cat("optimising direction", k, "out of", n.comp, "\n")

    r <- q - k + 1 # the dimension of the search space

    # step 1: Try random directions to
    # find a good starting point for the optimisation.
    n.trials <- min(10000, 100 * 2^(r-1))
    trials <- matrix(rnorm(n.trials*r), n.trials, r)
    trials <- trials / sqrt(rowSums(trials^2))
    trials.orig.space <- IC %*% rbind(matrix(0, k-1, n.trials), t(trials))
    best.value <- Inf
    best.dir <- NULL
    # TODO: is it possible/faster to get rid of this loop?
    for (i in 1:n.trials) {
    	value <- entropy(z %*% trials.orig.space[,i])
    	if (value < best.value) {
    		best.value <- value
    		best.dir <- trials[i,]
    	}
    }

    # step 2: use local optimisation to find the best solution in the
    # vincinity of best.dir
    res <- optim(best.dir, function(dir) {
    	dir <- dir / sqrt(sum(dir^2))
    	dir.orig.space <- IC %*% c(rep(0, k-1), dir)
    	entropy((z %*% dir.orig.space)[,1])
    }, control=list(maxit=maxit))
    if (res$convergence == 1) {
    	warning("optimisation did not converge, consider increasing maxit")
    } else if (res$convergence != 0) {
    	warning("optimisation did not converge (error ", res$convergence, ")")
    }
    stopifnot(res$value <= best.value)
    print(res$value)
    best.dir <- res$par / sqrt(sum(res$par^2))

    # Use a Householder reflection which maps e1 to best.dir to update IC.
    e1 <- c(1, rep(0, r-1))
    v <- best.dir - e1
    v <- v / sqrt(sum(v^2))
    P <- diag(r) - 2 * tcrossprod(v)
    IC[,k:q] <- IC[,k:q,drop=FALSE] %*% P
  }
  IC <- IC[, seq_len(n.comp), drop=FALSE]
  colnames(IC) <- paste0('IC', seq_len(n.comp))

  res <- list()
  res$trans <- trans
  res$IC <- IC
  res$x <- z %*% IC
  res
}

#s <- jvICA(iris[,1:4])
#pairs(s$x, col=iris$Species)
