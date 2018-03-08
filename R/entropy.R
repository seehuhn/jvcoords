# Estimate the entropy of a continuous distribution from data.  Based on
# the estimate (13) in 'Nonparametric entropy estimation: an overview' by
# J. Beirlant et al., 2001.
entropy <- function(x, m) {
    x <- sort(x)
    n <- length(x)
    if (missing(m)) {
        m <- as.integer(min(3 + n^(1/3), n/2))
    }
    d <- x[(m + 1):n] - x[1:(n - m)]

    # If data is rounded/integer-valued, the spacing between two
    # samples may be 0 and the estimate for the entropy is negative
    # infinity.  While this is the mathematically correct answer, it
    # confuses the optimisation algorithm, so we add 1e-6 before
    # taking the log.
    sum(log(n * d/m + 1e-06)) / n - digamma(m) + log(m)
}
