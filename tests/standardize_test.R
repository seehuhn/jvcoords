#! /usr/bin/env Rscript

library(jvcoords)

set.seed(1)
X1 <- matrix(rnorm(990), nrow = 33, ncol = 3)
X2 <- matrix(rnorm(990), nrow = 3, ncol = 33)
X3 <- matrix(rnorm(999), nrow = 999, ncol = 1)
X4 <- matrix(1:6, 3, 2)

for (x in list(X1, X2, X3, X4)) {
    s <- standardize(x)

    stopifnot(max(abs(colMeans(s$y))) < 1e-10)
    stopifnot(max(abs(apply(s$y, 2, var) - 1)) < 1e-10)

    stopifnot(max(abs(toCoords(s, x) - s$y)) < 1e-10)
    cat(dim(x), "to mat OK\n")

    stopifnot(max(abs(toCoords(s, x[1, ]) - s$y[1, ])) < 1e-10)
    cat(dim(x), "to vec OK\n")

    stopifnot(max(abs(fromCoords(s, s$y) - x)) < 1e-10)
    cat(dim(x), "from mat OK\n")

    stopifnot(max(abs(fromCoords(s, s$y[1, ]) - x[1, ])) < 1e-10)
    cat(dim(x), "from vec OK\n")

    Y <- matrix(rnorm(100 * s$q), 100, s$q)

    stopifnot(max(abs(Y - toCoords(s, fromCoords(s, Y)))) < 1e-10)
    cat(dim(x), "from&to mat OK\n")

    stopifnot(max(abs(Y[1, ] - toCoords(s, fromCoords(s, Y[1, ])))) < 1e-10)
    cat(dim(x), "from&to vec OK\n")
}
