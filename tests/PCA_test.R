#! /usr/bin/env Rscript

library(jvmulti)

set.seed(1)
X1 <- matrix(rnorm(990), nrow = 33, ncol = 3)
X2 <- matrix(rnorm(990), nrow = 3, ncol = 33)
X3 <- matrix(rnorm(999), nrow = 999, ncol = 1)
X4 <- matrix(1:6, 3, 2)

for (x in list(X1, X2, X3, X4)) {
    for (scale in c(TRUE, FALSE)) {
        pc <- jvPCA(x, scale = scale)
        stopifnot(max(abs(colMeans(pc$y))) < 1e-14)
        stopifnot(which.max(apply(pc$y, 2, var)) == 1)
        stopifnot(isTRUE(all.equal(toCoords(pc, x), pc$y)))
        stopifnot(isTRUE(all.equal(x, fromCoords(pc, pc$y))))
    }
}
