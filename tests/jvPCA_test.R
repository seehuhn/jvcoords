#! /usr/bin/env Rscript

library(jvPCA)

set.seed(1)
X1 <- matrix(rnorm(990), nrow=33, ncol=3)
X2 <- matrix(rnorm(990), nrow=3, ncol=33)
X3 <- matrix(rnorm(999), nrow=999, ncol=1)
X4 <- matrix(1:6, 3, 2)

for (x in list(X1, X2, X3, X4)) {
  pc <- jvPCA(x)
  stopifnot(max(abs(colMeans(pc$scores))) < 1e-10)
  stopifnot(isTRUE(all.equal(toPC(pc, x), pc$scores)))
  stopifnot(isTRUE(all.equal(x, fromPC(pc, pc$scores))))
}
