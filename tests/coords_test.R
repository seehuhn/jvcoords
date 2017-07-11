#! /usr/bin/env Rscript

library("jvmulti")

mu <- c(1, 2, 3)
sigma <- c(4, 5, 6)
A <- matrix(c(1, 0, 0, 0, 0, 1, 0, 1, 0), 3, 3, byrow=T)
m1 <- coords(loadings=A, pre.sub = mu, pre.div = sigma)

x <- c(7, 8, 9)
y1 <- toCoords(m1, x)
y2 <- as.numeric(((x - mu) / sigma) %*% A)
stopifnot(isTRUE(all.equal(y1, y2)))

A <- matrix(0, 100, 10)
diag(A) <- 1
m2 <- coords(loadings = A,
			 pre.sub = rnorm(nrow(A)),
			 pre.div = rexp(nrow(A)),
			 post.mul = rexp(ncol(A)),
			 post.add = rnorm(ncol(A)))

A <- La.svd(matrix(rnorm(400), 20, 20))$u
m3 <- coords(loadings = A,
			 pre.sub = rnorm(nrow(A)),
			 pre.div = rexp(nrow(A)),
			 post.mul = rexp(ncol(A)),
			 post.add = rnorm(ncol(A)))

mm <- list(m1, m2, m3)
for (i in seq_along(mm)) {
	m <- mm[[i]]
	p <- nrow(m$loadings)
	q <- ncol(m$loadings)

	n <- 100
	y <- matrix(rnorm(n*q), n, q)
	x <- fromCoords(m, y)
	x1 <- fromCoords(m, y[1,])
	stopifnot(isTRUE(all.equal(x1, x[1,])))
	z <- toCoords(m, x)
	z1 <- toCoords(m, x[1,])
	stopifnot(isTRUE(all.equal(z1, z[1,])))
	stopifnot(isTRUE(all.equal(y, z)))
}
