#! /usr/bin/env Rscript

# jvcoords - implement various coordinate transforms (e.g. PCA, whitening).
# https://github.com/seehuhn/jvcoords
#
# Copyright (C) 2018  Jochen Voss <voss@seehuhn.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

library("jvcoords")

# test addition
x <- rnorm(3)
y <- rnorm(3)
m <- coords(p = 3)
m <- AppendTrfm(m, "add", y)
z <- ToCoords(m, x)
stopifnot(isTRUE(all.equal(x + y, z)))

# test multiplication
m <- coords(p = 3)
m <- AppendTrfm(m, "mult", y)
z <- ToCoords(m, x)
stopifnot(isTRUE(all.equal(x * y, z)))

# test orthogonal transformation
m <- coords(p = 3)
U <- matrix(c(1, 0, 0, 0, 1, 0), 3, 2)
m <- AppendTrfm(m, "orth", U)
stopifnot(m$p == 3 && m$q == 2)
z <- ToCoords(m, x)
stopifnot(isTRUE(all.equal(x[1:2], z)))

# end-to-end tests
mu <- c(1, 2, 3)
sigma <- c(4, 5, 6)
A <- matrix(c(1, 0, 0, 0, 0, 1, 0, 1, 0), 3, 3, byrow = T)
m1 <- coords(p = 3)
m1 <- AppendTrfm(m1, "add", -mu)
m1 <- AppendTrfm(m1, "mult", 1/sigma)
m1 <- AppendTrfm(m1, "orth", A)

x <- c(7, 8, 9)
y1 <- ToCoords(m1, x)
y2 <- as.numeric(((x - mu)/sigma) %*% A)
stopifnot(isTRUE(all.equal(y1, y2)))

A <- diag(1, 100, 10)
m2 <- coords(p = nrow(A))
m2 <- AppendTrfm(m2, "add", rnorm(nrow(A)))
m2 <- AppendTrfm(m2, "mult", rexp(nrow(A)))
m2 <- AppendTrfm(m2, "orth", A)
m2 <- AppendTrfm(m2, "mult", rexp(ncol(A)))

A <- La.svd(matrix(rnorm(400), 20, 20))$u
m3 <- coords(p = nrow(A))
m3 <- AppendTrfm(m3, "add", rnorm(nrow(A)))
m3 <- AppendTrfm(m3, "mult", rexp(nrow(A)))
m3 <- AppendTrfm(m3, "orth", A)
m3 <- AppendTrfm(m3, "mult", rexp(ncol(A)))

mm <- list(m1, m2, m3)
for (i in seq_along(mm)) {
    m <- mm[[i]]
    p <- m$p
    q <- m$q

    n <- 100
    y <- matrix(rnorm(n * q), n, q)
    x <- FromCoords(m, y)
    x1 <- FromCoords(m, y[1, ])
    stopifnot(isTRUE(all.equal(x1, x[1, ])))
    z <- ToCoords(m, x)
    z1 <- ToCoords(m, x[1, ])
    stopifnot(isTRUE(all.equal(z1, z[1, ])))
    stopifnot(isTRUE(all.equal(y, z)))
}
