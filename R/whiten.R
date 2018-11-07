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

whiten <- function(x, compute.scores = TRUE) {
  x <- as.matrix(x)

  n <- nrow(x)
  p <- ncol(x)
  n.comp <- min(n - 1, p)

  trfm <- coords(p, "whiten")

  col.mean <- colMeans(x)
  xt <- t(x) - col.mean
  trfm <- appendTrfm(trfm, "add", -col.mean)

  s <- La.svd(xt, nu = n.comp, nv = 0)
  eps <- 1e-14
  cols <- which(s$d[seq_len(n.comp)] >= eps)

  loadings <- s$u[, cols, drop = FALSE]
  rownames(loadings) <- colnames(x)[cols]
  colnames(loadings) <- paste0("W", seq.int(ncol(loadings)))
  trfm <- appendTrfm(trfm, "orth", loadings)

  trfm$loadings <- loadings

  inv <- sqrt(n - 1) / s$d[cols]
  names(inv) <- rownames(loadings)
  trfm <- appendTrfm(trfm, "mult", inv)

  if (compute.scores) {
    trfm$y <- t(crossprod(loadings, xt) * inv)
  }

  trfm
}
