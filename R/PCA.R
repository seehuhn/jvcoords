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

PCA <- function(x, n.comp, scale = FALSE, compute.scores = TRUE) {
  x <- as.matrix(x)

  n <- nrow(x)
  p <- ncol(x)
  max.comp <- min(n - 1, p)
  if (missing(n.comp)) {
    n.comp <- max.comp
  } else {
    stopifnot(n.comp <= max.comp)
  }

  trfm <- coords(p, "PCA")

  col.mean <- colMeans(x)
  xt <- t(x) - col.mean
  trfm <- appendTrfm(trfm, "add", -col.mean)

  if (scale) {
    col.sd <- sqrt((colSums(x^2) - col.mean^2 * n) / (n - 1))
    xt <- xt / col.sd
    trfm <- appendTrfm(trfm, "mult", 1 / col.sd)
  }

  s <- La.svd(xt, nu = n.comp, nv = 0)

  loadings <- s$u
  rownames(loadings) <- colnames(x)
  colnames(loadings) <- paste0("PC", seq.int(ncol(loadings)))
  trfm <- appendTrfm(trfm, "orth", loadings)

  trfm$loadings <- loadings

  if (compute.scores) {
    trfm$y <- crossprod(xt, loadings)
  }

  pc.var <- s$d^2 / (n - 1)
  names(pc.var) <- paste0("PC", seq_along(pc.var))
  trfm$var <- pc.var[1:n.comp]
  trfm$total.var <- sum(pc.var)

  trfm
}
