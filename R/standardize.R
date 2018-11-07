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

standardize <- function(x, compute.scores = TRUE) {
  x <- as.matrix(x)

  n <- nrow(x)
  p <- ncol(x)

  trfm <- coords(p, "standardize")

  col.mean <- colMeans(x)
  trfm <- appendTrfm(trfm, "add", -col.mean)

  col.sd <- sqrt((colSums(x^2) - col.mean^2 * n) / (n - 1))
  trfm <- appendTrfm(trfm, "mult", 1 / col.sd)

  if (compute.scores) {
    trfm$y <- t((t(x) - col.mean) / col.sd)
  }

  trfm
}
