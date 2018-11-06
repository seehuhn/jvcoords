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

library(jvcoords)

set.seed(1)
X1 <- matrix(rnorm(990), nrow = 33, ncol = 3)
X2 <- matrix(rnorm(990), nrow = 3, ncol = 33)
X3 <- matrix(rnorm(999), nrow = 999, ncol = 1)
X4 <- matrix(1:6, 3, 2)

for (x in list(X1, X2, X3, X4)) {
    for (scale in c(TRUE, FALSE)) {
        pc <- PCA(x, scale = scale)
        stopifnot(max(abs(colMeans(pc$y))) < 1e-14)
        stopifnot(which.max(apply(pc$y, 2, var)) == 1)
        stopifnot(isTRUE(all.equal(ToCoords(pc, x), pc$y)))
        stopifnot(isTRUE(all.equal(x, FromCoords(pc, pc$y))))
    }
}
