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

standardize <- function(x) {
    x <- as.matrix(x)

    n <- nrow(x)
    p <- ncol(x)
    n.comp <- min(n - 1, p)

    col.mean <- colMeans(x)
    s <- sqrt((colSums(x^2) - col.mean^2*n) / (n-1))
    trans <- coords(diag(1, p, p), pre.sub = col.mean, pre.div = s)
    trans$name <- "standardize"
    trans$y <- t((t(x) - col.mean) / s)

    trans
}
