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

coords <- function(p, name=NULL, shift=0) {
  res <- list(name=name, p=p, q=p, shift=shift, cmds=list())
  class(res) <- "coords"
  res
}

appendTrfm <- function(trfm, op=c("diag", "orth"), val) {
  op <- match.arg(op)
  if (op == "orth") {
    trfm$q <- ncol(val)
    cmd <- list(op, val)
  } else {
    cmd <- list(op, val)
  }
  trfm$cmds[[length(trfm$cmds) + 1]] <- cmd
  trfm
}

toCoords <- function(trfm, x) {
  if (is.null(dim(x))) {
    # assume that x is a vector
    stopifnot(length(x) == trfm$p)
    x <- x - trfm$shift
    for (cmd in trfm$cmds) {
      val <- cmd[[2]]
      x <- switch(cmd[[1]],
                  diag = x * val,
                  orth = (x %*% val)[1, ])
    }
  } else {
    # x is a matrix
    xt <- t(x) - trfm$shift
    for (cmd in trfm$cmds) {
      val <- cmd[[2]]
      xt <- switch(cmd[[1]],
                   diag = xt * val,
                   orth = crossprod(val, xt))
    }
    x <- t(xt)
  }
  x
}

fromCoords <- function(trfm, y, apply.shift=TRUE) {
  if (is.null(dim(y))) {
    # assume that y is a vector
    stopifnot(length(y) == trfm$q)
    for (cmd in rev(trfm$cmds)) {
      val <- cmd[[2]]
      y <- switch(cmd[[1]],
                  diag = y / val,
                  orth = (val %*% y)[, 1])
    }
    if (apply.shift) {
      y <- y + trfm$shift
    }
  } else {
    # y is a matrix
    yt <- t(y)
    for (cmd in rev(trfm$cmds)) {
      val <- cmd[[2]]
      yt <- switch(cmd[[1]],
                   diag = yt / val,
                   orth = val %*% yt)
    }
    if (apply.shift) {
      yt <- yt + trfm$shift
    }
    y <- t(yt)
  }
  y
}

print.coords <- function(x, ...) {
  if (!is.null(x$name)) {
    cat(x$name, ": ", sep = "")
  }
  cat("mapping p =", x$p, "coordinates",
      "to q =", x$q, "coordinates\n")

  if (!is.null(x$var)) {
    cat("\n")
    info <- rbind(`standard deviation` = sqrt(x$var),
                  variance = x$var,
                  `cum. variance fraction` = cumsum(x$var / x$total.var))
    print(info)
  }

  invisible(x)
}
