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

coords <- function(p, name=NULL) {
    res <- list(name=name, p=p, q=p, y=NULL, cmds=list())
    class(res) <- "coords"
    res
}

AppendTrfm <- function(trfm, op=c("add", "mult", "orth"), val) {
    op <- match.arg(op)
    if (op == "orth") {
        trfm$q <- ncol(val)
        cmd <- list(op, val)
    } else {
        cmd <- list(op, val)
    }
    trfm$cmds[[length(trfm$cmds)+1]] <- cmd
    trfm
}

ToCoords <- function(trfm, x) {
    if (is.null(dim(x))) {
        stopifnot(length(x) == trfm$p)
        for (cmd in trfm$cmds) {
            val <- cmd[[2]]
            x <- switch(cmd[[1]],
                        add = x + val,
                        mult = x * val,
                        orth = (x %*% val)[1, ])
        }
    } else {
        xt <- t(x)
        for (cmd in trfm$cmds) {
            val <- cmd[[2]]
            xt <- switch(cmd[[1]],
                         add = xt + val,
                         mult = xt * val,
                         orth = crossprod(val, xt))
        }
        x <- t(xt)
    }
    x
}

FromCoords <- function(trfm, y) {
    if (is.null(dim(y))) {
        stopifnot(length(y) == trfm$q)
        for (cmd in rev(trfm$cmds)) {
            val <- cmd[[2]]
            y <- switch(cmd[[1]],
                        add = y - val,
                        mult = y / val,
                        orth = (val %*% y)[, 1])
        }
    } else {
        yt <- t(y)
        for (cmd in rev(trfm$cmds)) {
            val <- cmd[[2]]
            yt <- switch(cmd[[1]],
                         add = yt - val,
                         mult = yt / val,
                         orth = val %*% yt)
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
                      `cum. variance fraction` = cumsum(x$var/x$total.var))
        print(info)
    }

    invisible(x)
}
