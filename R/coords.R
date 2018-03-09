coords <- function(loadings=NULL, pre.sub = 0, pre.div = 1, post.mul = 1) {
    if (! is.null(loadings)) {
        p <- nrow(loadings)
        q <- ncol(loadings)
    } else {
        if (length(pre.sub) > 1) {
            p <- length(pre.sub)
        } else if (length(pre.div) > 1) {
            p <- length(pre.div)
        } else if (length(post.mul) > 1) {
            p <- length(post.mul)
        } else {
            stop("cannot infer p and q")
        }
        q <- p
    }
    stopifnot(length(pre.sub) == 1 || length(pre.sub) == p)
    stopifnot(length(pre.div) == 1 || length(pre.div) == p)
    stopifnot(length(post.mul) == 1 || length(post.mul) == q)
    res <- list(loadings = loadings, p = p, q = q,
        pre.sub = pre.sub, pre.div = pre.div, post.mul = post.mul)
    class(res) <- "coords"
    res
}

toCoords <- function(trans, x) {
    if (is.null(dim(x))) {
        stopifnot(length(x) == trans$p)
        x <- (x - trans$pre.sub) / trans$pre.div
        if (is.null(trans$loadings)) {
            y <- x
        } else {
            y <- (x %*% trans$loadings)[1,]
        }
        y * trans$post.mul
    } else {
        xt <- (t(x) - trans$pre.sub) / trans$pre.div
        if (is.null(trans$loadings)) {
            yt <- xt
        } else {
            yt <- crossprod(trans$loadings, xt)
        }
        t(yt * trans$post.mul)
    }
}

fromCoords <- function(trans, y) {
    if (is.null(dim(y))) {
        stopifnot(length(y) == trans$q)
        y <- y / trans$post.mul
        if (is.null(trans$loadings)) {
            x <- y
        } else {
            x <- (trans$loadings %*% y)[, 1]
        }
        x * trans$pre.div + trans$pre.sub
    } else {
        yt <- t(y) / trans$post.mul
        if (is.null(trans$loadings)) {
            xt <- yt
        } else {
            xt <- trans$loadings %*% yt
        }
        t(xt * trans$pre.div + trans$pre.sub)
    }
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
