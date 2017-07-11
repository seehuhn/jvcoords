coords <- function(loadings, pre.sub=0, pre.div=1, post.mul=1, post.add=0) {
	# nrow(loadings) == p
	# ncol(loadings) == q
	stopifnot(length(pre.sub) == 1 || length(pre.sub) == nrow(loadings))
	stopifnot(length(pre.div) == 1 || length(pre.div) == nrow(loadings))
	stopifnot(length(post.mul) == 1 || length(post.mul) == ncol(loadings))
	stopifnot(length(post.add) == 1 || length(post.add) == ncol(loadings))
	res <- list(loadings=loadings,
				pre.sub=pre.sub,
				pre.div=pre.div,
				post.add=post.add,
				post.mul=post.mul)
	class(res) <- "coords"
	res
}

toCoords <- function(trans, x) {
	if (is.null(dim(x))) {
		stopifnot(length(x) == nrow(trans$loadings))
		x <- (x - trans$pre.sub) / trans$pre.div
		y <- as.vector(x %*% trans$loadings)
		y * trans$post.mul + trans$post.add
	} else {
		xt <- (t(x) - trans$pre.sub) / trans$pre.div
		yt <- crossprod(trans$loadings, xt)
		t(yt * trans$post.mul + trans$post.add)
	}
}

fromCoords <- function(trans, y) {
	if (is.null(dim(y))) {
		stopifnot(length(x) == ncol(trans$loadings))
		y <- (y - trans$post.add) / trans$post.mul
		x <- as.vector(trans$loadings %*% y)
		x * trans$pre.div + trans$pre.sub
	} else {
		yt <- (t(y) - trans$post.add) / trans$post.mul
		xt <- trans$loadings %*% yt
		t(xt * trans$pre.div + trans$pre.sub)
	}
}

print.coords <- function(x, ...) {
	if (! is.null(x$name)) {
		cat(x$name, ": ", sep="")
	}
	cat("map p =", nrow(x$loadings), "coordinates to q =",
		ncol(x$loadings), "coordinates\n")

	if (! is.null(x$var)) {
		cat("\n")
		info <- rbind(
			`standard deviation`=sqrt(x$var),
			variance=x$var,
			`cum. variance fraction`=cumsum(x$var / x$total.var)
		)
		print(info)
	}

	invisible(x)
}
