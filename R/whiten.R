whiten <- function(x, compute.scores=TRUE) {
	x <- as.matrix(x)

	n <- nrow(x)
	p <- ncol(x)
	n.comp <- min(n-1, p)

	col.mean <- colMeans(x)

	zt <- t(x) - col.mean
	s <- La.svd(zt, nu=n.comp, nv=0)
	eps <- 1e-14
	cols <- which(s$d[seq_len(n.comp)] >= eps)

	loadings <- s$u[, cols, drop=FALSE]
	rownames(loadings) <- colnames(x)[cols]
	colnames(loadings) <- paste0('W', 1:ncol(loadings))

	inv <- sqrt(n-1) / s$d[cols]
	names(inv) <- rownames(loadings)

	trans <- coords(loadings, pre.sub = col.mean, post.mul = inv)
	trans$name <- "whiten"

	if (compute.scores) {
		trans$scores <- crossprod(zt, trans$loadings)
	}

	trans
}
