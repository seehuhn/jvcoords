jvPCA <- function(x, n.comp, scale=FALSE, compute.scores=TRUE) {
	x <- as.matrix(x)

	n <- nrow(x)
	p <- ncol(x)
	max.comp <- min(n-1, p)
	if (missing(n.comp)) {
		n.comp = max.comp
	} else {
		stopifnot(n.comp <= max.comp)
	}

	col.mean <- colMeans(x)
	if (scale) {
		col.sd <- sqrt((colSums(x^2) - col.mean^2*n) / (n-1))
	} else {
		col.sd <- 1
	}

	zt <- (t(x) - col.mean) / col.sd
	s <- La.svd(zt, nu=n.comp, nv=0)

	loadings <- s$u
	rownames(loadings) <- colnames(x)
	colnames(loadings) <- paste0('PC', 1:ncol(loadings))
	trans <- coords(loadings, pre.sub = col.mean, pre.div = col.sd)
	trans$name <- "PCA"

	pc.var <- s$d^2 / (n-1)
	trans$var <- pc.var[1:n.comp]
	names(trans$var) <- paste0('PC', 1:length(trans$var))
	trans$total.var <- sum(pc.var)

	if (compute.scores) {
		trans$scores <- crossprod(zt, trans$loadings)
	}

	trans
}
