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

	col.mean <- apply(x, 2, mean)
	if (scale) {
		col.sd <- apply(x, 2, sd)
	} else {
		col.sd <- 1
	}
	z <- t((t(x) - col.mean) / col.sd)
	s <- La.svd(z, nu=0, nv=n.comp)

	pca <- list()

	pca$mean <- col.mean
	pca$scale <- col.sd

	pc.var <- s$d^2 / (n-1)
	pca$var <- pc.var[1:n.comp]
	names(pca$var) <- paste0('PC', 1:length(pca$var))
	pca$total.var <- sum(pc.var)

	loadings <- t(s$vt)
	rownames(loadings) <- names(col.mean)
	colnames(loadings) <- paste0('PC', 1:ncol(loadings))
	pca$loadings <- loadings

	if (compute.scores) {
		pca$scores <- z %*% pca$loadings
	}

	class(pca) <- "jvPCA"
	pca
}

print.jvPCA <- function(x, ...) {
	cat("cumulative fraction of variance explained:\n")
	print(cumsum(x$var) / x$total.var)
	cat("\nloadings:\n")
	print(x$loadings)
}

toPC <- function(trans, x) {
	if (is.null(dim(x))) {
		stopifnot(length(x) == length(trans$mean))
		(t(trans$loadings) %*% ((x - trans$mean) / trans$scale))[,1]
	} else {
		t((t(x) - trans$mean) / trans$scale) %*% trans$loadings
	}
}

fromPC <- function(trans, y) {
	if (is.null(dim(y))) {
		(trans$loadings %*% y)[,1] * trans$scale + trans$mean
	} else {
		# y = t((t(x) - trans$mean) / trans$scale) %*% trans$loadings
		t((trans$loadings %*% t(y)) * trans$scale + trans$mean)
	}
}
