whiten <- function(x) {
	x <- as.matrix(x)

	n <- nrow(x)
	p <- ncol(x)
	stopifnot(n > 1 && p >= 1)
	n.comp <- min(n-1, p)

	col.mean <- apply(x, 2, mean)

	z <- t(t(x) - col.mean)
	s <- La.svd(z, nu=0, nv=n.comp)
	stopifnot(s$d[1] >= .Machine$double.eps)
	if (s$d[n.comp] < .Machine$double.eps) {
		while (s$d[n.comp] < .Machine$double.eps) {
			n.comp <- n.comp - 1
		}
		# protect against s$vt turning into a vector when n.comp == 1:
		s$vt <- s$vt[seq_len(n.comp), , drop=FALSE]
	}
	colnames(s$vt) <- colnames(x)
	rownames(s$vt) <- paste0('W', seq_len(n.comp))

	inv <- sqrt(n-1) / s$d[seq_len(n.comp)]
	names(inv) <- rownames(s$vt)

	w <- t(tcrossprod(s$vt, z) * inv)

	trans <- list()
	trans$mean <- col.mean
	trans$vt <- s$vt
	trans$inv <- inv
	trans$w <- w

	trans
}

toWhite <- function(trans, x) {
	if (is.null(dim(x))) {
		stopifnot(length(x) == length(trans$mean))
		(trans$vt %*% (x - trans$mean))[,1] * trans$inv
	} else {
		t(trans$vt %*% (t(x) - trans$mean) * trans$inv)
	}
}

fromWhite <- function(trans, w) {
	if (is.null(dim(w))) {
		stopifnot(length(w) == length(trans$inv))
		crossprod(trans$vt, w / trans$inv)[,1] + trans$mean
	} else {
		t(t(w %*% (trans$vt / trans$inv)) + trans$mean)
	}
}
