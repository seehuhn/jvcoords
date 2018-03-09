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
