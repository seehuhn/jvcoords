#! /usr/bin/env Rscript

library("jvmulti")

n <- 1000

mu <- 0
sigma <- 1
target <- 0.5 * log(2 * pi * exp(1) * sigma^2)

K <- 10000
e <- numeric(K)
for (k in 1:K) {
    x <- rnorm(n, mu, sigma)
    e[k] <- entropy(x) - target
}
print(mean(e[k]^2))
