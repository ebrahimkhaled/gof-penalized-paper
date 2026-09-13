# Verify my EDGE statistic reproduces ebrahim.gof::edge.gof exactly on an MLE fit.
suppressPackageStartupMessages({library(glmnet); library(ebrahim.gof)})
set.seed(1); G <- 10
ridge_fit <- function(X1, y, lam, D, tol=1e-9, maxit=60) {
  beta <- rep(0, ncol(X1))
  for (it in 1:maxit) {
    eta <- drop(X1 %*% beta); pi <- 1/(1+exp(-eta))
    w <- pmax(pi*(1-pi), 1e-8); z <- eta + (y-pi)/w
    bn <- drop(solve(crossprod(X1, w*X1) + lam*D, crossprod(X1, w*z)))
    if (max(abs(bn-beta)) < tol) { beta <- bn; break }
    beta <- bn
  }
  eta <- drop(X1 %*% beta); list(beta=beta, pi=1/(1+exp(-eta)))
}
edge_stat <- function(X1, y, fit) {
  pi <- pmin(pmax(fit$pi,1e-6),1-1e-6); w <- pmax(pi*(1-pi),1e-8)
  g <- pmin(ceiling(rank(pi, ties.method="first")/(length(y)/G)), G)
  idx <- split(seq_along(y), g)
  Vg <- vapply(idx,function(I) sum(w[I]),0.0); og <- vapply(idx,function(I) sum(y[I]),0.0)
  eg <- vapply(idx,function(I) sum(pi[I]),0.0); pbar <- vapply(idx,function(I) mean(pi[I]),0.0)
  r <- (og-eg)/sqrt(Vg); Z <- as.matrix(stats::poly(pbar,3))
  Zr <- crossprod(Z,r); as.numeric(t(Zr) %*% solve(crossprod(Z)) %*% Zr)
}
n<-500; p<-5; b<-c(0.5,-0.4,0.3,-0.3,0.2); D<-diag(c(0,rep(1,p)))
for (k in 1:3) {
  X <- matrix(rnorm(n*p),n,p); X1 <- cbind(1,X)
  y <- rbinom(n,1,1/(1+exp(-drop(X%*%b))))
  g <- glm(y ~ X, family=binomial)
  pk <- edge.gof(g)$Test_Statistic
  mine <- edge_stat(X1, y, list(pi=as.numeric(fitted(g))))
  cat(sprintf("rep %d  package=%.8f  mine=%.8f  diff=%.2e\n", k, pk, mine, abs(pk-mine)))
}
