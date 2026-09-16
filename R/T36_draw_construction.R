# T3.6 -- a referee-style challenge to T33: its tau = 0 column does not reproduce the glaucoma table
# cell for cell, and at lambda.min, G = 10 it sits above the eight-stream range that T35 measured
# with rbinom draws (0.076 and 0.166 against [0.036, 0.066] and [0.096, 0.136]). The two studies
# differ in how the bootstrap response is drawn -- T33 compares one shared n x NB matrix of uniforms
# against the generator, T35 and the published analysis call rbinom -- so the question is whether
# that construction shifts the p-value or whether 0.166 was simply an unlucky draw.
# Six seeds of T33's own construction on the same cell. Output: T36_draw_construction.csv
suppressPackageStartupMessages({library(glmnet); library(TH.data)})
data("GlaucomaM", package = "TH.data")
Xr <- as.matrix(GlaucomaM[, setdiff(names(GlaucomaM), "Class")])
y  <- as.numeric(GlaucomaM$Class == "glaucoma")
X  <- scale(Xr); X <- X[, apply(X, 2, function(c) all(is.finite(c)) & sd(c) > 0), drop = FALSE]
n <- nrow(X); p <- ncol(X); X1 <- cbind(1, X); D <- diag(c(0, rep(1, p)))
G <- 10; NB <- 499

ridge_fit <- function(X1, y, lam, D, tol = 1e-9, maxit = 200) {
  beta <- rep(0, ncol(X1))
  for (it in 1:maxit) { eta <- drop(X1 %*% beta); pi <- 1/(1 + exp(-eta))
    w <- pmax(pi*(1-pi), 1e-8); z <- eta + (y - pi)/w
    bn <- drop(solve(crossprod(X1, w*X1) + lam*D, crossprod(X1, w*z)))
    if (max(abs(bn - beta)) < tol) { beta <- bn; break }; beta <- bn }
  list(beta = beta, pi = 1/(1 + exp(-drop(X1 %*% beta))))
}
pieces <- function(X1, y, fit, lam, D, G) {
  pi <- pmin(pmax(fit$pi, 1e-6), 1-1e-6); w <- pmax(pi*(1-pi), 1e-8)
  g <- pmin(ceiling(rank(pi, ties.method = "first")/(length(y)/G)), G)
  idx <- split(seq_along(y), g); Vg <- vapply(idx, function(I) sum(w[I]), 0.0)
  r <- (vapply(idx, function(I) sum(y[I]), 0.0) - vapply(idx, function(I) sum(pi[I]), 0.0))/sqrt(Vg)
  U <- t(vapply(idx, function(I) colSums(w[I]*X1[I, , drop = FALSE]), numeric(ncol(X1))))/sqrt(Vg)
  Fm <- crossprod(X1, w*X1); K <- lam*D
  bt <- fit$beta + drop(solve(Fm, K %*% fit$beta)); mu <- drop(U %*% solve(Fm + K, K %*% bt))
  pbar <- vapply(idx, function(I) mean(pi[I]), 0.0); v <- r - mu
  Z <- tryCatch(as.matrix(stats::poly(pbar, 3)), error = function(e) NULL)
  list(Sc = sum(v^2), bt = bt,
       Sce = if (is.null(Z)) NA else { zr <- crossprod(Z, v)
                                      as.numeric(t(zr) %*% solve(crossprod(Z)) %*% zr) })
}
# T33's construction: one n x NB uniform matrix, compared against the generator
run_unif <- function(lam, seed) {
  ft <- ridge_fit(X1, y, lam, D); st <- pieces(X1, y, ft, lam, D, G)
  gen <- pmin(pmax(as.numeric(1/(1 + exp(-X1 %*% st$bt))), 1e-6), 1 - 1e-6)
  set.seed(seed); Umat <- matrix(runif(n*NB), n, NB)
  Sd <- Se <- numeric(NB)
  for (b in 1:NB) { ys <- as.numeric(Umat[, b] < gen)
    s2 <- pieces(X1, ys, ridge_fit(X1, ys, lam, D), lam, D, G); Sd[b] <- s2$Sc; Se[b] <- s2$Sce }
  c(dec = (1 + sum(Sd >= st$Sc))/(NB + 1), edge = (1 + sum(Se >= st$Sce, na.rm = TRUE))/(NB + 1))
}

set.seed(1); cv <- cv.glmnet(X, y, family = "binomial", alpha = 0)
lmin <- cv$lambda.min*n
# 7101 is the seed T33 itself used for this cell (7000 + 10*G + 1)
seeds <- c(7101L, 9001L, 9002L, 9003L, 9004L, 9005L)
out <- t(vapply(seeds, function(s) run_unif(lmin, s), c(dec = 0, edge = 0)))
res <- data.frame(seed = seeds, p_dec = out[, "dec"], p_edge = out[, "edge"])
write.csv(res, "T36_draw_construction.csv", row.names = FALSE)
cat(sprintf("GlaucomaM, lambda.min = %.1f, G = %d, NB = %d, uniform-matrix draws\n", lmin, G, NB))
print(res, row.names = FALSE)
cat(sprintf("\ndecile: mean %.3f sd %.3f range [%.3f, %.3f]\n", mean(res$p_dec), sd(res$p_dec),
            min(res$p_dec), max(res$p_dec)))
cat(sprintf("EDGE  : mean %.3f sd %.3f range [%.3f, %.3f]\n", mean(res$p_edge), sd(res$p_edge),
            min(res$p_edge), max(res$p_edge)))
cat("T35's rbinom eight-seed ranges on this cell: decile [0.036, 0.066], EDGE [0.096, 0.136]\n")
