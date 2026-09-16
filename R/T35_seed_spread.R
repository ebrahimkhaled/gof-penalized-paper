# T3.5 -- BJ revision (co-author comment 3, the application's half): how much does one corrected
# p-value on GlaucomaM move between bootstrap streams at NB = 499? Everything is held fixed --
# the data, lambda, G, the statistic -- and only the random stream changes, eight seeds per cell.
# This is the quantity a reader needs in order to know how far to trust a single entry of the
# glaucoma table, and it is also why an independent run of the same grid (T33 at tau = 0) does not
# reproduce that table cell for cell.
# Output: T35_seed_spread.csv (one row per seed) and the summary below.
suppressPackageStartupMessages({library(glmnet); library(TH.data)})
data("GlaucomaM", package = "TH.data")
Xr <- as.matrix(GlaucomaM[, setdiff(names(GlaucomaM), "Class")])
y  <- as.numeric(GlaucomaM$Class == "glaucoma")
X  <- scale(Xr); X <- X[, apply(X, 2, function(c) all(is.finite(c)) & sd(c) > 0), drop = FALSE]
n <- nrow(X); p <- ncol(X); X1 <- cbind(1, X); D <- diag(c(0, rep(1, p)))
NB <- 499; SEEDS <- 1:8

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
run <- function(lam, G, seed) {
  set.seed(seed)
  ft <- ridge_fit(X1, y, lam, D); st <- pieces(X1, y, ft, lam, D, G)
  gen <- pmin(pmax(as.numeric(1/(1 + exp(-X1 %*% st$bt))), 1e-6), 1 - 1e-6)
  Sd <- Se <- numeric(NB)
  for (b in 1:NB) { ys <- rbinom(n, 1, gen)
    s2 <- pieces(X1, ys, ridge_fit(X1, ys, lam, D), lam, D, G); Sd[b] <- s2$Sc; Se[b] <- s2$Sce }
  c(dec = (1 + sum(Sd >= st$Sc))/(NB + 1), edge = (1 + sum(Se >= st$Sce, na.rm = TRUE))/(NB + 1))
}

set.seed(1); cv <- cv.glmnet(X, y, family = "binomial", alpha = 0)
lams <- c("lambda.min" = cv$lambda.min*n, "lambda.1se" = cv$lambda.1se*n)
cat(sprintf("GlaucomaM: lambda.min %.1f, lambda.1se %.1f (theory scale), NB = %d\n\n",
            lams[1], lams[2], NB))

out <- NULL; t0 <- Sys.time()
for (nm in names(lams)) for (G in c(10)) {
  r <- t(vapply(SEEDS, function(s) run(lams[[nm]], G, s), c(dec = 0, edge = 0)))
  out <- rbind(out, data.frame(lambda = nm, lam = round(lams[[nm]], 1), G = G,
                               seed = SEEDS, p_dec = r[, "dec"], p_edge = r[, "edge"]))
  write.csv(out, "T35_seed_spread.csv", row.names = FALSE)
  cat(sprintf("%s (lambda %.1f), G = %d, %d seeds\n", nm, lams[[nm]], G, length(SEEDS)))
  for (b in c("dec", "edge")) {
    v <- r[, b]
    cat(sprintf("  %-4s mean %.3f  sd %.3f  range [%.3f, %.3f]  | binomial sd %.3f\n",
                b, mean(v), sd(v), min(v), max(v), sqrt(mean(v)*(1 - mean(v))/NB)))
  }
  cat("\n")
}
cat("runtime:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
