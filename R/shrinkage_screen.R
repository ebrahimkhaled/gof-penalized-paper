# W10(a): THE SHRINKAGE-INFLATION SCREEN.
#
# Every ingredient is already computed inside the correction, so this is a corollary of
# Proposition 1, not a new theorem. Define
#     N = ||mu.hat||^2 ,  mu.hat = U M^{-1} K beta.tilde     (the shrinkage non-centrality)
#     d = tr(I_G - U F^{-1} U')                              (effective residual dimension)
# and predict the UNCORRECTED test's size as
#     alpha.naive = P{ Q_K > q } ,  q = 95th pct of the central weighted chi-square with
#                                       weights eig(Omega_MLE),
#                                   Q_K = weighted NON-central form, weights eig(Omega_K).
# Davies/Imhof exact (CompQuadForm), with the Satterthwaite chi^2_d(N) version reported too.
#
# VALIDATION is the point: predict all 12 registered cells of manifest section 2 and compare
# against the MEASURED naive sizes. Reported honestly, including where it fails.
suppressPackageStartupMessages({library(CompQuadForm)})
G <- 10; NDRAW <- 200

ridge_fit <- function(X1, y, lam, D, tol = 1e-9, maxit = 100) {
  beta <- rep(0, ncol(X1))
  for (it in 1:maxit) {
    eta <- drop(X1 %*% beta); pi <- 1/(1+exp(-eta))
    w <- pmax(pi*(1-pi), 1e-8); z <- eta + (y-pi)/w
    bn <- drop(solve(crossprod(X1, w*X1) + lam*D, crossprod(X1, w*z)))
    if (max(abs(bn-beta)) < tol) { beta <- bn; break }; beta <- bn
  }
  list(beta = beta, pi = 1/(1+exp(-drop(X1 %*% beta))))
}

screen_one <- function(X1, y, lam, D) {
  ft <- ridge_fit(X1, y, lam, D)
  pi <- pmin(pmax(ft$pi, 1e-8), 1-1e-8); w <- pmax(pi*(1-pi), 1e-10)
  g <- pmin(ceiling(rank(pi, ties.method = "first")/(length(y)/G)), G)
  idx <- split(seq_along(y), g)
  Vg <- vapply(idx, function(I) sum(w[I]), 0.0)
  U  <- t(vapply(idx, function(I) colSums(w[I]*X1[I,,drop=FALSE]), numeric(ncol(X1))))/sqrt(Vg)
  F  <- crossprod(X1, w*X1); K <- lam*D; M <- F + K
  bt <- ft$beta + drop(solve(F, K %*% ft$beta))
  mu <- drop(U %*% solve(M, K %*% bt))
  Om_ML <- diag(G) - U %*% solve(F,  t(U))
  Om_K  <- diag(G) - U %*% solve(M, (F + 2*K) %*% solve(M, t(U)))
  N <- sum(mu^2); d <- sum(diag(Om_ML))
  lML <- pmax(eigen(Om_ML, symmetric = TRUE, only.values = TRUE)$values, 0)
  lK  <- pmax(eigen(Om_K,  symmetric = TRUE, only.values = TRUE)$values, 0)
  # critical value: 95th pct of the CENTRAL weighted chi-square with Omega_MLE weights
  q <- tryCatch(uniroot(function(x) davies(x, lML, rep(1, G))$Qq - 0.05,
                        c(1e-6, 500))$root, error = function(e) NA_real_)
  # size: the NON-central weighted form under Omega_K, non-centrality spread on eig(Om_K)
  dl <- rep(0, G); dl[which.max(lK)] <- N/max(lK)      # place N on the leading direction
  a_exact <- if (is.na(q)) NA_real_ else
    tryCatch(davies(q, lK, rep(1, G), delta = dl)$Qq, error = function(e) NA_real_)
  a_satt <- pchisq(qchisq(0.95, d), d, ncp = N, lower.tail = FALSE)
  c(N = N, d = d, ratio = N/d, q = q, alpha_exact = a_exact, alpha_satt = a_satt)
}

gen <- function(panel) {
  if (panel == "A") {
    n <- 500; p <- 5; b0 <- c(0.5,-0.4,0.3,-0.3,0.2); X <- matrix(rnorm(n*p), n, p)
  } else {
    n <- 400; p <- 100; S <- 0.7^abs(outer(1:p,1:p,"-")); Ch <- chol(S)
    b0 <- rep(c(.35,-.30,.25,-.20,.15), length.out = p)*0.8887
    X <- matrix(rnorm(n*p), n, p) %*% Ch
  }
  list(X1 = cbind(1,X), y = rbinom(n, 1, 1/(1+exp(-drop(X %*% b0)))),
       D = diag(c(0, rep(1, p))), n = n)
}

cells <- list(
  list(p="A", lg=c(0,0.02,0.05,0.10,0.20,0.274), meas=c(0.068,0.068,0.118,0.362,0.916,0.988)),
  list(p="B", lg=c(0,0.01,0.05,0.124,0.513,1.5), meas=c(0.178,0.432,0.990,0.998,1.000,1.000)))   # 0.124: T38, was typed 0.994

set.seed(20260812)
out <- NULL
for (ce in cells) for (k in seq_along(ce$lg)) {
  lg <- ce$lg[k]
  acc <- replicate(NDRAW, { d <- gen(ce$p); screen_one(d$X1, d$y, lg*d$n, d$D) })
  m <- rowMeans(acc, na.rm = TRUE)
  out <- rbind(out, data.frame(panel = ce$p, lambda_g = lg,
    N = m["N"], d = m["d"], ratio = m["ratio"],
    pred_exact = m["alpha_exact"], pred_satt = m["alpha_satt"], measured = ce$meas[k]))
}
rownames(out) <- NULL
print(out, digits = 3, row.names = FALSE)
ok <- !is.na(out$pred_exact)
cat(sprintf("\nmean |predicted - measured|, exact : %.3f  (A only: %.3f ; B only: %.3f)\n",
    mean(abs(out$pred_exact-out$measured)[ok]),
    mean(abs(out$pred_exact-out$measured)[ok & out$panel=="A"]),
    mean(abs(out$pred_exact-out$measured)[ok & out$panel=="B"])))
cat(sprintf("rule of thumb N > d  =>  size > 0.5 : correct in %d of %d cells\n",
    sum((out$ratio > 1) == (out$measured > 0.5)), nrow(out)))
write.csv(out, "shrinkage_screen.csv", row.names = FALSE)
cat("wrote shrinkage_screen.csv\n")
