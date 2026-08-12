# =====================================================================
# PROOF OF CONCEPT: the corrected GOF test of Proposition 1, built and
# verified end to end.
#
#   1. Own penalized-IRLS ridge fitter (penalty (lambda/2) b' D b,
#      D = diag(0,1,...,1)) -- no glmnet scale ambiguity; cross-checked
#      against glmnet once below.
#   2. Grouped standardized residuals r (HL-type deciles, G = 10).
#   3. Three references for S:
#        naive     : S = r'r        vs  N(0, Omega_MLE)   <- current practice
#        Omega-only: S = r'r        vs  N(0, Omega_K)     <- covariance fix alone
#        corrected : S_c=|r-mu.hat|^2 vs N(0, Omega_MLE)  <- full Prop.1 fix
#      where mu.hat = U M^{-1} K beta.tilde and
#      beta.tilde = beta.hat + F^{-1} K beta.hat  (one-step debias).
#      LEMMA (verified below): Var(r - mu.hat) = Omega_MLE to first order,
#      because M^{-1}(I + K F^{-1}) = F^{-1}.
#   4. Direct verification of Proposition 1 at a FIXED design:
#      empirical mean(r) vs mu_K and empirical Cov(r) vs Omega_K.
#   5. Size sweep (H0 true) + power check (omitted quadratic).
# =====================================================================
suppressPackageStartupMessages(library(glmnet))
set.seed(2026)

G  <- 10
ND <- 4000     # Monte Carlo draws per reference p-value
B  <- 500      # replications per cell

## ---- penalized IRLS ----------------------------------------------------
ridge_fit <- function(X1, y, lam, D, tol = 1e-10, maxit = 100) {
  beta <- rep(0, ncol(X1))
  for (it in 1:maxit) {
    eta <- drop(X1 %*% beta)
    pi  <- 1/(1 + exp(-eta))
    w   <- pmax(pi*(1-pi), 1e-8)
    z   <- eta + (y - pi)/w
    A   <- crossprod(X1, w * X1) + lam * D
    bn  <- drop(solve(A, crossprod(X1, w * z)))
    if (max(abs(bn - beta)) < tol) { beta <- bn; break }
    beta <- bn
  }
  eta <- drop(X1 %*% beta)
  list(beta = beta, pi = 1/(1 + exp(-eta)))
}

grp_of <- function(pi) ceiling(rank(pi, ties.method = "first") * G / length(pi))

## grouped residuals + all Proposition-1 matrices at the fitted values
stats_for <- function(X1, y, fit, lam, D, g = NULL) {
  pi <- fit$pi
  w  <- pmax(pi*(1-pi), 1e-8)
  if (is.null(g)) g <- grp_of(pi)
  Vg  <- as.numeric(tapply(w, g, sum))
  r   <- as.numeric(tapply(y - pi, g, sum)) / sqrt(Vg)
  U   <- rowsum(w * X1, g) / sqrt(Vg)
  Fm  <- crossprod(X1, w * X1)
  K   <- lam * D
  M   <- Fm + K
  MiU <- solve(M, t(U))                                   # M^{-1} U'
  Om_MLE <- diag(G) - U %*% solve(Fm, t(U))
  Om_nv  <- diag(G) - U %*% MiU
  Om_K   <- Om_nv - t(MiU) %*% K %*% MiU
  bt  <- fit$beta + drop(solve(Fm, K %*% fit$beta))       # one-step debias
  mu  <- drop(t(MiU) %*% (K %*% bt))
  list(r = r, Om_MLE = Om_MLE, Om_K = Om_K, mu = mu)
}

pmc <- function(S, Sig, nd = ND) {                        # MC ref p-value
  R <- chol(Sig + diag(1e-9, nrow(Sig)))
  Z <- matrix(rnorm(nd * nrow(Sig)), nd) %*% R
  mean(rowSums(Z^2) >= S)
}

one_rep <- function(X1, y, lam, D) {
  fit <- ridge_fit(X1, y, lam, D)
  st  <- stats_for(X1, y, fit, lam, D)
  S   <- sum(st$r^2)
  Sc  <- sum((st$r - st$mu)^2)
  c(naive = pmc(S,  st$Om_MLE),
    omega = pmc(S,  st$Om_K),
    corr  = pmc(Sc, st$Om_MLE))
}

## ---- 0. one-off glmnet cross-check ------------------------------------
n <- 500; p <- 5
beta_true <- c(0.5, -0.4, 0.3, -0.3, 0.2)
D  <- diag(c(0, rep(1, p)))
X  <- matrix(rnorm(n*p), n, p); X1 <- cbind(1, X)
y  <- rbinom(n, 1, 1/(1+exp(-drop(X %*% beta_true))))
f_own <- ridge_fit(X1, y, 25, D)
f_gn  <- glmnet(X, y, family="binomial", alpha=0, lambda=25/n,
                standardize=FALSE, thresh=1e-14)
cat("glmnet cross-check (lam_eq=25): max|coef diff| =",
    format(max(abs(f_own$beta - as.numeric(coef(f_gn)))), digits=3), "\n\n")

## ---- 1. Panel A size sweep under H0 ------------------------------------
lam_grid <- c(0, 10, 25, 50, 100, 137)
resA <- array(NA_real_, c(B, length(lam_grid), 3),
              dimnames = list(NULL, lam_grid, c("naive","omega","corr")))
t0 <- Sys.time()
for (b in 1:B) {
  X  <- matrix(rnorm(n*p), n, p); X1 <- cbind(1, X)
  y  <- rbinom(n, 1, 1/(1+exp(-drop(X %*% beta_true))))
  for (j in seq_along(lam_grid)) resA[b, j, ] <- one_rep(X1, y, lam_grid[j], D)
}
rejA <- apply(resA < 0.05, c(2,3), mean)
cat("=== PANEL A (n=500, p=5, H0 TRUE): rejection at alpha=0.05, B =", B, "===\n")
print(round(rejA, 3))
cat("MC SE at 0.05 ~", round(sqrt(.05*.95/B), 4), "\n",
    "elapsed:", round(as.numeric(difftime(Sys.time(), t0, units="mins")),1), "min\n\n")

## ---- 2. Proposition-1 verification at a FIXED design -------------------
## (groups fixed at pi0-deciles so r has a stable definition across reps)
set.seed(7)
Xf  <- matrix(rnorm(n*p), n, p); X1f <- cbind(1, Xf)
pi0 <- 1/(1+exp(-drop(Xf %*% beta_true)))
g0  <- grp_of(pi0)
w0  <- pi0*(1-pi0)
lamV <- 100
Vg0 <- as.numeric(tapply(w0, g0, sum))
U0  <- rowsum(w0 * X1f, g0) / sqrt(Vg0)
F0  <- crossprod(X1f, w0 * X1f)
K0  <- lamV * D
M0  <- F0 + K0
MiU0 <- solve(M0, t(U0))
mu0    <- drop(t(MiU0) %*% (K0 %*% c(0, beta_true)))
OmK0   <- diag(G) - U0 %*% MiU0 - t(MiU0) %*% K0 %*% MiU0
OmNV0  <- diag(G) - U0 %*% MiU0
OmMLE0 <- diag(G) - U0 %*% solve(F0, t(U0))

BV <- 4000
rmat  <- matrix(NA_real_, BV, G)
rcmat <- matrix(NA_real_, BV, G)          # bias-corrected r - mu.hat
for (b in 1:BV) {
  y  <- rbinom(n, 1, pi0)
  ft <- ridge_fit(X1f, y, lamV, D)
  st <- stats_for(X1f, y, ft, lamV, D, g = g0)
  rmat[b, ]  <- st$r
  rcmat[b, ] <- st$r - st$mu
}
frob <- function(A, Bm) sqrt(sum((A-Bm)^2))/sqrt(sum(Bm^2))
cat("=== PROPOSITION 1 VERIFICATION (fixed design, lam_eq=100, B=4000) ===\n")
cat("mean shift  : max|emp mean r  - mu_K|      =",
    round(max(abs(colMeans(rmat) - mu0)), 4),
    "   (||mu_K|| =", round(sqrt(sum(mu0^2)), 3), ")\n")
cat("covariance  : relFrob(emp Cov r, Omega_K)  =", round(frob(cov(rmat),  OmK0),  3), "\n")
cat("              relFrob(emp Cov r, naive)    =", round(frob(cov(rmat),  OmNV0), 3), "\n")
cat("              relFrob(emp Cov r, Omega_MLE)=", round(frob(cov(rmat),  OmMLE0),3), "\n")
cat("LEMMA check : max|emp mean (r-mu.hat)|     =", round(max(abs(colMeans(rcmat))), 4), "\n")
cat("              relFrob(emp Cov(r-mu.hat), Omega_MLE) =",
    round(frob(cov(rcmat), OmMLE0), 3), "\n\n")

## ---- 3. Panel B (n=400, p=100, AR(1) 0.7) ------------------------------
set.seed(11)
nB <- 400; pB <- 100; rho <- 0.7
Sig <- rho^abs(outer(1:pB, 1:pB, "-")); Ch <- chol(Sig)
betaB <- rep(c(0.35,-0.30,0.25,-0.20,0.15), length.out = pB) * 0.8887
DB  <- diag(c(0, rep(1, pB)))
lamB_grid <- c(20, 50, 416)               # lam_g = 0.05, ~cv lambda.min, ~cv lambda.1se
resB <- array(NA_real_, c(B, length(lamB_grid), 3),
              dimnames = list(NULL, lamB_grid, c("naive","omega","corr")))
for (b in 1:B) {
  Xb  <- matrix(rnorm(nB*pB), nB, pB) %*% Ch; X1b <- cbind(1, Xb)
  yb  <- rbinom(nB, 1, 1/(1+exp(-drop(Xb %*% betaB))))
  for (j in seq_along(lamB_grid)) resB[b, j, ] <- one_rep(X1b, yb, lamB_grid[j], DB)
}
rejB <- apply(resB < 0.05, c(2,3), mean)
cat("=== PANEL B (n=400, p=100, AR(1) 0.7, H0 TRUE): rejection at 0.05 ===\n")
print(round(rejB, 3))

## ---- 4. Power: omitted quadratic, lam_eq = 25 ---------------------------
set.seed(21)
pw_corr <- pw_mle <- numeric(B)
for (b in 1:B) {
  X  <- matrix(rnorm(n*p), n, p); X1 <- cbind(1, X)
  eta <- drop(X %*% beta_true) + 0.5*(X[,1]^2 - 1)   # omitted quadratic
  y  <- rbinom(n, 1, 1/(1+exp(-eta)))
  # corrected test on the ridge fit
  ft <- ridge_fit(X1, y, 25, D)
  st <- stats_for(X1, y, ft, 25, D)
  pw_corr[b] <- pmc(sum((st$r - st$mu)^2), st$Om_MLE) < 0.05
  # standard test on the MLE fit (comparator)
  fm <- ridge_fit(X1, y, 0, D)
  sm <- stats_for(X1, y, fm, 0, D)
  pw_mle[b] <- pmc(sum(sm$r^2), sm$Om_MLE) < 0.05
}
cat("\n=== POWER (omitted quadratic, gamma=0.5): ===\n")
cat("corrected test on ridge fit (lam_eq=25):", round(mean(pw_corr), 3), "\n")
cat("standard  test on MLE   fit           :", round(mean(pw_mle),  3), "\n")

saveRDS(list(rejA=rejA, rejB=rejB,
             verif=list(mu_err=max(abs(colMeans(rmat)-mu0)), mu_norm=sqrt(sum(mu0^2)),
                        fK=frob(cov(rmat),OmK0), fNV=frob(cov(rmat),OmNV0),
                        fMLE=frob(cov(rmat),OmMLE0),
                        lem_mu=max(abs(colMeans(rcmat))), lem_cov=frob(cov(rcmat),OmMLE0)),
             power=c(corr=mean(pw_corr), mle=mean(pw_mle))),
        "corrected_test_results.rds")
cat("\nSaved corrected_test_results.rds. Total",
    round(as.numeric(difftime(Sys.time(), t0, units="mins")),1), "min\n")
