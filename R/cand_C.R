# =====================================================================
# CAND-C "second-order-mu": add the next (quadratic) term of E[r].
#
# Expansion:  pi(beta.hat) - pi(beta0)
#   = diag(w0) X1 d + (1/2) diag(w0 (1-2 pi0)) (X1 d)^2 + ...,
#   d = beta.hat - beta0,  E[d] = -M^{-1} K beta0,  Var(d) = M^{-1} F M^{-1}.
# Since r = V^{-1/2} C (y - pi.hat),
#   E[r] = mu1 + mu_q,
#   mu1  =  U M^{-1} K beta0                                (first order, existing)
#   mu_q = -V^{-1/2} C (1/2) diag(w0 (1-2 pi0)) q,
#   q_i  = (X1 E[d])_i^2 + [X1 Var(d) X1']_ii.
# Plug-in: quadratic ingredients evaluated at the DEBIASED fit beta.tilde
# (pi_t, w_t, F_t, M_t); grouping and V stay from the penalized fit so the
# object under test r is unchanged. Test: S_c2 = ||r - mu2||^2 vs N(0, Om_MLE).
#
# Steps: (0) glmnet cross-check, (1) numeric verification of the algebra at a
# FIXED panel-B design lambda=416 (must beat the first-order error),
# (2) 6-cell size sweep B=300, (3) power spot-check.
# =====================================================================
suppressPackageStartupMessages(library(glmnet))
set.seed(2026)
t_all <- Sys.time()

G  <- 10
ND <- 4000     # MC draws per reference p-value
B  <- 300      # replications per cell (analytic candidate, per brief)

## ---- penalized IRLS (verbatim from corrected_test.R) --------------------
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

## grouped residuals + Prop-1 matrices + CAND-C second-order mean
stats_for2 <- function(X1, y, fit, lam, D, g = NULL) {
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
  bt  <- fit$beta + drop(solve(Fm, K %*% fit$beta))       # one-step debias
  mu1 <- drop(t(MiU) %*% (K %*% bt))                      # first-order mu.hat
  ## ---- CAND-C: quadratic term, ingredients at the debiased fit ----------
  eta_t <- drop(X1 %*% bt)
  pi_t  <- 1/(1 + exp(-eta_t))
  w_t   <- pmax(pi_t*(1-pi_t), 1e-8)
  F_t   <- crossprod(X1, w_t * X1)
  M_t   <- F_t + K
  Ed    <- -drop(solve(M_t, K %*% bt))                    # E[d] plug-in
  Mi_t  <- solve(M_t)
  Vd    <- Mi_t %*% F_t %*% Mi_t                          # Var(d) plug-in
  q     <- drop(X1 %*% Ed)^2 + rowSums((X1 %*% Vd) * X1)
  h     <- 0.5 * w_t * (1 - 2*pi_t) * q
  mu_q  <- -as.numeric(tapply(h, g, sum)) / sqrt(Vg)
  list(r = r, Om_MLE = Om_MLE, mu1 = mu1, mu2 = mu1 + mu_q, mu_q = mu_q)
}

## MC reference p-values for several statistics against the SAME Sigma
pmc_multi <- function(Svec, Sig, nd = ND) {
  R  <- chol(Sig + diag(1e-9, nrow(Sig)))
  Z  <- matrix(rnorm(nd * nrow(Sig)), nd) %*% R
  ss <- rowSums(Z^2)
  vapply(Svec, function(S) mean(ss >= S), 0.0)
}

one_rep <- function(X1, y, lam, D) {
  fit <- ridge_fit(X1, y, lam, D)
  st  <- stats_for2(X1, y, fit, lam, D)
  S    <- sum(st$r^2)
  Sc1  <- sum((st$r - st$mu1)^2)
  Sc2  <- sum((st$r - st$mu2)^2)
  p <- pmc_multi(c(S, Sc1, Sc2), st$Om_MLE)   # all three referred to Om_MLE
  c(naive = p[1], corr1 = p[2], corr2 = p[3])
}

## ---- 0. one-off glmnet cross-check -------------------------------------
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

## ---- 1. VERIFY THE ALGEBRA at a fixed panel-B design, lambda = 416 ------
## (groups fixed at pi0-deciles so r has a stable definition across reps)
set.seed(11)
nB <- 400; pB <- 100; rho <- 0.7
SigB <- rho^abs(outer(1:pB, 1:pB, "-")); Ch <- chol(SigB)
betaB <- rep(c(0.35,-0.30,0.25,-0.20,0.15), length.out = pB) * 0.8887
DB    <- diag(c(0, rep(1, pB)))
lamV  <- 416

Xf  <- matrix(rnorm(nB*pB), nB, pB) %*% Ch; X1f <- cbind(1, Xf)
b0  <- c(0, betaB)
pi0 <- 1/(1+exp(-drop(X1f %*% b0)))
w0  <- pi0*(1-pi0)
g0  <- grp_of(pi0)
Vg0 <- as.numeric(tapply(w0, g0, sum))
U0  <- rowsum(w0 * X1f, g0) / sqrt(Vg0)
F0  <- crossprod(X1f, w0 * X1f)
K0  <- lamV * DB
M0  <- F0 + K0
mu1_0 <- drop(U0 %*% solve(M0, K0 %*% b0))                # first-order oracle
Ed0   <- -drop(solve(M0, K0 %*% b0))
Mi0   <- solve(M0)
Vd0   <- Mi0 %*% F0 %*% Mi0
q0    <- drop(X1f %*% Ed0)^2 + rowSums((X1f %*% Vd0) * X1f)
h0    <- 0.5 * w0 * (1 - 2*pi0) * q0
muq_0 <- -as.numeric(tapply(h0, g0, sum)) / sqrt(Vg0)
mu2_0 <- mu1_0 + muq_0

BV <- 4000
rmat  <- matrix(NA_real_, BV, G)
c1mat <- matrix(NA_real_, BV, G)   # r - mu.hat1 (plug-in first order)
c2mat <- matrix(NA_real_, BV, G)   # r - mu.hat2 (plug-in second order)
t0 <- Sys.time()
for (b in 1:BV) {
  yv <- rbinom(nB, 1, pi0)
  ft <- ridge_fit(X1f, yv, lamV, DB)
  st <- stats_for2(X1f, yv, ft, lamV, DB, g = g0)
  rmat[b, ]  <- st$r
  c1mat[b, ] <- st$r - st$mu1
  c2mat[b, ] <- st$r - st$mu2
}
emp <- colMeans(rmat)
relerr <- function(mu) sqrt(sum((emp - mu)^2)) / sqrt(sum(emp^2))
cat("=== CAND-C ALGEBRA VERIFICATION (fixed panel-B design, lam=416, B=4000) ===\n")
cat("||emp mean r|| =", round(sqrt(sum(emp^2)), 4),
    "  ||mu1_oracle|| =", round(sqrt(sum(mu1_0^2)), 4),
    "  ||mu_q_oracle|| =", round(sqrt(sum(muq_0^2)), 4), "\n")
cat("oracle  first-order : max|emp - mu1| =", round(max(abs(emp - mu1_0)), 4),
    "  relerr =", round(relerr(mu1_0), 4), "\n")
cat("oracle  second-order: max|emp - mu2| =", round(max(abs(emp - mu2_0)), 4),
    "  relerr =", round(relerr(mu2_0), 4), "\n")
cat("plug-in first-order : max|mean(r - mu.hat1)| =",
    round(max(abs(colMeans(c1mat))), 4), "\n")
cat("plug-in second-order: max|mean(r - mu.hat2)| =",
    round(max(abs(colMeans(c2mat))), 4), "\n")
cat("MC SE of a component mean ~", round(1/sqrt(BV), 4),
    " | elapsed:", round(as.numeric(difftime(Sys.time(), t0, units="mins")),1),
    "min\n\n")

## ---- 2. SIZE SWEEP: all 6 target cells, H0 TRUE, B = 300 ----------------
run_panel <- function(gen, lam_grid, D, label) {
  res <- array(NA_real_, c(B, length(lam_grid), 3),
               dimnames = list(NULL, lam_grid, c("naive","corr1","corr2")))
  t0 <- Sys.time()
  for (b in 1:B) {
    dat <- gen()
    for (j in seq_along(lam_grid))
      res[b, j, ] <- one_rep(dat$X1, dat$y, lam_grid[j], D)
  }
  rej <- apply(res < 0.05, c(2,3), mean)
  cat("=== ", label, " (H0 TRUE, alpha=0.05, B =", B, ") ===\n")
  print(round(rej, 3))
  cat("elapsed:", round(as.numeric(difftime(Sys.time(), t0, units="mins")),1),
      "min\n\n")
  rej
}

set.seed(2026)
genA <- function() {
  X <- matrix(rnorm(n*p), n, p)
  list(X1 = cbind(1, X), y = rbinom(n, 1, 1/(1+exp(-drop(X %*% beta_true)))))
}
rejA <- run_panel(genA, c(100, 137, 200), D, "PANEL A (n=500, p=5, iid)")

set.seed(31)
genB <- function() {
  Xb <- matrix(rnorm(nB*pB), nB, pB) %*% Ch
  list(X1 = cbind(1, Xb), y = rbinom(nB, 1, 1/(1+exp(-drop(Xb %*% betaB)))))
}
rejB <- run_panel(genB, c(50, 416, 1000), DB, "PANEL B (n=400, p=100, AR(1) 0.7)")

## ---- 3. POWER spot-check: omitted quadratic, lam_eq = 25 ----------------
set.seed(21)
pw2 <- pw_mle <- numeric(B)
t0 <- Sys.time()
for (b in 1:B) {
  X  <- matrix(rnorm(n*p), n, p); X1 <- cbind(1, X)
  eta <- drop(X %*% beta_true) + 0.5*(X[,1]^2 - 1)
  y  <- rbinom(n, 1, 1/(1+exp(-eta)))
  ft <- ridge_fit(X1, y, 25, D)
  st <- stats_for2(X1, y, ft, 25, D)
  pw2[b] <- pmc_multi(sum((st$r - st$mu2)^2), st$Om_MLE) < 0.05
  fm <- ridge_fit(X1, y, 0, D)
  sm <- stats_for2(X1, y, fm, 0, D)
  pw_mle[b] <- pmc_multi(sum(sm$r^2), sm$Om_MLE) < 0.05
}
cat("=== POWER (omitted quadratic, gamma=0.5, B =", B, ") ===\n")
cat("CAND-C corrected test on ridge fit (lam_eq=25):", round(mean(pw2), 3), "\n")
cat("standard test on MLE fit (comparator)        :", round(mean(pw_mle), 3), "\n")
cat("elapsed:", round(as.numeric(difftime(Sys.time(), t0, units="mins")),1),
    "min\n\n")

## ---- 4. verdict + save --------------------------------------------------
mcse <- function(x) sqrt(x*(1-x)/B)
cells <- rbind(
  data.frame(panel="A", lambda=c(100,137,200), rej=rejA[, "corr2"]),
  data.frame(panel="B", lambda=c(50,416,1000), rej=rejB[, "corr2"]))
cells$mcse <- mcse(cells$rej)
cells$pass <- cells$rej >= 0.03 & cells$rej <= 0.08
cat("=== CAND-C PER-CELL VERDICT (corr2, target [0.03, 0.08]) ===\n")
print(cells, row.names = FALSE)
cat("OVERALL:", if (all(cells$pass)) "PASS" else "FAIL", "\n")

saveRDS(list(rejA = rejA, rejB = rejB, cells = cells, B = B,
             verif = list(emp = emp, mu1_0 = mu1_0, mu2_0 = mu2_0,
                          plug1 = colMeans(c1mat), plug2 = colMeans(c2mat)),
             power = c(cand_C = mean(pw2), mle = mean(pw_mle))),
        "cand_C_results.rds")
cat("\nSaved cand_C_results.rds. TOTAL",
    round(as.numeric(difftime(Sys.time(), t_all, units="mins")),1), "min\n")
warns <- warnings()
if (length(warns)) { cat("WARNINGS:\n"); print(warns) } else cat("No warnings.\n")
