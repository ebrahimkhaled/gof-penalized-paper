# =====================================================================
# Pilot study: the uncorrected test along the ridge penalty path.
#
# Goal: replace the two-arbitrary-lambda table with a lambda SWEEP that
#   (a) shows size distortion as a smooth function of the penalty,
#   (b) marks where cross-validation actually lands (lambda.min / 1se),
#   (c) separates the two mechanisms:
#         - reference-distribution invalidation (statistic distribution shifts)
#         - genuine shrinkage-induced miscalibration (fitted probs biased),
#   (d) reports Monte Carlo standard errors.
#
# H0 IS TRUE THROUGHOUT: the linear logistic model is correctly specified.
# =====================================================================
suppressPackageStartupMessages({ library(ebrahim.gof); library(glmnet) })

set.seed(2026)

R      <- 500          # replications
n      <- 500
p      <- 5
beta   <- c(0.5, -0.4, 0.3, -0.3, 0.2)
alpha  <- 0.05

# glmnet lambda grid. glmnet minimises -(1/n) loglik + lambda * P(beta),
# so the equivalent penalty in the proposal's  l(beta) - lambda_eq * P(beta)
# parameterisation is lambda_eq = lambda_glmnet * n.
lam_grid <- c(0, 0.005, 0.01, 0.02, 0.0268, 0.05, 0.10, 0.15, 0.20, 0.274, 0.40)
L <- length(lam_grid)

# ---- Hosmer-Lemeshow (same implementation as the original pilot, with a
# ---- degeneracy guard so we can COUNT how often ridge compression breaks it)
hl_p <- function(y, pi_hat, g = 10) {
  br <- quantile(pi_hat, probs = seq(0, 1, length.out = g + 1))
  degenerate <- anyDuplicated(br) > 0
  cuts <- tryCatch(cut(pi_hat, breaks = br, include.lowest = TRUE),
                   error = function(e) NULL)
  if (is.null(cuts)) return(c(p = NA_real_, stat = NA_real_, degen = 1))
  obs <- tapply(y, cuts, sum); exp <- tapply(pi_hat, cuts, sum)
  ng  <- tapply(y, cuts, length)
  den <- exp * (1 - exp / ng); den[den <= 0] <- 1e-6
  st  <- sum((obs - exp)^2 / den, na.rm = TRUE)
  c(p = pchisq(st, df = g - 2, lower.tail = FALSE), stat = st,
    degen = as.numeric(degenerate))
}

# storage
p_hl   <- matrix(NA_real_, R, L); s_hl   <- matrix(NA_real_, R, L)
p_edge <- matrix(NA_real_, R, L); s_edge <- matrix(NA_real_, R, L)
degen  <- matrix(NA_real_, R, L)
sd_pi  <- matrix(NA_real_, R, L)   # spread of fitted probabilities
mae_pi <- matrix(NA_real_, R, L)   # |pi_hat - pi_TRUE| : direct bias measure
lam_min <- numeric(R); lam_1se <- numeric(R)

cat("Running", R, "replications over", L, "penalty levels...\n")
t0 <- Sys.time()

for (b in 1:R) {
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("X", 1:p)
  pi_true <- as.numeric(1 / (1 + exp(-(X %*% beta))))
  y <- rbinom(n, 1, pi_true)
  d <- data.frame(y = y, X)

  fit_mle <- glm(y ~ ., data = d, family = binomial)

  cv <- cv.glmnet(X, y, family = "binomial", alpha = 0)
  lam_min[b] <- cv$lambda.min; lam_1se[b] <- cv$lambda.1se

  for (j in 1:L) {
    lam <- lam_grid[j]
    if (lam == 0) {
      pi_hat <- as.numeric(fit_mle$fitted.values)
      eta_hat <- as.numeric(fit_mle$linear.predictors)
    } else {
      g <- glmnet(X, y, family = "binomial", alpha = 0, lambda = lam)
      pi_hat  <- as.numeric(predict(g, newx = X, type = "response"))
      eta_hat <- as.numeric(predict(g, newx = X, type = "link"))
    }

    sd_pi[b, j]  <- sd(pi_hat)
    mae_pi[b, j] <- mean(abs(pi_hat - pi_true))

    h <- hl_p(y, pi_hat)
    p_hl[b, j] <- h["p"]; s_hl[b, j] <- h["stat"]; degen[b, j] <- h["degen"]

    # "naive practice": feed the penalised fitted probabilities to a test
    # whose null distribution was derived under MLE.
    fp <- fit_mle
    fp$fitted.values <- pi_hat
    fp$linear.predictors <- eta_hat
    e <- tryCatch(edge.gof(fp), error = function(err) NULL)
    if (!is.null(e)) { p_edge[b, j] <- e$p_value; s_edge[b, j] <- e$Test_Statistic }
  }
  if (b %% 50 == 0) cat("  rep", b, "/", R,
                        sprintf(" (%.1f min elapsed)\n",
                                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

rej <- function(m) apply(m, 2, function(v) mean(na.omit(v) < alpha))
mcse <- function(r, k) sqrt(r * (1 - r) / k)
nvalid <- function(m) apply(m, 2, function(v) sum(!is.na(v)))

r_hl <- rej(p_hl); r_ed <- rej(p_edge)

out <- data.frame(
  lambda_glmnet = lam_grid,
  lambda_eq     = lam_grid * n,
  rej_HL        = round(r_hl, 3),
  mcse_HL       = round(mcse(r_hl, nvalid(p_hl)), 4),
  rej_EDGE      = round(r_ed, 3),
  mcse_EDGE     = round(mcse(r_ed, nvalid(p_edge)), 4),
  mean_stat_EDGE= round(colMeans(s_edge, na.rm = TRUE), 2),
  mean_stat_HL  = round(colMeans(s_hl,   na.rm = TRUE), 2),
  sd_fitted_pi  = round(colMeans(sd_pi),  4),
  MAE_vs_true   = round(colMeans(mae_pi), 4),
  HL_degenerate = round(colMeans(degen, na.rm = TRUE), 3)
)

cat("\n===================== LAMBDA SWEEP (H0 TRUE) =====================\n")
print(out, row.names = FALSE)
cat("\nCross-validated penalties across", R, "replications:\n")
cat(sprintf("  lambda.min : median %.4f   IQR [%.4f, %.4f]\n",
            median(lam_min), quantile(lam_min, .25), quantile(lam_min, .75)))
cat(sprintf("  lambda.1se : median %.4f   IQR [%.4f, %.4f]\n",
            median(lam_1se), quantile(lam_1se, .25), quantile(lam_1se, .75)))
cat(sprintf("\nEDGE reference null is chi-square with about 2 df (mean ~2).\n"))
cat("==================================================================\n")

write.csv(out, "lambda_sweep_results.csv", row.names = FALSE)
saveRDS(list(out = out, lam_min = lam_min, lam_1se = lam_1se,
             s_edge = s_edge, p_edge = p_edge, p_hl = p_hl,
             lam_grid = lam_grid, R = R, n = n, p = p, beta = beta),
        "lambda_sweep_full.rds")
cat("\nSaved lambda_sweep_results.csv and lambda_sweep_full.rds\n")
cat("Total runtime:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
sessionInfo()$R.version$version.string
