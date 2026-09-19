# =====================================================================
# HIGH-DIMENSIONAL ARM: n = 400, p = 100 (kappa = p/n = 0.25),
# AR(1) correlated design with rho = 0.7 -- the regime where ridge is
# the RIGHT estimator.
# H0 IS TRUE: the linear logistic model is correctly specified.
# =====================================================================
suppressPackageStartupMessages({ library(ebrahim.gof); library(glmnet) })

set.seed(2026)

R <- 500; n <- 400; p <- 100; rho <- 0.7; alpha <- 0.05
Sig  <- rho^abs(outer(1:p, 1:p, "-"))
Ch   <- chol(Sig)
beta <- rep(c(0.35, -0.30, 0.25, -0.20, 0.15), length.out = p) * 0.8887  # sd(eta) ~ 1.5

lam_grid <- c(0, 0.01, 0.025, 0.05, 0.064, 0.10, 0.20, 0.40, 0.513, 0.80, 1.50)
L <- length(lam_grid)

hl_p <- function(y, pi_hat, g = 10) {
  br <- quantile(pi_hat, probs = seq(0, 1, length.out = g + 1))
  if (anyDuplicated(br) > 0) return(c(p = NA_real_, stat = NA_real_, degen = 1))
  cuts <- cut(pi_hat, breaks = br, include.lowest = TRUE)
  obs <- tapply(y, cuts, sum); exp <- tapply(pi_hat, cuts, sum)
  ng  <- tapply(y, cuts, length)
  den <- exp * (1 - exp / ng)
  if (any(den <= 0, na.rm = TRUE)) return(c(p = NA_real_, stat = NA_real_, degen = 1))
  st  <- sum((obs - exp)^2 / den, na.rm = TRUE)
  c(p = pchisq(st, df = g - 2, lower.tail = FALSE), stat = st, degen = 0)
}

# calibration slope: logistic regression of y on logit(pi_hat). 1 = perfect.
cal_slope <- function(y, pi_hat) {
  e <- pmin(pmax(pi_hat, 1e-8), 1 - 1e-8)
  lp <- log(e / (1 - e))
  f <- tryCatch(glm(y ~ lp, family = binomial), error = function(err) NULL)
  if (is.null(f)) return(c(slope = NA_real_, citl = NA_real_))
  f0 <- tryCatch(glm(y ~ offset(lp), family = binomial), error = function(err) NULL)
  c(slope = unname(coef(f)[2]),
    citl  = if (is.null(f0)) NA_real_ else unname(coef(f0)[1]))
}

p_hl <- s_hl <- p_ed <- s_ed <- degen <- sdpi <- mae <- slp <- ctl <-
  matrix(NA_real_, R, L)
lam_min <- lam_1se <- numeric(R)

cat("HIGH-DIM ARM:", R, "reps, n =", n, "p =", p, "kappa =", p/n,
    "AR(1) rho =", rho, "\n")
t0 <- Sys.time()

for (b in 1:R) {
  X <- matrix(rnorm(n * p), n, p) %*% Ch; colnames(X) <- paste0("X", 1:p)
  pi_true <- as.numeric(1 / (1 + exp(-(X %*% beta))))
  y <- rbinom(n, 1, pi_true)
  d <- data.frame(y = y, X)

  fit_mle <- tryCatch(glm(y ~ ., data = d, family = binomial),
                      error = function(e) NULL)
  if (is.null(fit_mle)) next

  cv <- cv.glmnet(X, y, family = "binomial", alpha = 0)
  lam_min[b] <- cv$lambda.min; lam_1se[b] <- cv$lambda.1se

  for (j in 1:L) {
    lam <- lam_grid[j]
    if (lam == 0) {
      pi_hat  <- as.numeric(fit_mle$fitted.values)
      eta_hat <- as.numeric(fit_mle$linear.predictors)
    } else {
      g <- glmnet(X, y, family = "binomial", alpha = 0, lambda = lam)
      pi_hat  <- as.numeric(predict(g, newx = X, type = "response"))
      eta_hat <- as.numeric(predict(g, newx = X, type = "link"))
    }
    sdpi[b, j] <- sd(pi_hat); mae[b, j] <- mean(abs(pi_hat - pi_true))
    cs <- cal_slope(y, pi_hat); slp[b, j] <- cs["slope"]; ctl[b, j] <- cs["citl"]

    h <- hl_p(y, pi_hat)
    p_hl[b, j] <- h["p"]; s_hl[b, j] <- h["stat"]; degen[b, j] <- h["degen"]

    fp <- fit_mle
    fp$fitted.values <- pi_hat; fp$linear.predictors <- eta_hat
    e <- tryCatch(edge.gof(fp), error = function(err) NULL)
    if (!is.null(e)) { p_ed[b, j] <- e$p_value; s_ed[b, j] <- e$Test_Statistic }
  }
  if (b %% 50 == 0) cat("  rep", b, "/", R,
    sprintf(" (%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

rej    <- function(m) apply(m, 2, function(v) mean(na.omit(v) < alpha))
nvalid <- function(m) apply(m, 2, function(v) sum(!is.na(v)))
r_hl <- rej(p_hl); r_ed <- rej(p_ed)

out <- data.frame(
  lambda        = lam_grid,
  rej_HL        = round(r_hl, 3),
  mcse_HL       = round(sqrt(r_hl * (1 - r_hl) / nvalid(p_hl)), 4),
  rej_EDGE      = round(r_ed, 3),
  mcse_EDGE     = round(sqrt(r_ed * (1 - r_ed) / nvalid(p_ed)), 4),
  mean_stat_EDGE= round(colMeans(s_ed, na.rm = TRUE), 2),
  cal_slope     = round(colMeans(slp, na.rm = TRUE), 3),
  CITL          = round(colMeans(ctl, na.rm = TRUE), 3),
  sd_fitted_pi  = round(colMeans(sdpi, na.rm = TRUE), 4),
  MAE_vs_true   = round(colMeans(mae, na.rm = TRUE), 4),
  n_valid_EDGE  = nvalid(p_ed),
  HL_degenerate = round(colMeans(degen, na.rm = TRUE), 3)
)

cat("\n============ HIGH-DIM SWEEP (n=400, p=100, AR(1) 0.7, H0 TRUE) ============\n")
print(out, row.names = FALSE)
cat(sprintf("\nlambda.min : median %.4f  IQR [%.4f, %.4f]\n",
            median(lam_min), quantile(lam_min,.25), quantile(lam_min,.75)))
cat(sprintf("lambda.1se : median %.4f  IQR [%.4f, %.4f]\n",
            median(lam_1se), quantile(lam_1se,.25), quantile(lam_1se,.75)))
cat(sprintf("P(lambda.1se > 0.20) = %.3f ;  P(lambda.1se > 0.05) = %.3f\n",
            mean(lam_1se > 0.20), mean(lam_1se > 0.05)))
cat("=========================================================================\n")

write.csv(out, "hd_sweep_results.csv", row.names = FALSE)
saveRDS(list(out=out, lam_min=lam_min, lam_1se=lam_1se, lam_grid=lam_grid,
             R=R, n=n, p=p, rho=rho, beta=beta), "hd_sweep_full.rds")
cat("Saved. Runtime:",
    round(as.numeric(difftime(Sys.time(), t0, units="mins")),1), "min\n")
