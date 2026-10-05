# T44 stage-1 check: can the bootstrap recover the true signal strength sd(eta0) = 1.5 from one dataset?
# Model the debiased fit as eta_tilde = a * eta + noise. In the bootstrap world (truth = eta_tilde) the same
# relation gives a = cov(eta*, eta_tilde) / var(eta_tilde) and noise variance nu2 = var(eta* - a eta_tilde),
# so the signal sd is gamma_hat = sqrt(var(eta_tilde) - nu2) / a. Report gamma_hat / 1.5 against the naive
# sd(eta_tilde) / 1.5 (T43), on the cells where the test fails and on one where it does not.
suppressPackageStartupMessages(library(glmnet))
src <- readLines("T44_rescaled_generator_pilot.R")
eval(parse(text = src[seq_len(grep("^rescaled_rep <- function", src) - 1)]))
K0 <- 50L
res <- do.call(rbind, lapply(c(15L, 72L, 92L, 205L, 209L), function(cc) {
  ce <- spec[spec$cell == cc, ]
  do.call(rbind, lapply(1:12, function(r) {
    d <- make_data_a(ce, r); X <- d$X; y <- d$y; p <- ncol(X)
    X1 <- cbind(1, X); D <- diag(c(0, rep(1, p)))
    lam <- cv.glmnet(X, y, family = "binomial", alpha = 0, nfolds = 10)$lambda.1se * N
    st <- pieces(X1, y, ridge_fit(X1, y, lam, D), lam, D)
    et <- drop(X %*% st$bt[-1]); gen0 <- pmin(pmax(plogis(drop(X1 %*% st$bt)), 1e-6), 1 - 1e-6)
    E <- vapply(seq_len(K0), function(k) {
      ys <- rbinom(N, 1, gen0); drop(X %*% pieces(X1, ys, ridge_fit(X1, ys, lam, D), lam, D)$bt[-1]) },
      numeric(N))
    a <- mean(apply(E, 2, function(e) cov(e, et))) / var(et)
    nu2 <- mean(apply(E, 2, function(e) var(e - a * et)))
    g <- sqrt(max(var(et) - nu2, 0)) / a
    data.frame(cell = cc, p = p, naive = sd(et) / GAMMA, a = a, gamma_hat = g / GAMMA)
  }))
}))
print(aggregate(cbind(naive, a, gamma_hat) ~ cell + p, data = res, FUN = function(x) round(median(x), 3)))
print(aggregate(gamma_hat ~ cell, data = res, FUN = function(x) round(sd(x), 3)))
