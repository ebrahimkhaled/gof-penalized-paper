# T43 (diagnostic, not a pre-declared study): does the bootstrap generator pi(beta_tilde) over-state the signal
# at large p/n, as the Sur-Candes inflation of the unpenalized MLE would predict? The debiased beta_tilde is the
# first-order MLE, so in the proportional regime its linear predictor should be too spread out.
# For null datasets of T40 cells (same generator and seeds as T40), report the ratio of sd(X beta_tilde) and
# sd(X beta_hat) to the true sd(eta0) = 1.5. Uses T40's own functions, sourced without running the driver.
suppressPackageStartupMessages(library(glmnet))
src <- readLines("T40_highdim_grid.R")
eval(parse(text = src[seq_len(grep("^# ---- driver", src) - 1)]))   # everything above the driver section

R_PER_CELL <- 40L
out <- do.call(rbind, lapply(c(1L, 15L, 22L, 29L, 44L, 52L, 64L, 72L, 84L, 88L, 92L), function(cc) {
  ce <- cells[cells$cell == cc, ]
  do.call(rbind, lapply(seq_len(R_PER_CELL), function(r) {
    d <- make_data(ce, r); X <- d$X; y <- d$y; p <- ncol(X)
    X1 <- cbind(1, X); D <- diag(c(0, rep(1, p)))
    cv <- cv.glmnet(X, y, family = "binomial", alpha = 0, nfolds = 10)
    lam <- cv$lambda.1se * N
    fit <- ridge_fit(X1, y, lam, D); st <- pieces(X1, y, fit, lam, D)
    data.frame(cell = cc, panel = ce$panel, rho = ce$rho, p = p, rep = r,
               ratio_tilde = sd(drop(X1 %*% st$bt)) / GAMMA, ratio_hat = sd(drop(X1 %*% fit$beta)) / GAMMA)
  }))
}))
write.csv(out, file.path("..", "data", "T43_generator_inflation.csv"), row.names = FALSE)
print(aggregate(cbind(ratio_tilde, ratio_hat) ~ panel + rho + p, data = out, FUN = function(x) round(median(x), 2)))
