# H6: the fragility discussion says the procedure depends on F^{-1}, and the simulation
# measures the overshoot as ||beta.tilde||/||beta_0||. On REAL data beta_0 does not exist,
# so the analogous quantity is ||beta.tilde||/||beta.hat||, reported alongside cond(F).
# Registered fits: glaucoma_deep.R, set.seed(1), lambda in {19.9, 147.65, 295.3, 590.6}.
suppressPackageStartupMessages({library(glmnet); library(TH.data)})
data("GlaucomaM", package = "TH.data")
X <- scale(as.matrix(GlaucomaM[, setdiff(names(GlaucomaM), "Class")]))
X <- X[, apply(X, 2, function(c) all(is.finite(c)) & sd(c) > 0), drop = FALSE]
y <- as.numeric(GlaucomaM$Class == "glaucoma"); n <- nrow(X); p <- ncol(X)
X1 <- cbind(1, X); D <- diag(c(0, rep(1, p)))

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

set.seed(1); cv <- cv.glmnet(X, y, family = "binomial", alpha = 0)
lams <- c(lambda.min = cv$lambda.min*n, half.1se = cv$lambda.1se*n/2,
          lambda.1se = cv$lambda.1se*n, twice.1se = 2*cv$lambda.1se*n)
out <- NULL
for (nm in names(lams)) {
  lam <- lams[[nm]]
  ft <- ridge_fit(X1, y, lam, D)
  pi <- pmin(pmax(ft$pi, 1e-8), 1-1e-8); w <- pmax(pi*(1-pi), 1e-10)
  F <- crossprod(X1, w*X1)                      # F = X' W-hat X at the ridge fit
  bt <- ft$beta + drop(solve(F, lam*D %*% ft$beta))
  ev <- eigen(F, symmetric = TRUE, only.values = TRUE)$values
  out <- rbind(out, data.frame(
    lambda = nm, lam = round(lam, 1),
    condF  = max(ev)/min(ev),
    nbh    = sqrt(sum(ft$beta^2)),
    nbt    = sqrt(sum(bt^2)),
    ratio  = sqrt(sum(bt^2))/sqrt(sum(ft$beta^2))))
}
rownames(out) <- NULL
out$condF_fmt <- sprintf("%.2e", out$condF)
print(out[, c("lambda","lam","condF_fmt","nbh","nbt","ratio")], digits = 4, row.names = FALSE)
cat(sprintf("\ncond(F) range: %.2e to %.2e\n", min(out$condF), max(out$condF)))
cat(sprintf("||beta.tilde||/||beta.hat|| range: %.1f to %.1f\n", min(out$ratio), max(out$ratio)))
write.csv(out, "glaucoma_conditioning.csv", row.names = FALSE)
cat("wrote glaucoma_conditioning.csv\n")
