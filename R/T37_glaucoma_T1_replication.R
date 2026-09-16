# T3.7 -- the paper's own rule, applied to its own new result. T32 returned 0.080 for the EDGE
# basis under truth T1 on the glaucoma design, which is 4.35 null standard errors above nominal and
# has an exact interval of [0.064, 0.099]. Section S4.6 of the supplement states that any cell
# falling outside the pre-specified range must be re-estimated on an independent replication block
# before anything is written about it, so this repeats the T1 arm on fresh seeds (420001+ against
# T32's 400001+). Design, truth, penalty rule, G and NB are exactly T32's.
# Output: T37_glaucoma_T1_pvalues.csv and T37_glaucoma_T1_results.rds
suppressPackageStartupMessages({library(glmnet); library(TH.data); library(parallel)})
data("GlaucomaM", package = "TH.data")
Xr <- as.matrix(GlaucomaM[, setdiff(names(GlaucomaM), "Class")])
yobs <- as.numeric(GlaucomaM$Class == "glaucoma")
X  <- scale(Xr); X <- X[, apply(X, 2, function(c) all(is.finite(c)) & sd(c) > 0), drop = FALSE]
n <- nrow(X); p <- ncol(X); X1 <- cbind(1, X); D <- diag(c(0, rep(1, p)))
G <- 10; NB <- 499; REPS <- 1000

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
  qf <- function(u, Z) { zr <- crossprod(Z, u); as.numeric(t(zr) %*% solve(crossprod(Z)) %*% zr) }
  list(Sc = sum(v^2), Sce = if (is.null(Z)) NA else qf(v, Z), bt = bt)
}

set.seed(1); cv0 <- cv.glmnet(X, yobs, family = "binomial", alpha = 0)
PI0 <- ridge_fit(X1, yobs, cv0$lambda.min*n, D)$pi      # truth T1, exactly as in T32
cat(sprintf("truth T1: probabilities in [%.2e, %.4f], sd %.3f, expected events %.1f\n",
            min(PI0), max(PI0), sd(PI0), sum(PI0)))

worker <- function(seed) {
  set.seed(seed)
  y <- rbinom(n, 1, PI0)
  if (sum(y) < 5 || sum(y) > n - 5) return(c(seed = seed, lam = NA, p_dec = NA, p_edge = NA))
  cv <- try(cv.glmnet(X, y, family = "binomial", alpha = 0), silent = TRUE)
  if (inherits(cv, "try-error")) return(c(seed = seed, lam = NA, p_dec = NA, p_edge = NA))
  lam <- cv$lambda.1se*n
  ft <- ridge_fit(X1, y, lam, D); st <- pieces(X1, y, ft, lam, D, G)
  gen <- pmin(pmax(as.numeric(1/(1 + exp(-X1 %*% st$bt))), 1e-6), 1 - 1e-6)
  Sd <- Se <- numeric(NB)
  for (b in 1:NB) { ys <- rbinom(n, 1, gen)
    s2 <- pieces(X1, ys, ridge_fit(X1, ys, lam, D), lam, D, G); Sd[b] <- s2$Sc; Se[b] <- s2$Sce }
  c(seed = seed, lam = lam, p_dec = (1 + sum(Sd >= st$Sc))/(NB + 1),
    p_edge = (1 + sum(Se >= st$Sce, na.rm = TRUE))/(NB + 1))
}

nc <- max(1, detectCores() - 2)
cl <- makeCluster(nc); on.exit(stopCluster(cl))
clusterEvalQ(cl, suppressPackageStartupMessages(library(glmnet)))
clusterExport(cl, c("X", "X1", "D", "n", "G", "NB", "PI0", "ridge_fit", "pieces"))
clusterSetRNGStream(cl, 20260919)
cat("T3.7 independent replication of the T1 arm on", nc, "workers, NB =", NB, "\n")

t0 <- Sys.time()
R <- as.data.frame(t(parSapply(cl, 420000L + seq_len(REPS), worker)))
write.csv(R, "T37_glaucoma_T1_pvalues.csv", row.names = FALSE)
saveRDS(R, "T37_glaucoma_T1_results.rds")
ok <- !is.na(R$p_edge)
for (b in c("p_dec", "p_edge")) {
  x <- R[[b]][ok]; r <- mean(x < 0.05); ci <- binom.test(sum(x < 0.05), length(x))$conf.int
  cat(sprintf("%-7s %d reps: rate %.4f  se %.4f  exact 95%% [%.4f, %.4f]  z vs 0.05 %+.2f\n",
              b, length(x), r, sqrt(r*(1-r)/length(x)), ci[1], ci[2],
              (r - 0.05)/sqrt(0.05*0.95/length(x))))
}
cat(sprintf("T32 reported 0.042 (decile) and 0.080 (EDGE) on seeds 400001+\n"))
cat("dropped replicates:", sum(!ok), "| runtime:",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
