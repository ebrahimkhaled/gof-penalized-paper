# =====================================================================
# shrink.gof() -- the shrinkage-corrected Hosmer-Lemeshow test
#
#   SC.HL    the correction on the decile (Hosmer-Lemeshow) basis
#   SC.EDGE  the correction on the EDGE basis
#
# This is the reference implementation for
#
#   Ebrahim, E.K. (2026). Shrinkage invalidates the Hosmer-Lemeshow test:
#   goodness of fit for penalized logistic regression, with an application
#   to glaucoma diagnosis.
#
# It is SELF-CONTAINED: it depends only on base R and stats, so reproducing
# the paper does not depend on any package version. It is the code that will
# be released as shrink.gof() in a future version of the 'ebrahim.gof'
# package; that package does NOT currently export it.
#
# The procedure (Sections 3.1 and 3.2 of the paper):
#   1. fit the ridge model by penalized IRLS;
#   2. form the grouped standardized residuals r;
#   3. subtract the estimated shrinkage non-centrality mu-hat (Proposition 2);
#   4. refer the result to a bootstrap built from the debiased generator
#      pi(beta-tilde) -- Beran prepivoting.
#
# lambda is on the THEORY scale: lambda = n * lambda_glmnet.
# =====================================================================

# ---- penalized IRLS ---------------------------------------------------
.sg_ridge <- function(X1, y, lambda, D, tol = 1e-10, maxit = 100) {
  beta <- rep(0, ncol(X1))
  for (it in seq_len(maxit)) {
    eta <- drop(X1 %*% beta)
    pi  <- 1 / (1 + exp(-eta))
    w   <- pmax(pi * (1 - pi), 1e-8)
    z   <- eta + (y - pi) / w
    bn  <- drop(solve(crossprod(X1, w * X1) + lambda * D,
                      crossprod(X1, w * z)))
    if (max(abs(bn - beta)) < tol) { beta <- bn; break }
    beta <- bn
  }
  list(beta = beta, pi = 1 / (1 + exp(-drop(X1 %*% beta))))
}

# ---- grouped residuals, mu-hat, and both statistics -------------------
.sg_pieces <- function(X1, y, fit, lambda, D, G) {
  pi <- pmin(pmax(fit$pi, 1e-6), 1 - 1e-6)
  w  <- pmax(pi * (1 - pi), 1e-8)
  g  <- pmin(ceiling(rank(pi, ties.method = "first") / (length(y) / G)), G)
  idx <- split(seq_along(y), g)
  Vg <- vapply(idx, function(I) sum(w[I]), 0.0)
  r  <- (vapply(idx, function(I) sum(y[I]),  0.0) -
         vapply(idx, function(I) sum(pi[I]), 0.0)) / sqrt(Vg)
  U  <- t(vapply(idx, function(I) colSums(w[I] * X1[I, , drop = FALSE]),
                 numeric(ncol(X1)))) / sqrt(Vg)
  Fm <- crossprod(X1, w * X1)
  K  <- lambda * D
  # one-step debias, then the shrinkage non-centrality it implies
  bt <- fit$beta + drop(solve(Fm, K %*% fit$beta))
  mu <- drop(U %*% solve(Fm + K, K %*% bt))
  v  <- r - mu                                   # the corrected residual
  pbar <- vapply(idx, function(I) mean(pi[I]), 0.0)
  Z <- tryCatch(as.matrix(stats::poly(pbar, 3)), error = function(e) NULL)
  qf <- function(u, Z) {
    zr <- crossprod(Z, u)
    as.numeric(t(zr) %*% solve(crossprod(Z)) %*% zr)
  }
  list(S_dec_unc = sum(r^2),
       S_dec     = sum(v^2),
       S_edge_unc = if (is.null(Z)) NA_real_ else qf(r, Z),
       S_edge     = if (is.null(Z)) NA_real_ else qf(v, Z),
       Om = diag(G) - U %*% solve(Fm, t(U)),     # the MLE covariance
       Z = Z, beta_tilde = bt)
}

# ---- the uncorrected reference, by direct Monte Carlo -----------------
.sg_mc_p <- function(S, Om, Z = NULL, nd = 20000) {
  if (is.na(S)) return(NA_real_)
  R <- chol(Om + diag(1e-9, nrow(Om)))
  W <- matrix(rnorm(nd * nrow(Om)), nd) %*% R
  if (is.null(Z)) return(mean(rowSums(W^2) >= S))
  A <- solve(crossprod(Z))
  mean(apply(W, 1, function(u) {
    zr <- crossprod(Z, u); as.numeric(t(zr) %*% A %*% zr)
  }) >= S)
}

#' The shrinkage-corrected Hosmer-Lemeshow test
#'
#' @param X   numeric matrix of covariates, WITHOUT an intercept column.
#' @param y   0/1 response.
#' @param lambda ridge penalty on the theory scale, lambda = n * lambda_glmnet.
#' @param G   number of groups (default 10).
#' @param basis "edge" (SC.EDGE, the default and the more powerful),
#'   "decile" (SC.HL), or both.
#' @param B   bootstrap replicates for the prepivoted reference.
#' @param seed optional integer; set it for a reproducible p-value.
#' @param penalize logical vector of length ncol(X): which columns are
#'   penalized. Defaults to all of them; the intercept is never penalized.
#' @param uncorrected if TRUE, also return the uncorrected p-value, which is
#'   the quantity the paper shows to be invalid. For comparison only.
#'
#' @return a list with the statistic and p-value for each basis requested.
shrink.gof <- function(X, y, lambda, G = 10,
                       basis = c("edge", "decile"), B = 999,
                       seed = NULL, penalize = NULL, uncorrected = FALSE) {
  basis <- match.arg(basis, c("edge", "decile"), several.ok = TRUE)
  X <- as.matrix(X); y <- as.numeric(y)
  stopifnot(all(y %in% c(0, 1)), nrow(X) == length(y), lambda >= 0, G >= 3)
  if (!is.null(seed)) set.seed(seed)
  n  <- length(y)
  X1 <- cbind(1, X)
  if (is.null(penalize)) penalize <- rep(TRUE, ncol(X))
  D  <- diag(c(0, as.numeric(penalize)))       # intercept never penalized

  fit <- .sg_ridge(X1, y, lambda, D)
  st  <- .sg_pieces(X1, y, fit, lambda, D, G)
  if ("edge" %in% basis && is.null(st$Z))
    stop("the EDGE basis needs at least 4 distinct group means; raise G")

  # prepivot: resample from the DEBIASED generator pi(beta-tilde), not pi(beta-hat)
  gen <- pmin(pmax(as.numeric(1 / (1 + exp(-X1 %*% st$beta_tilde))), 1e-6), 1 - 1e-6)
  Sd <- Se <- numeric(B)
  for (b in seq_len(B)) {
    ys <- rbinom(n, 1, gen)
    s2 <- .sg_pieces(X1, ys, .sg_ridge(X1, ys, lambda, D), lambda, D, G)
    Sd[b] <- s2$S_dec; Se[b] <- s2$S_edge
  }

  out <- list(lambda = lambda, G = G, B = B, n = n, p = ncol(X))
  if ("decile" %in% basis) {
    out$SC.HL <- list(statistic = st$S_dec,
                      p.value = (1 + sum(Sd >= st$S_dec)) / (B + 1))
    if (uncorrected)
      out$SC.HL$p.uncorrected <- .sg_mc_p(st$S_dec_unc, st$Om)
  }
  if ("edge" %in% basis) {
    out$SC.EDGE <- list(statistic = st$S_edge,
                        p.value = (1 + sum(Se >= st$S_edge, na.rm = TRUE)) / (B + 1))
    if (uncorrected)
      out$SC.EDGE$p.uncorrected <- .sg_mc_p(st$S_edge_unc, st$Om, st$Z)
  }
  class(out) <- "shrink.gof"
  out
}

print.shrink.gof <- function(x, ...) {
  cat("\nShrinkage-corrected Hosmer-Lemeshow test\n")
  cat(sprintf("n = %d, p = %d, lambda = %.4g, G = %d, B = %d\n\n",
              x$n, x$p, x$lambda, x$G, x$B))
  for (nm in c("SC.HL", "SC.EDGE")) if (!is.null(x[[nm]])) {
    cat(sprintf("  %-8s statistic = %8.4f   p = %.4f", nm,
                x[[nm]]$statistic, x[[nm]]$p.value))
    if (!is.null(x[[nm]]$p.uncorrected))
      cat(sprintf("   (uncorrected p = %.5f)", x[[nm]]$p.uncorrected))
    cat("\n")
  }
  cat("\n")
  invisible(x)
}
