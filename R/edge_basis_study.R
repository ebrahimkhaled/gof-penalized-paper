# EDGE-poly3 vs decile basis, under the SAME prepivoted correction.
# Paired design: both statistics computed from the SAME ridge fits and the SAME
# bootstrap refits, so the comparison is free of between-run noise.
#
# Basis definitions follow ebrahim.gof::def.gof exactly:
#   groups = equal-frequency deciles of pi.hat; r = (O_g - E_g)/sqrt(V_g)
#   decile : S = ||r||^2                       (10 directions)
#   EDGE   : Z = poly(pbar, 3) (10x3);  S = (Z'r)'(Z'Z)^{-1}(Z'r)
# Correction and prepivoting are identical for both.
suppressPackageStartupMessages(library(glmnet))
set.seed(7302)
G <- 10; NB <- 149; ND <- 4000

ridge_fit <- function(X1, y, lam, D, tol=1e-9, maxit=60) {
  beta <- rep(0, ncol(X1))
  for (it in 1:maxit) {
    eta <- drop(X1 %*% beta); pi <- 1/(1+exp(-eta))
    w <- pmax(pi*(1-pi), 1e-8); z <- eta + (y-pi)/w
    bn <- drop(solve(crossprod(X1, w*X1) + lam*D, crossprod(X1, w*z)))
    if (max(abs(bn-beta)) < tol) { beta <- bn; break }
    beta <- bn
  }
  eta <- drop(X1 %*% beta); list(beta=beta, pi=1/(1+exp(-eta)))
}

# returns BOTH statistics (corrected) + the pieces needed for a MC reference
pieces <- function(X1, y, fit, lam, D) {
  pi <- pmin(pmax(fit$pi, 1e-6), 1-1e-6); w <- pmax(pi*(1-pi),1e-8)
  g  <- pmin(ceiling(rank(pi, ties.method="first")/(length(y)/G)), G)
  idx <- split(seq_along(y), g)
  Vg  <- vapply(idx, function(I) sum(w[I]), 0.0)
  og  <- vapply(idx, function(I) sum(y[I]), 0.0)
  eg  <- vapply(idx, function(I) sum(pi[I]), 0.0)
  pbar<- vapply(idx, function(I) mean(pi[I]), 0.0)
  r   <- (og - eg)/sqrt(Vg)
  U   <- t(vapply(idx, function(I) colSums(w[I]*X1[I,,drop=FALSE]), numeric(ncol(X1))))/sqrt(Vg)
  Fm  <- crossprod(X1, w*X1); K <- lam*D
  bt  <- fit$beta + drop(solve(Fm, K %*% fit$beta))
  mu  <- drop(U %*% solve(Fm+K, K %*% bt))
  rc  <- r - mu
  # EDGE basis
  Z <- tryCatch(as.matrix(stats::poly(pbar, 3)), error=function(e) NULL)
  qf <- function(v, Z) { Zr <- crossprod(Z, v); as.numeric(t(Zr) %*% solve(crossprod(Z)) %*% Zr) }
  list(dec_c = sum(rc^2), edge_c = if (is.null(Z)) NA_real_ else qf(rc, Z),
       dec_p = sum(r^2),  edge_p = if (is.null(Z)) NA_real_ else qf(r, Z),
       Z = Z, bt = bt, Om = diag(G) - U %*% solve(Fm, t(U)))
}

# prepivoted p-values for BOTH bases from ONE set of refits
p_prepivot2 <- function(X1, y, lam, D) {
  ft <- ridge_fit(X1,y,lam,D); st <- pieces(X1,y,ft,lam,D)
  gen <- pmin(pmax(as.numeric(1/(1+exp(-X1 %*% st$bt))),1e-6),1-1e-6)
  Sd <- Se <- numeric(NB)
  for (bb in 1:NB) {
    ys <- rbinom(nrow(X1),1,gen)
    s2 <- pieces(X1,ys,ridge_fit(X1,ys,lam,D),lam,D)
    Sd[bb] <- s2$dec_c; Se[bb] <- s2$edge_c
  }
  c(dec = (1+sum(Sd>=st$dec_c))/(NB+1),
    edge= (1+sum(Se>=st$edge_c, na.rm=TRUE))/(NB+1))
}
# classical comparator on the MLE fit, both bases, MC reference
p_mle2 <- function(X1, y, D) {
  ft <- ridge_fit(X1,y,0,D); st <- pieces(X1,y,ft,0,D)
  R <- chol(st$Om + diag(1e-9,G)); Zr <- matrix(rnorm(ND*G),ND) %*% R
  pd <- mean(rowSums(Zr^2) >= st$dec_p)
  Zb <- st$Z; A <- solve(crossprod(Zb))
  se <- apply(Zr, 1, function(v){ zz <- crossprod(Zb, v); as.numeric(t(zz)%*%A%*%zz) })
  c(dec = pd, edge = mean(se >= st$edge_p))
}

B <- 200; out <- NULL
run_cell <- function(panel, gam, lam, gen_fn, D, do_mle) {
  hp <- matrix(NA, B, 2); hm <- matrix(NA, B, 2)
  for (b in 1:B) {
    d <- gen_fn(gam)
    hp[b,] <- p_prepivot2(d$X1, d$y, lam, D) < 0.05
    if (do_mle) hm[b,] <- p_mle2(d$X1, d$y, D) < 0.05
  }
  r <- data.frame(panel=panel, gamma=gam,
        pre_dec=mean(hp[,1]), pre_edge=mean(hp[,2]),
        mle_dec=if(do_mle) mean(hm[,1]) else NA,
        mle_edge=if(do_mle) mean(hm[,2]) else NA,
        se=round(sqrt(0.25/B),4))
  cat(sprintf("%s gam=%.2f | prepivot dec=%.3f EDGE=%.3f | MLE dec=%s EDGE=%s\n",
      panel, gam, r$pre_dec, r$pre_edge,
      ifelse(do_mle, sprintf("%.3f", r$mle_dec), "--"),
      ifelse(do_mle, sprintf("%.3f", r$mle_edge), "--")))
  r
}

## PANEL A: n=500, p=5, lambda=100
nA<-500; pA<-5; bA<-c(0.5,-0.4,0.3,-0.3,0.2); DA<-diag(c(0,rep(1,pA)))
genA <- function(gam) { X <- matrix(rnorm(nA*pA),nA,pA)
  list(X1=cbind(1,X), y=rbinom(nA,1,1/(1+exp(-(drop(X%*%bA)+gam*(X[,1]^2-1)))))) }
for (gam in c(0, 0.25, 0.5, 0.75, 1.0))
  out <- rbind(out, run_cell("A", gam, 100, genA, DA, TRUE))

## PANEL B: n=400, p=100, AR(1) 0.7, lambda=416
nB<-400; pB<-100; Sg <- 0.7^abs(outer(1:pB,1:pB,"-")); Ch <- chol(Sg)
bB <- rep(c(0.35,-0.30,0.25,-0.20,0.15), length.out=pB)*0.8887; DB <- diag(c(0,rep(1,pB)))
genB <- function(gam) { X <- matrix(rnorm(nB*pB),nB,pB) %*% Ch
  list(X1=cbind(1,X), y=rbinom(nB,1,1/(1+exp(-(drop(X%*%bB)+gam*(X[,1]^2-1)))))) }
for (gam in c(0, 0.5, 1.0, 1.5))
  out <- rbind(out, run_cell("B", gam, 416, genB, DB, FALSE))

cat("\n===== EDGE-poly3 vs DECILE, prepivoted correction, B=200 =====\n")
print(out, row.names=FALSE)
cat("\ngamma=0 rows are SIZE (target 0.05). Panel B: no valid MLE comparator (size 0.284).\n")
saveRDS(out, "edge_basis_results.rds")
