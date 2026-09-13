# CONFIRMATION: P1 "prepivoted parametric bootstrap at the debiased fit".
# Procedure (one rule, applied everywhere):
#   1. Fit ridge at lambda; compute the corrected statistic S_c = ||r - mu.hat||^2.
#   2. beta.tilde = one-step debias; generator pi(beta.tilde).
#   3. NB times: y* ~ Bern(pi(beta.tilde)); refit ridge at same lambda;
#      recompute the SAME corrected statistic S_c* in the bootstrap world.
#   4. p = (1 + #{S_c* >= S_c}) / (NB + 1).
# SIZE: six cells, fresh seed 4051, B=300, NB=149.
# POWER: panel A omitted quadratic 0.5*(x1^2-1), lambda in {25, 100}, B=200.
suppressPackageStartupMessages(library(glmnet))
set.seed(4051)
G <- 10; NB <- 149

ridge_fit <- function(X1, y, lam, D, tol=1e-9, maxit=60) {
  beta <- rep(0, ncol(X1))
  for (it in 1:maxit) {
    eta <- drop(X1 %*% beta); pi <- 1/(1+exp(-eta))
    w <- pmax(pi*(1-pi), 1e-8); z <- eta + (y-pi)/w
    bn <- drop(solve(crossprod(X1, w*X1) + lam*D, crossprod(X1, w*z)))
    if (max(abs(bn-beta)) < tol) { beta <- bn; break }
    beta <- bn
  }
  eta <- drop(X1 %*% beta)
  list(beta=beta, pi=1/(1+exp(-eta)))
}
Sc_of <- function(X1, y, fit, lam, D) {
  pi <- fit$pi; w <- pmax(pi*(1-pi),1e-8)
  g  <- ceiling(rank(pi, ties.method="first")*G/length(pi))
  Vg <- as.numeric(tapply(w,g,sum))
  r  <- as.numeric(tapply(y-pi,g,sum))/sqrt(Vg)
  U  <- rowsum(w*X1, g)/sqrt(Vg)
  Fm <- crossprod(X1, w*X1); K <- lam*D
  bt <- fit$beta + drop(solve(Fm, K %*% fit$beta))
  mu <- drop(U %*% solve(Fm + K, K %*% bt))
  list(Sc = sum((r-mu)^2), bt = bt)
}
p1_test <- function(X1, y, lam, D) {
  ft <- ridge_fit(X1,y,lam,D)
  st <- Sc_of(X1,y,ft,lam,D)
  gen <- pmin(pmax(as.numeric(1/(1+exp(-X1 %*% st$bt))),1e-6),1-1e-6)
  Ss <- numeric(NB)
  for (bb in 1:NB) {
    ys <- rbinom(nrow(X1),1,gen)
    Ss[bb] <- Sc_of(X1,ys,ridge_fit(X1,ys,lam,D),lam,D)$Sc
  }
  (1+sum(Ss >= st$Sc))/(NB+1)
}

## panel definitions
nA <- 500; pA <- 5; betaA <- c(0.5,-0.4,0.3,-0.3,0.2); DA <- diag(c(0,rep(1,pA)))
nB_ <- 400; pB <- 100; rho <- 0.7
SigB <- rho^abs(outer(1:pB,1:pB,"-")); ChB <- chol(SigB)
betaB <- rep(c(0.35,-0.30,0.25,-0.20,0.15), length.out=pB)*0.8887
DB <- diag(c(0,rep(1,pB)))

B <- 300
cells <- list(list(p="A",lam=100), list(p="A",lam=137), list(p="A",lam=200),
              list(p="B",lam=50),  list(p="B",lam=416), list(p="B",lam=1000))
t0 <- Sys.time(); out <- NULL
for (cl in cells) {
  pv <- numeric(B)
  for (b in 1:B) {
    if (cl$p=="A") { X <- matrix(rnorm(nA*pA),nA,pA); X1 <- cbind(1,X)
      y <- rbinom(nA,1,1/(1+exp(-drop(X%*%betaA)))); pv[b] <- p1_test(X1,y,cl$lam,DA)
    } else { X <- matrix(rnorm(nB_*pB),nB_,pB)%*%ChB; X1 <- cbind(1,X)
      y <- rbinom(nB_,1,1/(1+exp(-drop(X%*%betaB)))); pv[b] <- p1_test(X1,y,cl$lam,DB) }
  }
  rej <- mean(pv<0.05)
  out <- rbind(out, data.frame(panel=cl$p, lambda=cl$lam, rejection=round(rej,3),
                               mcse=round(sqrt(rej*(1-rej)/B),4), B=B))
  cat(sprintf("cell %s lam=%d done: rej=%.3f (%.1f min)\n", cl$p, cl$lam, rej,
              as.numeric(difftime(Sys.time(),t0,units="mins"))))
}
cat("\n===== P1 CONFIRMATION SIZE (target [0.03,0.08]) =====\n"); print(out, row.names=FALSE)

## POWER
Bp <- 200; pw <- setNames(numeric(2), c("25","100"))
for (lam in c(25,100)) {
  hit <- numeric(Bp)
  for (b in 1:Bp) {
    X <- matrix(rnorm(nA*pA),nA,pA); X1 <- cbind(1,X)
    eta <- drop(X%*%betaA) + 0.5*(X[,1]^2 - 1)
    y <- rbinom(nA,1,1/(1+exp(-eta)))
    hit[b] <- p1_test(X1,y,lam,DA) < 0.05
  }
  pw[as.character(lam)] <- mean(hit)
  cat(sprintf("power lam=%d: %.3f\n", lam, mean(hit)))
}
cat("\n===== P1 POWER (omitted quadratic, must be > 0.10 at lam=25) =====\n"); print(round(pw,3))
saveRDS(list(size=out, power=pw), "confirm_P1_results.rds")
cat("TOTAL:", round(as.numeric(difftime(Sys.time(),t0,units="mins")),1), "min\n")
