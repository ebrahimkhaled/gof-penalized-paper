# ROUND 2 (inline): three candidate fixes at the two failing panel-B cells.
#   P1 prepivot-DEB : y* ~ Bern(pi(beta.tilde)); statistic = CORRECTED S_c in
#                     both worlds (refit ridge, re-debias, re-subtract). The
#                     reference inherits S_c's own finite-sample drift.
#   P2 prepivot-PEN : same but generator = pi.hat_K (penalized fit itself).
#   P3 wild         : multiplier bootstrap of the corrected residual
#                     r*(g) = V^{-1/2}C(eps*g) - U F^{-1} X1'(eps*g),
#                     eps = y - pi.hat_K, g ~ Rademacher. No refits.
# Cells: panel B (n=400, p=100, AR1 0.7) lambda in {50, 416}. B=120, NB=149.
suppressPackageStartupMessages(library(glmnet))
set.seed(2026)
G <- 10

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

# corrected statistic S_c and the pieces needed for the wild bootstrap
Sc_of <- function(X1, y, fit, lam, D, need_wild=FALSE) {
  pi <- fit$pi; w <- pmax(pi*(1-pi),1e-8)
  g  <- ceiling(rank(pi, ties.method="first")*G/length(pi))
  Vg <- as.numeric(tapply(w,g,sum))
  r  <- as.numeric(tapply(y-pi,g,sum))/sqrt(Vg)
  U  <- rowsum(w*X1, g)/sqrt(Vg)
  Fm <- crossprod(X1, w*X1)
  K  <- lam*D
  bt <- fit$beta + drop(solve(Fm, K %*% fit$beta))
  mu <- drop(U %*% solve(Fm + K, K %*% bt))
  out <- list(Sc = sum((r-mu)^2), bt = bt)
  if (need_wild) { out$grp <- g; out$Vg <- Vg; out$U <- U; out$Fm <- Fm
                   out$eps <- y - pi; out$X1 <- X1 }
  out
}

wild_ref <- function(st, nb=299) {
  # r*(g) = V^{-1/2} C (eps*g) - U F^{-1} X1'(eps*g), g Rademacher
  n <- length(st$eps)
  FinvXt <- solve(st$Fm, t(st$X1))          # (p1 x n)
  Sst <- numeric(nb)
  for (b in 1:nb) {
    gm <- sample(c(-1,1), n, TRUE)
    e  <- st$eps * gm
    r1 <- as.numeric(tapply(e, st$grp, sum))/sqrt(st$Vg)
    r2 <- drop(st$U %*% (FinvXt %*% e))
    Sst[b] <- sum((r1 - r2)^2)
  }
  Sst
}

n <- 400; p <- 100; rho <- 0.7
Sig <- rho^abs(outer(1:p,1:p,"-")); Ch <- chol(Sig)
beta0 <- rep(c(0.35,-0.30,0.25,-0.20,0.15), length.out=p)*0.8887
D <- diag(c(0, rep(1,p)))

B <- 120; NB <- 149
res <- list()
t0 <- Sys.time()
for (lam in c(50, 416)) {
  p1 <- p2 <- p3 <- numeric(B)
  for (b in 1:B) {
    X <- matrix(rnorm(n*p),n,p)%*%Ch; X1 <- cbind(1,X)
    y <- rbinom(n,1,1/(1+exp(-drop(X%*%beta0))))
    ft <- ridge_fit(X1,y,lam,D)
    st <- Sc_of(X1,y,ft,lam,D,need_wild=TRUE)

    # P3 wild (cheap, do first)
    p3[b] <- (1+sum(wild_ref(st) >= st$Sc))/300

    # generators
    pi_deb <- pmin(pmax(as.numeric(1/(1+exp(-X1 %*% st$bt))),1e-6),1-1e-6)
    pi_pen <- pmin(pmax(ft$pi,1e-6),1-1e-6)
    S1 <- S2 <- numeric(NB)
    for (bb in 1:NB) {
      y1 <- rbinom(n,1,pi_deb)
      S1[bb] <- Sc_of(X1,y1,ridge_fit(X1,y1,lam,D),lam,D)$Sc
      y2 <- rbinom(n,1,pi_pen)
      S2[bb] <- Sc_of(X1,y2,ridge_fit(X1,y2,lam,D),lam,D)$Sc
    }
    p1[b] <- (1+sum(S1 >= st$Sc))/(NB+1)
    p2[b] <- (1+sum(S2 >= st$Sc))/(NB+1)
    if (b%%20==0) cat("lam",lam," rep",b,"/",B,
      " P1:",round(mean(p1[1:b]<0.05),3),
      " P2:",round(mean(p2[1:b]<0.05),3),
      " P3:",round(mean(p3[1:b]<0.05),3),
      " (",round(as.numeric(difftime(Sys.time(),t0,units="mins")),1),"min)\n")
  }
  res[[as.character(lam)]] <- c(P1_prepivot_deb=mean(p1<0.05),
                                P2_prepivot_pen=mean(p2<0.05),
                                P3_wild=mean(p3<0.05))
}
cat("\n===== ROUND 2, panel B, B=",B," (MC SE ~",round(sqrt(.05*.95/B),3),") =====\n")
print(round(do.call(rbind,res),3))
cat("Target: [0.03, 0.08]. Analytic-corrected was 0.026 (lam=50), 0.016 (lam=416).\n")
saveRDS(res,"round2_results.rds")
cat("Runtime:",round(as.numeric(difftime(Sys.time(),t0,units="mins")),1),"min\n")
