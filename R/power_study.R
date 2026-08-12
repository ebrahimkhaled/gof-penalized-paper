# POWER STUDY for the prepivoted corrected test.
# Q1 power CURVE vs effect size (not one point)
# Q2 does the correction cost power vs the standard test on an MLE fit?
# Q3 does it detect anything at all in HIGH DIMENSION (p=100)? -- no evidence yet
# Q4 does the alternative TYPE matter (quadratic vs interaction)?
suppressPackageStartupMessages(library(glmnet))
set.seed(9001)
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
pieces <- function(X1, y, fit, lam, D) {
  pi <- fit$pi; w <- pmax(pi*(1-pi),1e-8)
  g <- ceiling(rank(pi, ties.method="first")*G/length(pi))
  Vg <- as.numeric(tapply(w,g,sum))
  r <- as.numeric(tapply(y-pi,g,sum))/sqrt(Vg)
  U <- rowsum(w*X1,g)/sqrt(Vg); Fm <- crossprod(X1,w*X1); K <- lam*D
  bt <- fit$beta + drop(solve(Fm, K %*% fit$beta))
  mu <- drop(U %*% solve(Fm+K, K %*% bt))
  list(S=sum(r^2), Sc=sum((r-mu)^2), bt=bt,
       Om_MLE=diag(G) - U %*% solve(Fm, t(U)))
}
p_prepivot <- function(X1,y,lam,D) {                # the proposed test
  ft <- ridge_fit(X1,y,lam,D); st <- pieces(X1,y,ft,lam,D)
  gen <- pmin(pmax(as.numeric(1/(1+exp(-X1 %*% st$bt))),1e-6),1-1e-6)
  Ss <- numeric(NB)
  for (bb in 1:NB) { ys <- rbinom(nrow(X1),1,gen)
    Ss[bb] <- pieces(X1,ys,ridge_fit(X1,ys,lam,D),lam,D)$Sc }
  (1+sum(Ss>=st$Sc))/(NB+1)
}
p_mle_std <- function(X1,y,D) {                     # comparator: classical test, MLE fit
  ft <- ridge_fit(X1,y,0,D); st <- pieces(X1,y,ft,0,D)
  R <- chol(st$Om_MLE + diag(1e-9,G)); Z <- matrix(rnorm(ND*G),ND) %*% R
  mean(rowSums(Z^2) >= st$S)
}

B <- 200
## ---------- PANEL A: n=500, p=5, lambda=100 ----------
nA<-500; pA<-5; bA<-c(0.5,-0.4,0.3,-0.3,0.2); DA<-diag(c(0,rep(1,pA)))
altA <- function(X, kind, gam) {
  e <- drop(X %*% bA)
  if (kind=="quad") e + gam*(X[,1]^2 - 1)
  else if (kind=="inter") e + gam*(X[,1]*X[,2])
  else e
}
resA <- NULL
for (kind in c("quad","inter")) for (gam in c(0, 0.25, 0.5, 0.75, 1.0)) {
  if (kind=="inter" && gam==0) next
  hp <- hm <- numeric(B)
  for (b in 1:B) {
    X <- matrix(rnorm(nA*pA),nA,pA); X1 <- cbind(1,X)
    y <- rbinom(nA,1,1/(1+exp(-altA(X,kind,gam))))
    hp[b] <- p_prepivot(X1,y,100,DA) < 0.05
    hm[b] <- p_mle_std(X1,y,DA) < 0.05
  }
  resA <- rbind(resA, data.frame(panel="A", alt=kind, gamma=gam,
                 prepivot=mean(hp), mle_std=mean(hm),
                 se_pre=round(sqrt(mean(hp)*(1-mean(hp))/B),4)))
  cat(sprintf("A %-5s gam=%.2f  prepivot=%.3f  mle_std=%.3f\n", kind, gam, mean(hp), mean(hm)))
}
cat("\n===== PANEL A POWER (n=500, p=5, lambda=100, B=200) =====\n"); print(resA, row.names=FALSE)

## ---------- PANEL B: n=400, p=100, AR(1) 0.7, lambda=416 ----------
## NOTE: no valid comparator here -- the classical test on the MLE fit has size 0.284.
nB<-400; pB<-100; rho<-0.7
Sg <- rho^abs(outer(1:pB,1:pB,"-")); Ch <- chol(Sg)
bB <- rep(c(0.35,-0.30,0.25,-0.20,0.15), length.out=pB)*0.8887
DB <- diag(c(0,rep(1,pB)))
resB <- NULL
for (gam in c(0, 0.5, 1.0, 1.5)) {
  hp <- numeric(B)
  for (b in 1:B) {
    X <- matrix(rnorm(nB*pB),nB,pB) %*% Ch; X1 <- cbind(1,X)
    e <- drop(X %*% bB) + gam*(X[,1]^2 - 1)
    y <- rbinom(nB,1,1/(1+exp(-e)))
    hp[b] <- p_prepivot(X1,y,416,DB) < 0.05
  }
  resB <- rbind(resB, data.frame(panel="B", alt="quad", gamma=gam,
                 prepivot=mean(hp), mle_std=NA,
                 se_pre=round(sqrt(mean(hp)*(1-mean(hp))/B),4)))
  cat(sprintf("B quad  gam=%.2f  prepivot=%.3f\n", gam, mean(hp)))
}
cat("\n===== PANEL B POWER (n=400, p=100, lambda=416, B=200) =====\n"); print(resB, row.names=FALSE)
saveRDS(rbind(resA,resB), "power_study_results.rds")
cat("\nNOTE: gamma=0 rows are SIZE checks under this seed; panel B has no valid comparator\n")
