# T3.4 -- BJ revision, co-author comment 4: size of the prepivoted corrected test under a smooth,
# convex, non-quadratic penalty (log-cosh; definition and implementation checks in
# T34a_logcosh_check.R). delta is set so that the median |beta.hat_j| / delta is about 2 at the
# penalty used, so most coefficients sit where the penalty acts like a lasso rather than a ridge:
#   design A: lambda = 137, delta = 0.1       design B: lambda = 416, delta = 0.01
# The correction replaces K b by the penalty gradient a(b):
#   beta.tilde = beta.hat + F^{-1} a(beta.hat),  mu.hat = U M^{-1} a(beta.tilde),  M = F + diag(K(beta.hat)).
# The bootstrap refits the same penalty at the same lambda and delta. NB = 499. Also reported: the
# uncorrected tests referred to the maximum likelihood covariance by Monte Carlo (4000 draws).
# Output: T34_logcosh_pvalues.csv and T34_logcosh_results.rds.
suppressPackageStartupMessages(library(parallel))
G <- 10; NB <- 499; BA <- 2000; BB <- 1000

worker <- function(seed, panel, lam, delta) {
  set.seed(seed)
  lcosh <- function(x) { a <- abs(x); ifelse(a < 1, log1p(2*sinh(a/2)^2), a + log1p(exp(-2*a)) - log(2)) }
  lc_fit <- function(X1, y, pen, tol=1e-9, maxit=200) {
    obj <- function(b) { eta <- drop(X1%*%b)
      sum(y*eta - (pmax(eta,0) + log1p(exp(-abs(eta))))) - lam*delta^2*sum(pen*lcosh(b/delta)) }
    beta <- rep(0, ncol(X1)); conv <- FALSE
    for (it in 1:maxit) {
      eta <- drop(X1%*%beta); pi <- 1/(1+exp(-eta)); w <- pmax(pi*(1-pi), 1e-8)
      g <- drop(crossprod(X1, y-pi)) - pen*lam*delta*tanh(beta/delta)
      H <- crossprod(X1, w*X1) + diag(pen*lam/cosh(beta/delta)^2)
      step <- drop(solve(H, g))
      if (max(abs(step)) < tol) { conv <- TRUE; break }
      f0 <- obj(beta); t <- 1
      for (h in 1:30) { if (obj(beta + t*step) >= f0 - 1e-10*max(1, abs(f0))) break; t <- t/2 }
      beta <- beta + t*step
    }
    list(beta=beta, pi=1/(1+exp(-drop(X1%*%beta))), conv=conv)
  }
  pieces <- function(X1, y, fit, pen) {
    pi<-pmin(pmax(fit$pi,1e-6),1-1e-6); w<-pmax(pi*(1-pi),1e-8)
    g<-pmin(ceiling(rank(pi,ties.method="first")/(length(y)/G)),G)
    idx<-split(seq_along(y),g); Vg<-vapply(idx,function(I) sum(w[I]),0.0)
    r<-(vapply(idx,function(I) sum(y[I]),0.0)-vapply(idx,function(I) sum(pi[I]),0.0))/sqrt(Vg)
    U<-t(vapply(idx,function(I) colSums(w[I]*X1[I,,drop=FALSE]),numeric(ncol(X1))))/sqrt(Vg)
    Fm<-crossprod(X1,w*X1)
    bt<-fit$beta+drop(solve(Fm, pen*lam*delta*tanh(fit$beta/delta)))
    mu<-drop(U%*%solve(Fm+diag(pen*lam/cosh(fit$beta/delta)^2), pen*lam*delta*tanh(bt/delta)))
    pbar<-vapply(idx,function(I) mean(pi[I]),0.0); v<-r-mu
    Z<-tryCatch(as.matrix(stats::poly(pbar,3)),error=function(e) NULL)
    qf<-function(u,Z){zr<-crossprod(Z,u); as.numeric(t(zr)%*%solve(crossprod(Z))%*%zr)}
    list(Sp=sum(r^2), Sc=sum(v^2), Spe=if(is.null(Z)) NA else qf(r,Z), Sce=if(is.null(Z)) NA else qf(v,Z),
         Om=diag(G)-U%*%solve(Fm,t(U)), Z=Z, bt=bt)
  }
  mc_p <- function(S,Om,Z=NULL,nd=4000){
    R<-chol(Om+diag(1e-9,nrow(Om))); W<-matrix(rnorm(nd*nrow(Om)),nd)%*%R
    if(is.null(Z)) return(mean(rowSums(W^2)>=S))
    A<-solve(crossprod(Z)); mean(apply(W,1,function(u){zr<-crossprod(Z,u); as.numeric(t(zr)%*%A%*%zr)})>=S)
  }
  if (panel=="A") {
    n<-500; p<-5; b0<-c(0.5,-0.4,0.3,-0.3,0.2)
    X<-matrix(rnorm(n*p),n,p); X1<-cbind(1,X); y<-rbinom(n,1,1/(1+exp(-drop(X%*%b0))))
  } else {
    n<-400; p<-100; S<-0.7^abs(outer(1:p,1:p,"-")); Ch<-chol(S)
    b0<-rep(c(.35,-.30,.25,-.20,.15),length.out=p)*0.8887
    X<-matrix(rnorm(n*p),n,p)%*%Ch; X1<-cbind(1,X); y<-rbinom(n,1,1/(1+exp(-drop(X%*%b0))))
  }
  pen <- c(0, rep(1, p))
  ft <- lc_fit(X1, y, pen); st <- pieces(X1, y, ft, pen)
  gen <- pmin(pmax(as.numeric(1/(1+exp(-X1%*%st$bt))),1e-6),1-1e-6)
  Sd <- Se <- numeric(NB); nconv <- 0
  for (bb in 1:NB) { ys <- rbinom(n,1,gen); fb <- lc_fit(X1, ys, pen); nconv <- nconv + !fb$conv
    s2 <- pieces(X1, ys, fb, pen); Sd[bb] <- s2$Sc; Se[bb] <- s2$Sce }
  c(seed=seed, dec=(1+sum(Sd>=st$Sc))/(NB+1), edge=(1+sum(Se>=st$Sce,na.rm=TRUE))/(NB+1),
    pu_dec=mc_p(st$Sp,st$Om), pu_edge=if(is.null(st$Z)) NA_real_ else mc_p(st$Spe,st$Om,st$Z),
    ratio=median(abs(ft$beta[-1]))/delta, conv=ft$conv, boot_nonconv=nconv)
}

nc <- max(1, detectCores() - 2)
cl <- makeCluster(nc); on.exit(stopCluster(cl))
clusterExport(cl, c("G","NB"))
clusterSetRNGStream(cl, 20260918)          # RNG kind as in the other studies; each replicate is set.seed()
cat("T3.4 log-cosh penalty on", nc, "workers, NB =", NB, "\n")

cells <- list(list(p="A", lam=137, delta=0.1,  seeds=500000L + seq_len(BA)),
              list(p="B", lam=416, delta=0.01, seeds=510000L + seq_len(BB)))
all <- NULL; summ <- NULL; t0 <- Sys.time()
rate <- function(x) { x <- x[!is.na(x)]; m <- mean(x < 0.05); c(round(m,4), round(sqrt(m*(1-m)/length(x)),4)) }
for (ce in cells) {
  R <- t(parSapply(cl, ce$seeds, worker, panel=ce$p, lam=ce$lam, delta=ce$delta))
  d <- data.frame(panel=ce$p, lambda=ce$lam, delta=ce$delta, R); all <- rbind(all, d)
  write.csv(all, "T34_logcosh_pvalues.csv", row.names=FALSE)
  s <- data.frame(panel=ce$p, lambda=ce$lam, delta=ce$delta, B=nrow(d),
    corr_dec=rate(d$dec)[1], se_dec=rate(d$dec)[2], corr_edge=rate(d$edge)[1], se_edge=rate(d$edge)[2],
    unc_dec=rate(d$pu_dec)[1], unc_edge=rate(d$pu_edge)[1],
    ratio_median=round(median(d$ratio),2), fit_nonconv=sum(d$conv==0), boot_nonconv=sum(d$boot_nonconv))
  summ <- rbind(summ, s); print(s, row.names=FALSE)
  cat(sprintf("  %.0f min\n", as.numeric(difftime(Sys.time(),t0,units="mins"))))
  saveRDS(list(pvalues=all, summary=summ), "T34_logcosh_results.rds")
}
cat("\n===== T3.4 SIZE UNDER THE LOG-COSH PENALTY =====\n"); print(summ, row.names=FALSE)
cat("runtime:", round(as.numeric(difftime(Sys.time(),t0,units="mins")),1), "min\n")
