# ORACLE-INDEX EXPERIMENT (ranked #1). Panel B, n=400, p=100, AR(1) 0.7, lambda=416.
#
# QUESTION: is the high-dimensional power loss caused by having to ESTIMATE the
# index (=> rho^{2k} attenuation), or is it something else?
#
# DESIGN: for every replicate compute the prepivoted corrected test TWICE from the
# SAME fit and the SAME bootstrap refits, differing only in the grouping variable:
#    FITTED grouping : deciles of pi.hat_K            (what a practitioner can do)
#    ORACLE grouping : deciles of the TRUE index eta0 (infeasible; isolates the noise)
# In the bootstrap world the generator's truth is beta.tilde, so its "true index"
# is X1 %*% beta.tilde -- that is the consistent oracle analogue there.
#
# PREDICTION (verdict sec.6 item 1):
#    cubic (index-aligned, k=3): oracle power ~0.40 vs fitted ~0.05  -> attenuation
#    quad  (coordinate-wise)   : BOTH stay ~0.05                    -> dilution, not attenuation
# A confirmed prediction converts "the test has no power" into "we lose exactly rho^{2k}".
suppressPackageStartupMessages(library(glmnet))
set.seed(5150)
G <- 10; NB <- 149; B <- 400

ridge_fit <- function(X1,y,lam,D,tol=1e-9,maxit=60){
  beta <- rep(0,ncol(X1))
  for(it in 1:maxit){ eta<-drop(X1%*%beta); pi<-1/(1+exp(-eta))
    w<-pmax(pi*(1-pi),1e-8); z<-eta+(y-pi)/w
    bn<-drop(solve(crossprod(X1,w*X1)+lam*D, crossprod(X1,w*z)))
    if(max(abs(bn-beta))<tol){beta<-bn;break}; beta<-bn }
  list(beta=beta, pi=1/(1+exp(-drop(X1%*%beta))))
}

# corrected statistic under a SUPPLIED grouping vector g
stat_g <- function(X1, y, fit, lam, D, g) {
  pi <- pmin(pmax(fit$pi,1e-6),1-1e-6); w <- pmax(pi*(1-pi),1e-8)
  idx <- split(seq_along(y), g)
  Vg <- vapply(idx,function(I) sum(w[I]),0.0)
  r  <- (vapply(idx,function(I) sum(y[I]),0.0) - vapply(idx,function(I) sum(pi[I]),0.0))/sqrt(Vg)
  U  <- t(vapply(idx,function(I) colSums(w[I]*X1[I,,drop=FALSE]), numeric(ncol(X1))))/sqrt(Vg)
  Fm <- crossprod(X1,w*X1); K <- lam*D
  bt <- fit$beta + drop(solve(Fm, K %*% fit$beta))
  mu <- drop(U %*% solve(Fm+K, K %*% bt))
  pbar <- vapply(idx,function(I) mean(pi[I]),0.0)
  v <- r - mu
  Z <- tryCatch(as.matrix(stats::poly(pbar,3)), error=function(e) NULL)
  edge <- if(is.null(Z)) NA_real_ else { Zr<-crossprod(Z,v); as.numeric(t(Zr)%*%solve(crossprod(Z))%*%Zr) }
  list(dec=sum(v^2), edge=edge, bt=bt)
}
deciles <- function(x) pmin(ceiling(rank(x, ties.method="first")/(length(x)/G)), G)

# one replicate: BOTH groupings, ONE set of refits
one_rep <- function(X1, y, eta0, lam, D) {
  ft <- ridge_fit(X1,y,lam,D)
  gF <- deciles(ft$pi)      # fitted grouping
  gO <- deciles(eta0)       # oracle grouping (true index)
  sF <- stat_g(X1,y,ft,lam,D,gF)
  sO <- stat_g(X1,y,ft,lam,D,gO)
  gen <- pmin(pmax(as.numeric(1/(1+exp(-X1 %*% sF$bt))),1e-6),1-1e-6)
  etaB <- drop(X1 %*% sF$bt)          # bootstrap world's TRUE index
  gOB  <- deciles(etaB)
  SFd<-SFe<-SOd<-SOe<-numeric(NB)
  for(bb in 1:NB){
    ys <- rbinom(nrow(X1),1,gen)
    f2 <- ridge_fit(X1,ys,lam,D)
    s1 <- stat_g(X1,ys,f2,lam,D,deciles(f2$pi))
    s2 <- stat_g(X1,ys,f2,lam,D,gOB)
    SFd[bb]<-s1$dec; SFe[bb]<-s1$edge; SOd[bb]<-s2$dec; SOe[bb]<-s2$edge
  }
  c(fit_dec =(1+sum(SFd>=sF$dec))/(NB+1),  fit_edge=(1+sum(SFe>=sF$edge,na.rm=TRUE))/(NB+1),
    orc_dec =(1+sum(SOd>=sO$dec))/(NB+1),  orc_edge=(1+sum(SOe>=sO$edge,na.rm=TRUE))/(NB+1))
}

n<-400; p<-100; lam<-416
Sig <- 0.7^abs(outer(1:p,1:p,"-")); Ch <- chol(Sig)
b0  <- rep(c(.35,-.30,.25,-.20,.15), length.out=p)*0.8887
D   <- diag(c(0,rep(1,p)))
vN  <- as.numeric(t(b0) %*% Sig %*% b0)            # var(eta0)

cells <- list(list(k="null" ,g=0),
              list(k="quad" ,g=0.5), list(k="quad" ,g=1.0),
              list(k="cubic",g=0.5), list(k="cubic",g=1.0), list(k="cubic",g=1.5))
out <- NULL; t0 <- Sys.time()
for (cl in cells) {
  P <- matrix(NA,B,4)
  for (b in 1:B) {
    X <- matrix(rnorm(n*p),n,p) %*% Ch; X1 <- cbind(1,X)
    e0 <- drop(X %*% b0)
    eta <- if (cl$k=="cubic") e0 + cl$g*(e0^3 - 3*vN*e0)/sqrt(6*vN^3) else
           if (cl$k=="quad")  e0 + cl$g*(X[,1]^2 - 1) else e0
    y <- rbinom(n,1,1/(1+exp(-eta)))
    P[b,] <- one_rep(X1, y, e0, lam, D)
  }
  rej <- colMeans(P < 0.05)
  out <- rbind(out, data.frame(alt=cl$k, gamma=cl$g,
      fit_dec=rej[1], fit_edge=rej[2], orc_dec=rej[3], orc_edge=rej[4],
      se=round(sqrt(0.25/B),4)))
  cat(sprintf("%-5s g=%.1f | FITTED dec=%.3f edge=%.3f | ORACLE dec=%.3f edge=%.3f | %.0f min\n",
      cl$k, cl$g, rej[1],rej[2],rej[3],rej[4],
      as.numeric(difftime(Sys.time(),t0,units="mins"))))
}
cat("\n===== ORACLE-INDEX EXPERIMENT (panel B, B=",B,", MC SE ~",round(sqrt(.25/B),3),") =====\n",sep="")
print(out, row.names=FALSE)
cat("\nPREDICTION: cubic ORACLE >> cubic FITTED (attenuation); quad BOTH ~0.05 (dilution).\n")
saveRDS(out,"oracle_index_results.rds")
cat("runtime:",round(as.numeric(difftime(Sys.time(),t0,units="mins")),1),"min\n")
