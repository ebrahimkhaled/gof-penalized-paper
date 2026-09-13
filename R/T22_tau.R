# =====================================================================
# T2.2 -- POWER RE-PARAMETERISED BY THE VISIBLE EFFECT SIZE tau
#
#   tau = gamma * sd[ E{g(X) | etahat} ] / sd(etahat)
#         estimated by regressing g(X) on poly(etahat,3) and taking
#         sd(fitted)/sd(etahat)                      (brief's definition)
#
# PREDICTION UNDER TEST: the two panels' power curves approximately collapse
# onto ONE curve when plotted against tau.
#
# COST DECISION (budget = shared machine).  The prepivoted bootstrap test costs
# ~3.5 s per test in panel B, so a 20-cell power surface at B=1000 would need
# ~19 h.  Instead every arm here is calibrated by DIRECT known-null Monte Carlo
# exactly as in `oracle_knownnull.R` (manifest 7.5b): the null distribution of
# the statistic is built from the TRUE model (B0 draws), its 95th percentile is
# the critical value, and rejection is measured under each alternative (B1
# draws).  Every cell is then EXACTLY 5% under H0 by construction, which is what
# a like-for-like curve comparison needs, and it is ~50x cheaper.  The manifest
# records that this construction and the prepivoted bootstrap agree on the
# fitted arm in panel B (0.040-0.092 vs 0.047-0.085).
# =====================================================================
set.seed(2026)
G <- 10; B0 <- 2000; B1 <- 1000
t_start <- Sys.time()

## ---- verbatim from oracle_knownnull.R / edge_basis_study.R ----------
ridge_fit <- function(X1,y,lam,D,tol=1e-9,maxit=60){
  beta <- rep(0,ncol(X1))
  for(it in 1:maxit){ eta<-drop(X1%*%beta); pi<-1/(1+exp(-eta))
    w<-pmax(pi*(1-pi),1e-8); z<-eta+(y-pi)/w
    bn<-drop(solve(crossprod(X1,w*X1)+lam*D, crossprod(X1,w*z)))
    if(max(abs(bn-beta))<tol){beta<-bn;break}; beta<-bn }
  list(beta=beta, pi=1/(1+exp(-drop(X1%*%beta))))
}
stat_g <- function(X1,y,fit,lam,D,g){
  pi<-pmin(pmax(fit$pi,1e-6),1-1e-6); w<-pmax(pi*(1-pi),1e-8)
  idx<-split(seq_along(y),g); Vg<-vapply(idx,function(I) sum(w[I]),0.0)
  r<-(vapply(idx,function(I) sum(y[I]),0.0)-vapply(idx,function(I) sum(pi[I]),0.0))/sqrt(Vg)
  U<-t(vapply(idx,function(I) colSums(w[I]*X1[I,,drop=FALSE]),numeric(ncol(X1))))/sqrt(Vg)
  Fm<-crossprod(X1,w*X1); K<-lam*D
  bt<-fit$beta+drop(solve(Fm,K%*%fit$beta)); mu<-drop(U%*%solve(Fm+K,K%*%bt))
  pbar<-vapply(idx,function(I) mean(pi[I]),0.0); v<-r-mu
  Z<-tryCatch(as.matrix(stats::poly(pbar,3)),error=function(e) NULL)
  edge<-if(is.null(Z)) NA_real_ else {zr<-crossprod(Z,v); as.numeric(t(zr)%*%solve(crossprod(Z))%*%zr)}
  c(dec=sum(v^2), edge=edge)
}
dec10 <- function(x) pmin(ceiling(rank(x,ties.method="first")/(length(x)/G)),G)

## ---- tau estimator --------------------------------------------------
tau_hat <- function(gx, etahat, gam) {
  n <- length(gx); Zp <- cbind(1, stats::poly(etahat,3))
  ft <- lm.fit(Zp, gx); fv <- ft$fitted.values
  ssf <- sum((fv-mean(fv))^2); s2 <- sum(ft$residuals^2)/(n-4)
  sdE <- sd(etahat)
  c(raw = gam*sqrt(ssf/(n-1))/sdE,                        # brief's plug-in
    dfc = gam*sqrt(max(0, ssf-3*s2)/(n-1))/sdE)           # df-corrected
}

## ---- designs (frozen) ----------------------------------------------
nA<-500; pA<-5; bA<-c(0.5,-0.4,0.3,-0.3,0.2); DA<-diag(c(0,rep(1,pA))); vA<-sum(bA^2)
nB<-400; pB<-100; SigB<-0.7^abs(outer(1:pB,1:pB,"-")); ChB<-chol(SigB)
bB<-rep(c(.35,-.30,.25,-.20,.15),length.out=pB)*0.8887; DB<-diag(c(0,rep(1,pB)))
vB<-as.numeric(t(bB)%*%SigB%*%bB)

mk <- function(panel,g,k){
  if(panel=="A"){ X<-matrix(rnorm(nA*pA),nA,pA); e0<-drop(X%*%bA); v<-vA }
  else          { X<-matrix(rnorm(nB*pB),nB,pB)%*%ChB; e0<-drop(X%*%bB); v<-vB }
  gx <- if(k=="cubic") (e0^3-3*v*e0)/sqrt(6*v^3) else if(k=="quad") X[,1]^2-1 else rep(0,length(e0))
  list(X1=cbind(1,X), y=rbinom(length(e0),1,1/(1+exp(-(e0+g*gx)))), gx=gx, e0=e0)
}
draw <- function(panel,g,k,lam,D){
  d <- mk(panel,g,k); ft <- ridge_fit(d$X1,d$y,lam,D)
  s <- stat_g(d$X1,d$y,ft,lam,D,dec10(ft$pi))
  eh <- drop(d$X1%*%ft$beta); th <- tau_hat(d$gx, eh, g)
  c(dec=unname(s[1]), edge=unname(s[2]), tau=unname(th[1]), tau_dfc=unname(th[2]),
    sd_eta=sd(eh), rho=cor(eh,d$e0))
}
bank <- function(panel,g,k,lam,D,B) t(replicate(B, draw(panel,g,k,lam,D)))

## ---- gamma grids chosen from pilot c = tau/gamma to hit a COMMON tau grid ----
TAU_TARGET <- c(0.4, 0.8, 1.2, 1.8, 2.5)
cpil <- list("A.quad"=1.176, "A.cubic"=3.245, "B.quad"=0.358, "B.cubic"=0.762)
gam_grid <- lapply(cpil, function(cc) round(pmin(TAU_TARGET/cc, 6), 3))
cat("gamma grids chosen (capped at 6):\n"); print(gam_grid)

res <- NULL
for (panel in c("A","B")) {
  lam <- if(panel=="A") 100 else 416; D <- if(panel=="A") DA else DB
  cat(sprintf("\n### PANEL %s (lambda=%d) null bank B0=%d\n", panel, lam, B0))
  N <- bank(panel,0,"null",lam,D,B0)
  cv <- apply(N[,1:2],2,quantile,0.95,na.rm=TRUE)
  cat(sprintf("  cv: dec=%.2f edge=%.2f | sd(etahat)=%.3f rho=%.3f | %.1f min\n",
      cv[1],cv[2],mean(N[,"sd_eta"]),mean(N[,"rho"]),
      as.numeric(difftime(Sys.time(),t_start,units="mins"))))
  res <- rbind(res, data.frame(panel=panel,alt="null",gamma=0,
        tau=0,tau_se=0,tau_dfc=0,pow_dec=0.05,pow_edge=0.05,se=0,B=B0,
        sd_eta=mean(N[,"sd_eta"]),rho=mean(N[,"rho"])))
  for (k in c("quad","cubic")) for (g in gam_grid[[paste0(panel,".",k)]]) {
    A <- bank(panel,g,k,lam,D,B1)
    pd <- mean(A[,1]>cv[1]); pe <- mean(A[,2]>cv[2],na.rm=TRUE)
    res <- rbind(res, data.frame(panel=panel,alt=k,gamma=g,
      tau=mean(A[,"tau"]), tau_se=sd(A[,"tau"])/sqrt(B1), tau_dfc=mean(A[,"tau_dfc"]),
      pow_dec=pd, pow_edge=pe, se=round(sqrt(0.25/B1),4), B=B1,
      sd_eta=mean(A[,"sd_eta"]), rho=mean(A[,"rho"])))
    cat(sprintf("  %-5s g=%.3f | tau=%.3f (dfc %.3f) | dec=%.3f edge=%.3f | %.1f min\n",
        k,g,mean(A[,"tau"]),mean(A[,"tau_dfc"]),pd,pe,
        as.numeric(difftime(Sys.time(),t_start,units="mins"))))
    saveRDS(res,"T22_tau_results.rds")
  }
}
cat("\n===== T2.2  POWER vs VISIBLE EFFECT SIZE tau  (known-null calibrated, both arms exactly 5%) =====\n")
print(res[,c("panel","alt","gamma","tau","tau_dfc","pow_dec","pow_edge","se","B")], row.names=FALSE)
saveRDS(res,"T22_tau_results.rds")
cat(sprintf("\ntotal %.1f min\n", as.numeric(difftime(Sys.time(),t_start,units="mins"))))
