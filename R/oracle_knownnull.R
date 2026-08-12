# T2.1c — ORACLE-INDEX COMPARISON DONE CORRECTLY (known-null Monte Carlo).
#
# WHY THE BOOTSTRAP VERSIONS FAILED. An "oracle" test cannot be built on a
# prepivoted bootstrap reference, because the bootstrap world's truth is
# beta.tilde, not beta0, so the two worlds can never be grouped symmetrically:
#   * group the bootstrap by X %*% beta.tilde  -> its oracle is NOISY while the
#     real world's oracle is EXACT  -> reference inflated -> conservative
#     (measured: cubic gamma=0.5 gave oracle EDGE = 0.005, below nominal)
#   * group the bootstrap by eta0 as well      -> eta0 is NOT the bootstrap's
#     truth, so the bootstrap statistic is attenuated while the real one is not
#     -> reference deflated -> ANTI-conservative
#     (measured: size 0.100 / 0.150 at gamma = 0)
# Neither is a valid test. The oracle is not a procedure a practitioner could run;
# it is a BENCHMARK. So benchmark it properly:
#
# THE CORRECT DESIGN. For each grouping rule separately, build the null
# distribution of the statistic by direct Monte Carlo from the TRUE model
# (knowledge no practitioner has - that is the point), take its 95th percentile
# as the critical value, then measure the rejection rate under each alternative.
# Both arms are then EXACTLY calibrated by construction, so any difference in
# power is attributable to the grouping alone - which is precisely the
# rho^{2k} attenuation we want to isolate. No bootstrap, no reference confound,
# and ~50x cheaper.
suppressPackageStartupMessages(library(glmnet))
set.seed(6060)
G <- 10; B0 <- 2000; B1 <- 1000

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

# one dataset -> both statistics under BOTH groupings
draw <- function(gen, lam, D) {
  d <- gen(); X1 <- d$X1; y <- d$y
  ft <- ridge_fit(X1,y,lam,D)
  sF <- stat_g(X1,y,ft,lam,D,dec10(ft$pi))
  sO <- stat_g(X1,y,ft,lam,D,dec10(d$eta0))
  c(fit_dec=sF["dec"], fit_edge=sF["edge"], orc_dec=sO["dec"], orc_edge=sO["edge"])
}
bank <- function(gen, lam, D, B) t(replicate(B, draw(gen,lam,D)))

run_panel <- function(nm, mk, lam, D, alts) {
  cat("\n### ", nm, " (lambda=",lam,", B0=",B0,", B1=",B1,")\n",sep="")
  t0 <- Sys.time()
  N <- bank(function() mk(0,"null"), lam, D, B0)
  cv <- apply(N, 2, quantile, 0.95, na.rm=TRUE)      # exact critical values
  cat(sprintf("  null bank done (%.1f min). 95%% critical values: fit_dec=%.2f fit_edge=%.2f orc_dec=%.2f orc_edge=%.2f\n",
      as.numeric(difftime(Sys.time(),t0,units="mins")), cv[1],cv[2],cv[3],cv[4]))
  out <- NULL
  for (a in alts) {
    A <- bank(function() mk(a$g, a$k), lam, D, B1)
    rej <- colMeans(sweep(A, 2, cv, ">"), na.rm=TRUE)
    out <- rbind(out, data.frame(panel=nm, alt=a$k, gamma=a$g,
      fit_dec=round(rej[1],3), fit_edge=round(rej[2],3),
      orc_dec=round(rej[3],3), orc_edge=round(rej[4],3),
      ratio_edge=round(rej[4]/max(rej[2],1e-9),2)))
    cat(sprintf("  %-5s g=%.1f | FITTED dec=%.3f edge=%.3f | ORACLE dec=%.3f edge=%.3f | EDGE ratio %.1fx\n",
        a$k,a$g,rej[1],rej[2],rej[3],rej[4],rej[4]/max(rej[2],1e-9)))
  }
  out
}

## PANEL B: n=400, p=100, AR(1) 0.7, lambda=416  (kappa=0.25)
nB<-400; pB<-100; SigB<-0.7^abs(outer(1:pB,1:pB,"-")); ChB<-chol(SigB)
bB<-rep(c(.35,-.30,.25,-.20,.15),length.out=pB)*0.8887; DB<-diag(c(0,rep(1,pB)))
vB<-as.numeric(t(bB)%*%SigB%*%bB)
mkB <- function(g,k){ X<-matrix(rnorm(nB*pB),nB,pB)%*%ChB; e0<-drop(X%*%bB)
  eta <- if(k=="cubic") e0+g*(e0^3-3*vB*e0)/sqrt(6*vB^3) else if(k=="quad") e0+g*(X[,1]^2-1) else e0
  list(X1=cbind(1,X), y=rbinom(nB,1,1/(1+exp(-eta))), eta0=e0) }

## PANEL A: n=500, p=5, lambda=100 (kappa=0.01) -- rho ~ 0.97, expect oracle ~ fitted
nA<-500; pA<-5; bA<-c(0.5,-0.4,0.3,-0.3,0.2); DA<-diag(c(0,rep(1,pA))); vA<-sum(bA^2)
mkA <- function(g,k){ X<-matrix(rnorm(nA*pA),nA,pA); e0<-drop(X%*%bA)
  eta <- if(k=="cubic") e0+g*(e0^3-3*vA*e0)/sqrt(6*vA^3) else if(k=="quad") e0+g*(X[,1]^2-1) else e0
  list(X1=cbind(1,X), y=rbinom(nA,1,1/(1+exp(-eta))), eta0=e0) }

alts <- list(list(k="cubic",g=0.5), list(k="cubic",g=1.0), list(k="cubic",g=1.5),
             list(k="quad", g=0.5), list(k="quad", g=1.0))
res <- rbind(run_panel("B (kappa=0.25)", mkB, 416, DB, alts),
             run_panel("A (kappa=0.01)", mkA, 100, DA, alts))
cat("\n===== T2.1c ORACLE vs FITTED GROUPING, exactly calibrated both arms =====\n")
print(res, row.names=FALSE)
cat("\nPREDICTION: panel B cubic -> ORACLE >> FITTED (rho^{2k} attenuation, rho=0.723);\n")
cat("            panel B quad  -> ORACLE ~ FITTED (dilution: signal not in the index at all);\n")
cat("            panel A both  -> ORACLE ~ FITTED (rho=0.969, almost nothing to recover).\n")
saveRDS(res,"oracle_knownnull_results.rds")
