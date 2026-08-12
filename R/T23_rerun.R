# T2.3 RERUN — headline SIZE table at high precision.
# The first attempt overran its budget and produced nothing citable; §5's B=300
# numbers are still the authority. This run is PARALLELISED (PSOCK, Windows) so it
# actually finishes: measured cost is 0.26 s/test in DGP A and 3.3 s/test in DGP B
# at NB=149, so serially DGP B at B=1000 would be ~55 min per cell.
#
# Targets: DGP A  lambda in {100,137,200}  at B = 2000  (MC SE 0.0049 at p=0.05)
#          DGP B  lambda in {50,416,1000}  at B = 1000  (MC SE 0.0069)
# Both bases (decile + EDGE-poly3), NB = 149, prepivoted corrected test, H0 TRUE.
suppressPackageStartupMessages({library(glmnet); library(parallel)})
G <- 10; NB <- 149
BA <- 2000; BB <- 1000

worker <- function(seed, panel, lam) {
  set.seed(seed)
  ridge_fit <- function(X1,y,lam,D,tol=1e-9,maxit=60){
    beta <- rep(0,ncol(X1))
    for(it in 1:maxit){ eta<-drop(X1%*%beta); pi<-1/(1+exp(-eta))
      w<-pmax(pi*(1-pi),1e-8); z<-eta+(y-pi)/w
      bn<-drop(solve(crossprod(X1,w*X1)+lam*D, crossprod(X1,w*z)))
      if(max(abs(bn-beta))<tol){beta<-bn;break}; beta<-bn }
    list(beta=beta, pi=1/(1+exp(-drop(X1%*%beta))))
  }
  pieces <- function(X1,y,fit,lam,D){
    pi<-pmin(pmax(fit$pi,1e-6),1-1e-6); w<-pmax(pi*(1-pi),1e-8)
    g<-pmin(ceiling(rank(pi,ties.method="first")/(length(y)/G)),G)
    idx<-split(seq_along(y),g); Vg<-vapply(idx,function(I) sum(w[I]),0.0)
    r<-(vapply(idx,function(I) sum(y[I]),0.0)-vapply(idx,function(I) sum(pi[I]),0.0))/sqrt(Vg)
    U<-t(vapply(idx,function(I) colSums(w[I]*X1[I,,drop=FALSE]),numeric(ncol(X1))))/sqrt(Vg)
    Fm<-crossprod(X1,w*X1); K<-lam*D
    bt<-fit$beta+drop(solve(Fm,K%*%fit$beta)); mu<-drop(U%*%solve(Fm+K,K%*%bt))
    pbar<-vapply(idx,function(I) mean(pi[I]),0.0); v<-r-mu
    Z<-tryCatch(as.matrix(stats::poly(pbar,3)),error=function(e) NULL)
    ed<-if(is.null(Z)) NA_real_ else {zr<-crossprod(Z,v); as.numeric(t(zr)%*%solve(crossprod(Z))%*%zr)}
    list(dec=sum(v^2), edge=ed, bt=bt)
  }
  if (panel=="A") {
    n<-500; p<-5; b0<-c(0.5,-0.4,0.3,-0.3,0.2); D<-diag(c(0,rep(1,p)))
    X<-matrix(rnorm(n*p),n,p); X1<-cbind(1,X); y<-rbinom(n,1,1/(1+exp(-drop(X%*%b0))))
  } else {
    n<-400; p<-100; S<-0.7^abs(outer(1:p,1:p,"-")); Ch<-chol(S)
    b0<-rep(c(.35,-.30,.25,-.20,.15),length.out=p)*0.8887; D<-diag(c(0,rep(1,p)))
    X<-matrix(rnorm(n*p),n,p)%*%Ch; X1<-cbind(1,X); y<-rbinom(n,1,1/(1+exp(-drop(X%*%b0))))
  }
  ft<-ridge_fit(X1,y,lam,D); st<-pieces(X1,y,ft,lam,D)
  gen<-pmin(pmax(as.numeric(1/(1+exp(-X1%*%st$bt))),1e-6),1-1e-6)
  Sd<-Se<-numeric(NB)
  for(bb in 1:NB){ ys<-rbinom(n,1,gen); s2<-pieces(X1,ys,ridge_fit(X1,ys,lam,D),lam,D)
    Sd[bb]<-s2$dec; Se[bb]<-s2$edge }
  c(dec=(1+sum(Sd>=st$dec))/(NB+1), edge=(1+sum(Se>=st$edge,na.rm=TRUE))/(NB+1))
}

nc <- max(1, detectCores() - 2)
cl <- makeCluster(nc); on.exit(stopCluster(cl))
clusterEvalQ(cl, suppressPackageStartupMessages(library(glmnet)))
clusterExport(cl, c("G","NB"))
clusterSetRNGStream(cl, 20260807)
cat("T2.3 rerun on", nc, "workers\n")

cells <- list(list(p="A",lam=100,B=BA), list(p="A",lam=137,B=BA), list(p="A",lam=200,B=BA),
              list(p="B",lam=50, B=BB), list(p="B",lam=416,B=BB), list(p="B",lam=1000,B=BB))
out <- NULL; t0 <- Sys.time()
for (cl_ in cells) {
  seeds <- seq_len(cl_$B) + 1000L
  R <- parSapply(cl, seeds, worker, panel=cl_$p, lam=cl_$lam)
  rd <- mean(R["dec",] < 0.05); re <- mean(R["edge",] < 0.05, na.rm=TRUE)
  out <- rbind(out, data.frame(panel=cl_$p, lambda=cl_$lam, B=cl_$B,
    dec=round(rd,4), se_dec=round(sqrt(rd*(1-rd)/cl_$B),4),
    edge=round(re,4), se_edge=round(sqrt(re*(1-re)/cl_$B),4),
    ok = rd>=0.03 & rd<=0.08 & re>=0.03 & re<=0.08))
  cat(sprintf("  %s lam=%-5g B=%d | dec %.4f (%.4f)  EDGE %.4f (%.4f)  %s | %.0f min\n",
      cl_$p, cl_$lam, cl_$B, rd, sqrt(rd*(1-rd)/cl_$B), re, sqrt(re*(1-re)/cl_$B),
      ifelse(rd>=0.03&rd<=0.08&re>=0.03&re<=0.08,"PASS","** OUTSIDE [0.03,0.08] **"),
      as.numeric(difftime(Sys.time(),t0,units="mins"))))
  saveRDS(out,"T23_rerun_results.rds")
}
cat("\n===== T2.3 HEADLINE SIZE, HIGH PRECISION =====\n"); print(out,row.names=FALSE)
cat("\nCompare to the B=300 values in manifest section 5:\n")
cat("  A: 0.060 / 0.043 / 0.027   B: 0.070 / 0.060 / 0.043 (decile)\n")
cat("runtime:",round(as.numeric(difftime(Sys.time(),t0,units="mins")),1),"min\n")
