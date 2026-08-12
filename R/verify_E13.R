# VERIFY E13. Two runs of the SAME cell (design B, lambda=50, EDGE basis, NB=149)
# disagree:
#   T2.3 (8f), B = 1000, seeds 1..1000       -> 0.0270 (SE 0.0051)   "conservative"
#   overshoot_test, 400 reps, seeds 1..400   -> 0.0500 (SE 0.0109)   "nominal"
# Both scripts seed with set.seed(1..N) after clusterSetRNGStream, so if the random
# draws were consumed identically the second would be a SUBSET of the first and could
# not differ this much. Something differs, so neither run can be trusted for this cell
# until an INDEPENDENT high-precision estimate is taken.
#
# This script takes that estimate on a disjoint seed block (100001..101500) and also
# re-runs the ORIGINAL seed block (1..1500) with identical code, so the two are
# comparable under one implementation. If they agree, E13 was a seed artefact of the
# old script; if the original block reproduces 0.027, the effect is real but the seed
# block is unrepresentative and the cell needs far more replication.
suppressPackageStartupMessages({library(parallel)})
G <- 10; NB <- 149; REPS <- 1500

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
# byte-identical to T23_rerun.R's design B block, in the same order of RNG consumption
worker <- function(seed, lam) {
  set.seed(seed)
  n<-400; p<-100; S<-0.7^abs(outer(1:p,1:p,"-")); Ch<-chol(S)
  b0<-rep(c(.35,-.30,.25,-.20,.15),length.out=p)*0.8887; D<-diag(c(0,rep(1,p)))
  X<-matrix(rnorm(n*p),n,p)%*%Ch; X1<-cbind(1,X); y<-rbinom(n,1,1/(1+exp(-drop(X%*%b0))))
  ft<-ridge_fit(X1,y,lam,D); st<-pieces(X1,y,ft,lam,D)
  gen<-pmin(pmax(as.numeric(1/(1+exp(-X1%*%st$bt))),1e-6),1-1e-6)
  Sd<-Se<-numeric(NB)
  for(b in 1:NB){ ys<-rbinom(n,1,gen); s2<-pieces(X1,ys,ridge_fit(X1,ys,lam,D),lam,D)
    Sd[b]<-s2$dec; Se[b]<-s2$edge }
  c(dec=(1+sum(Sd>=st$dec))/(NB+1), edge=(1+sum(Se>=st$edge,na.rm=TRUE))/(NB+1))
}

nc <- max(1, detectCores()-2)
cl <- makeCluster(nc); on.exit(stopCluster(cl))
clusterExport(cl, c("G","NB","ridge_fit","pieces"))
clusterSetRNGStream(cl, 20260807)   # same stream setup as T2.3
cat("verify E13 on", nc, "workers; REPS =", REPS, "per block; NB =", NB, "\n\n")

blk <- function(seeds, lam, tag) {
  t0 <- Sys.time()
  M <- simplify2array(parLapply(cl, seeds, worker, lam=lam))
  rd <- mean(M["dec",]<=0.05); re <- mean(M["edge",]<=0.05); N <- length(seeds)
  cat(sprintf("%-40s | decile %.4f (SE %.4f) | EDGE %.4f (SE %.4f)  (%.0fs)\n",
    tag, rd, sqrt(rd*(1-rd)/N), re, sqrt(re*(1-re)/N),
    as.numeric(difftime(Sys.time(),t0,units="secs"))))
  c(dec=rd, edge=re, n=N)
}

o <- list()
o$orig <- blk(1:REPS,                lam=50, "B lam=50  ORIGINAL seed block 1..1500")
o$ind  <- blk(100001:(100000+REPS),  lam=50, "B lam=50  INDEPENDENT block 100001..")
o$ctrl <- blk(200001:(200000+REPS),  lam=1000,"B lam=1000 control, block 200001..")
saveRDS(o, "verify_E13_results.rds")

p <- prop.test(c(round(o$orig["edge"]*REPS), round(o$ind["edge"]*REPS)), c(REPS,REPS))
cat(sprintf("\nEDGE, original vs independent block: p = %.3f%s\n", p$p.value,
            if (p$p.value > 0.05) "  -> blocks AGREE" else "  -> blocks DIFFER"))
cat("saved verify_E13_results.rds\n")
