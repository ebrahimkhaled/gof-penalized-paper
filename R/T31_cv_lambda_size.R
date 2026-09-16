# T3.1 -- BJ revision, co-author comment 1: size of the prepivoted corrected test when lambda
# is chosen by 10-fold cross-validation on every simulated dataset instead of being fixed.
# Designs A and B exactly as in T23_rerun.R. Two rules, lambda.1se and lambda.min, on the theory
# scale lambda = n * lambda_glmnet. lambda is chosen once per dataset and then held fixed inside
# the bootstrap, which is what a practitioner does with a cross-validated penalty. NB = 499.
# Output: T31_cv_lambda_pvalues.csv (one row per replicate) and T31_cv_lambda_results.rds.
suppressPackageStartupMessages({library(glmnet); library(parallel)})
G <- 10; NB <- 499
BA <- 2000; BB <- 1000

worker <- function(seed, panel) {
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
  prepiv <- function(X1,y,lam,D,n){
    ft<-ridge_fit(X1,y,lam,D); st<-pieces(X1,y,ft,lam,D)
    gen<-pmin(pmax(as.numeric(1/(1+exp(-X1%*%st$bt))),1e-6),1-1e-6)
    Sd<-Se<-numeric(NB)
    for(bb in 1:NB){ ys<-rbinom(n,1,gen); s2<-pieces(X1,ys,ridge_fit(X1,ys,lam,D),lam,D)
      Sd[bb]<-s2$dec; Se[bb]<-s2$edge }
    c(dec=(1+sum(Sd>=st$dec))/(NB+1), edge=(1+sum(Se>=st$edge,na.rm=TRUE))/(NB+1))
  }
  if (panel=="A") {
    n<-500; p<-5; b0<-c(0.5,-0.4,0.3,-0.3,0.2); D<-diag(c(0,rep(1,p)))
    X<-matrix(rnorm(n*p),n,p); X1<-cbind(1,X); y<-rbinom(n,1,1/(1+exp(-drop(X%*%b0))))
  } else {
    n<-400; p<-100; S<-0.7^abs(outer(1:p,1:p,"-")); Ch<-chol(S)
    b0<-rep(c(.35,-.30,.25,-.20,.15),length.out=p)*0.8887; D<-diag(c(0,rep(1,p)))
    X<-matrix(rnorm(n*p),n,p)%*%Ch; X1<-cbind(1,X); y<-rbinom(n,1,1/(1+exp(-drop(X%*%b0))))
  }
  cv <- cv.glmnet(X, y, family="binomial", alpha=0, nfolds=10)
  l1 <- cv$lambda.1se*n; l0 <- cv$lambda.min*n
  p1 <- prepiv(X1,y,l1,D,n); p0 <- prepiv(X1,y,l0,D,n)
  c(seed=seed, lam_1se=l1, lam_min=l0, dec_1se=p1[["dec"]], edge_1se=p1[["edge"]],
    dec_min=p0[["dec"]], edge_min=p0[["edge"]])
}

nc <- max(1, detectCores() - 2)
cl <- makeCluster(nc); on.exit(stopCluster(cl))
clusterEvalQ(cl, suppressPackageStartupMessages(library(glmnet)))
clusterExport(cl, c("G","NB"))
clusterSetRNGStream(cl, 20260916)          # RNG kind as in the other studies; each replicate is set.seed()
cat("T3.1 cross-validated lambda on", nc, "workers, NB =", NB, "\n")

t0 <- Sys.time(); all <- NULL; summ <- NULL
for (ce in list(list(p="A", seeds=300000L + seq_len(BA)), list(p="B", seeds=310000L + seq_len(BB)))) {
  R <- t(parSapply(cl, ce$seeds, worker, panel=ce$p))
  d <- data.frame(panel=ce$p, R)
  all <- rbind(all, d)
  write.csv(all, "T31_cv_lambda_pvalues.csv", row.names=FALSE)
  for (rule in c("1se","min")) {
    pd <- d[[paste0("dec_",rule)]]; pe <- d[[paste0("edge_",rule)]]
    rd <- mean(pd < 0.05); re <- mean(pe < 0.05, na.rm=TRUE); B <- nrow(d)
    lam <- d[[paste0("lam_",rule)]]
    summ <- rbind(summ, data.frame(panel=ce$p, rule=rule, B=B,
      lam_q25=round(quantile(lam,.25),1), lam_med=round(median(lam),1), lam_q75=round(quantile(lam,.75),1),
      dec=round(rd,4), se_dec=round(sqrt(rd*(1-rd)/B),4),
      edge=round(re,4), se_edge=round(sqrt(re*(1-re)/sum(!is.na(pe))),4),
      edge_na=sum(is.na(pe))))
    cat(sprintf("  %s lambda.%-3s B=%d | median lambda %.1f | dec %.4f (%.4f)  EDGE %.4f (%.4f) | %.0f min\n",
        ce$p, rule, B, median(lam), rd, sqrt(rd*(1-rd)/B), re, sqrt(re*(1-re)/B),
        as.numeric(difftime(Sys.time(),t0,units="mins"))))
  }
  saveRDS(list(pvalues=all, summary=summ), "T31_cv_lambda_results.rds")
}
cat("\n===== T3.1 SIZE WITH A CROSS-VALIDATED PENALTY =====\n"); print(summ, row.names=FALSE)
cat("runtime:", round(as.numeric(difftime(Sys.time(),t0,units="mins")),1), "min\n")
