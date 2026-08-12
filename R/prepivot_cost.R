# DIAGNOSTIC: does PREPIVOTING itself cost power under H1?
# Hypothesis: the bootstrap generator pi(beta.tilde) is fitted to the OBSERVED data,
# so under H1 it has already absorbed part of the alternative; the bootstrap world is
# therefore "closer to" the data than a true null world, S_c* is inflated, and power falls.
# Compare, on the SAME data, in panel A (where BOTH references are valid):
#    ANALYTIC : S_c vs a Monte Carlo draw from N(0, Omega_MLE)     (Remark A.1)
#    PREPIVOT : S_c vs the bootstrap distribution of S_c*          (the proposed test)
suppressPackageStartupMessages(library(glmnet)); set.seed(8801)
G <- 10; NB <- 149; ND <- 4000; B <- 300

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
  bt<-fit$beta+drop(solve(Fm,K%*%fit$beta))
  mu<-drop(U%*%solve(Fm+K,K%*%bt))
  list(Sc=sum((r-mu)^2), bt=bt, Om=diag(G)-U%*%solve(Fm,t(U)))
}
n<-500;p<-5;bA<-c(0.5,-0.4,0.3,-0.3,0.2);D<-diag(c(0,rep(1,p))); lam<-100
res<-NULL; t0<-Sys.time()
for(gam in c(0,0.5,1.0)){
  hA<-hP<-numeric(B)
  for(b in 1:B){
    X<-matrix(rnorm(n*p),n,p); X1<-cbind(1,X)
    eta<-drop(X%*%bA)+gam*(X[,1]^2-1)
    y<-rbinom(n,1,1/(1+exp(-eta)))
    ft<-ridge_fit(X1,y,lam,D); st<-pieces(X1,y,ft,lam,D)
    # ANALYTIC reference
    R<-chol(st$Om+diag(1e-9,G)); Z<-matrix(rnorm(ND*G),ND)%*%R
    hA[b] <- mean(rowSums(Z^2)>=st$Sc) < 0.05
    # PREPIVOT reference
    gen<-pmin(pmax(as.numeric(1/(1+exp(-X1%*%st$bt))),1e-6),1-1e-6)
    Ss<-numeric(NB)
    for(bb in 1:NB){ ys<-rbinom(n,1,gen)
      Ss[bb]<-pieces(X1,ys,ridge_fit(X1,ys,lam,D),lam,D)$Sc }
    hP[b] <- (1+sum(Ss>=st$Sc))/(NB+1) < 0.05
  }
  res<-rbind(res,data.frame(gamma=gam, analytic=mean(hA), prepivot=mean(hP),
                            se=round(sqrt(0.25/B),4)))
  cat(sprintf("gamma=%.1f | analytic=%.3f  prepivot=%.3f  (%.0f min)\n",
      gam, mean(hA), mean(hP), as.numeric(difftime(Sys.time(),t0,units="mins"))))
}
cat("\n===== COST OF PREPIVOTING (panel A, lambda=100, B=",B,") =====\n",sep="")
print(res,row.names=FALSE)
cat("\ngamma=0 row is SIZE (both must be ~0.05). gamma>0 rows: any gap = the power price of prepivoting.\n")
saveRDS(res,"prepivot_cost_results.rds")
