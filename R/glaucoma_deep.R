# DEEP ROBUSTNESS CHECK on the primary Series C application: GlaucomaM.
# The single-shot result was NAIVE p=0.0000 -> CORRECTED p=0.085 (conclusion flips).
# Before that goes in a paper it must survive: choice of lambda, choice of G,
# more bootstrap replicates, and the seed.
suppressPackageStartupMessages({library(glmnet); library(TH.data)})
data("GlaucomaM", package="TH.data")
Xr <- as.matrix(GlaucomaM[, setdiff(names(GlaucomaM),"Class")])
y  <- as.numeric(GlaucomaM$Class == "glaucoma")
X  <- scale(Xr); X <- X[, apply(X,2,function(c) all(is.finite(c)) & sd(c)>0), drop=FALSE]
n <- nrow(X); p <- ncol(X); X1 <- cbind(1,X); D <- diag(c(0,rep(1,p)))
cat(sprintf("GlaucomaM: n=%d p=%d kappa=%.3f events=%d cond=%.3g max|cor|=%.4f\n",
    n,p,p/n,sum(y),kappa(crossprod(X)),max(abs(cor(X)[upper.tri(cor(X))]))))
cat("MLE: "); print(tryCatch({glm(y~X,family=binomial); "converged"},
   warning=function(w) paste("WARNING:",conditionMessage(w)), error=function(e) "ERROR"))

ridge_fit <- function(X1,y,lam,D,tol=1e-9,maxit=200){
  beta <- rep(0,ncol(X1))
  for(it in 1:maxit){ eta<-drop(X1%*%beta); pi<-1/(1+exp(-eta))
    w<-pmax(pi*(1-pi),1e-8); z<-eta+(y-pi)/w
    bn<-drop(solve(crossprod(X1,w*X1)+lam*D, crossprod(X1,w*z)))
    if(max(abs(bn-beta))<tol){beta<-bn;break}; beta<-bn }
  list(beta=beta, pi=1/(1+exp(-drop(X1%*%beta))))
}
pieces <- function(X1,y,fit,lam,D,G){
  pi<-pmin(pmax(fit$pi,1e-6),1-1e-6); w<-pmax(pi*(1-pi),1e-8)
  g<-pmin(ceiling(rank(pi,ties.method="first")/(length(y)/G)),G)
  idx<-split(seq_along(y),g); Vg<-vapply(idx,function(I) sum(w[I]),0.0)
  r<-(vapply(idx,function(I) sum(y[I]),0.0)-vapply(idx,function(I) sum(pi[I]),0.0))/sqrt(Vg)
  U<-t(vapply(idx,function(I) colSums(w[I]*X1[I,,drop=FALSE]),numeric(ncol(X1))))/sqrt(Vg)
  Fm<-crossprod(X1,w*X1); K<-lam*D
  bt<-fit$beta+drop(solve(Fm,K%*%fit$beta)); mu<-drop(U%*%solve(Fm+K,K%*%bt))
  pbar<-vapply(idx,function(I) mean(pi[I]),0.0); v<-r-mu
  Z<-tryCatch(as.matrix(stats::poly(pbar,3)),error=function(e) NULL)
  qf<-function(u,Z){zr<-crossprod(Z,u); as.numeric(t(zr)%*%solve(crossprod(Z))%*%zr)}
  list(Sp=sum(r^2), Sc=sum(v^2),
       Spe=if(is.null(Z)) NA else qf(r,Z), Sce=if(is.null(Z)) NA else qf(v,Z),
       Om=diag(G)-U%*%solve(Fm,t(U)), Z=Z, bt=bt, G=G)
}
mc_p <- function(S,Om,Z=NULL,nd=20000){
  R<-chol(Om+diag(1e-9,nrow(Om))); W<-matrix(rnorm(nd*nrow(Om)),nd)%*%R
  if(is.null(Z)) return(mean(rowSums(W^2)>=S))
  A<-solve(crossprod(Z)); mean(apply(W,1,function(u){zr<-crossprod(Z,u); as.numeric(t(zr)%*%A%*%zr)})>=S)
}
run <- function(lam,G,NB,seed){
  set.seed(seed)
  ft<-ridge_fit(X1,y,lam,D); st<-pieces(X1,y,ft,lam,D,G)
  pn_d<-mc_p(st$Sp,st$Om); pn_e<-if(is.null(st$Z)) NA else mc_p(st$Spe,st$Om,st$Z)
  gen<-pmin(pmax(as.numeric(1/(1+exp(-X1%*%st$bt))),1e-6),1-1e-6)
  Sd<-Se<-numeric(NB)
  for(b in 1:NB){ ys<-rbinom(n,1,gen); s2<-pieces(X1,ys,ridge_fit(X1,ys,lam,D),lam,D,G)
    Sd[b]<-s2$Sc; Se[b]<-s2$Sce }
  c(naive_dec=pn_d, corr_dec=(1+sum(Sd>=st$Sc))/(NB+1),
    naive_edge=pn_e, corr_edge=(1+sum(Se>=st$Sce,na.rm=TRUE))/(NB+1))
}
set.seed(1); cv <- cv.glmnet(X,y,family="binomial",alpha=0)
lmin <- cv$lambda.min*n; l1se <- cv$lambda.1se*n
cat(sprintf("\ncv: lambda.min=%.4g  lambda.1se=%.4g  (theory scale: %.4g / %.4g)\n",
    cv$lambda.min, cv$lambda.1se, lmin, l1se))

cat("\n=== ROBUSTNESS GRID (NB=499) ===\n")
res <- NULL
for (lab in c("lambda.min","lambda.1se","0.5*1se","2*1se")) {
  lam <- switch(lab, "lambda.min"=lmin, "lambda.1se"=l1se, "0.5*1se"=l1se/2, "2*1se"=2*l1se)
  for (G in c(5,10,20)) {
    r <- run(lam,G,499,seed=100+G)
    res <- rbind(res, data.frame(lambda=lab, lam_theory=round(lam,1), G=G,
      naive_dec=round(r["naive_dec"],4), corr_dec=round(r["corr_dec"],4),
      naive_edge=round(r["naive_edge"],4), corr_edge=round(r["corr_edge"],4)))
    cat(sprintf("  %-10s G=%2d | naive dec=%.4f edge=%.4f | corrected dec=%.4f edge=%.4f\n",
        lab,G,r["naive_dec"],r["naive_edge"],r["corr_dec"],r["corr_edge"]))
  }
}
cat("\n=== SEED STABILITY at lambda.1se, G=10, NB=999 (5 seeds) ===\n")
for (s in 1:5) { r<-run(l1se,10,999,seed=1000+s)
  cat(sprintf("  seed %d | corrected dec=%.4f edge=%.4f\n", s, r["corr_dec"], r["corr_edge"])) }
print(res,row.names=FALSE); saveRDS(res,"glaucoma_deep_results.rds")
