# T39 -- the GlaucomaMVF grid of Table S11 without lora.
# The class of GlaucomaMVF is reproduced in all 170 eyes by the rule in the ipred documentation,
# which uses clv, cs and lora; lora is the one of those three that T26_app.R kept as a predictor.
# Identical to part (b) of T26_app.R (same functions, seeds, penalties by cross-validation,
# NB = 499 grid and the five-seed NB = 999 check) with lora removed, so p = 62.
suppressPackageStartupMessages({library(glmnet); library(TH.data); library(ipred)})

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
mk_run <- function(X1,y,D,n) function(lam,G,NB,seed){
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

## ============================ (b) GlaucomaMVF ==============================
data("GlaucomaMVF", package="ipred")
dfv <- GlaucomaMVF
Xr <- as.matrix(dfv[, setdiff(names(dfv), c("Class", "lora"))])   # lora enters the class rule
yv <- as.numeric(dfv$Class == "glaucoma")
Xv <- scale(Xr); Xv <- Xv[, apply(Xv,2,function(c) all(is.finite(c)) & sd(c)>0), drop=FALSE]
nv <- nrow(Xv); pv <- ncol(Xv); X1v <- cbind(1,Xv); Dv <- diag(c(0,rep(1,pv)))
cat(sprintf("GlaucomaMVF: n=%d p=%d kappa=%.3f events=%d cond=%.3g max|cor|=%.4f\n",
    nv,pv,pv/nv,sum(yv),kappa(crossprod(Xv)),max(abs(cor(Xv)[upper.tri(cor(Xv))]))))
cat("MLE: "); print(tryCatch({fm<-glm(yv~Xv,family=binomial)
    if(max(abs(coef(fm)),na.rm=TRUE)>50) "SEPARATION (|beta| huge)" else "converged"},
   warning=function(w) paste("WARNING:",conditionMessage(w)), error=function(e) "ERROR"))

set.seed(1); cvv <- cv.glmnet(Xv,yv,family="binomial",alpha=0)
lminV <- cvv$lambda.min*nv; l1seV <- cvv$lambda.1se*nv
cat(sprintf("cv: lambda.min=%.4g lambda.1se=%.4g (theory scale %.4g / %.4g)\n",
    cvv$lambda.min, cvv$lambda.1se, lminV, l1seV))
runV <- mk_run(X1v,yv,Dv,nv)

cat("\n=== (b) GlaucomaMVF ROBUSTNESS GRID (NB=499) ===\n")
resV <- NULL
for (lab in c("lambda.min","0.5*1se","lambda.1se","2*1se")) {
  lam <- switch(lab,"lambda.min"=lminV,"0.5*1se"=l1seV/2,"lambda.1se"=l1seV,"2*1se"=2*l1seV)
  for (G in c(5,10,20)) {
    r <- runV(lam,G,499,seed=100+G)
    resV <- rbind(resV, data.frame(dataset="GlaucomaMVF", lambda=lab,
      lam_theory=round(lam,1), G=G,
      naive_dec=round(r["naive_dec"],4), corr_dec=round(r["corr_dec"],4),
      naive_edge=round(r["naive_edge"],4), corr_edge=round(r["corr_edge"],4)))
    cat(sprintf("  %-10s G=%2d | naive dec=%.4f edge=%.4f | corrected dec=%.4f edge=%.4f\n",
        lab,G,r["naive_dec"],r["naive_edge"],r["corr_dec"],r["corr_edge"]))
  }
}
cat("\n=== (b) SEED STABILITY at lambda.1se, G=10, NB=999 (5 seeds) ===\n")
stab <- NULL
for (s in 1:5) { r<-runV(l1seV,10,999,seed=2000+s)
  stab <- rbind(stab, data.frame(seed=s, corr_dec=r["corr_dec"], corr_edge=r["corr_edge"]))
  cat(sprintf("  seed %d | corrected dec=%.4f edge=%.4f\n", s, r["corr_dec"], r["corr_edge"])) }
print(resV, row.names=FALSE)
# AUC + fitted-probability spread, to guard against the wpbc artefact (E7)
for (lab in c("lambda.min","lambda.1se","2*1se")) {
  lam <- switch(lab,"lambda.min"=lminV,"lambda.1se"=l1seV,"2*1se"=2*l1seV)
  ft <- ridge_fit(X1v,yv,lam,Dv); ph <- ft$pi
  au <- mean(outer(ph[yv==1], ph[yv==0], ">") + 0.5*outer(ph[yv==1], ph[yv==0], "=="))
  cat(sprintf("  %-10s lam=%8.2f | AUC=%.3f sd(pi.hat)=%.3f range=[%.3f,%.3f] max|beta|=%.4f\n",
      lab, lam, au, sd(ph), min(ph), max(ph), max(abs(ft$beta[-1]))))
}
saveRDS(list(grid=resV, stab=stab), file.path("..", "data", "T39_mvf_nolora_results.rds"))
write.csv(resV, file.path("..", "data", "T39_mvf_nolora_grid.csv"), row.names = FALSE)
write.csv(stab, file.path("..", "data", "T39_mvf_nolora_stability.csv"), row.names = FALSE)
stopifnot(ncol(Xv) == 62, !("lora" %in% colnames(Xv)))

