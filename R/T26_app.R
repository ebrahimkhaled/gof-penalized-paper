# T2.6(b)+(c) -- APPLICATION COMPLETION
#
# (b) Full lambda x G grid (12 configurations, NB=499, naive + corrected, both
#     bases) for the SECONDARY real dataset.  Chosen: ipred::GlaucomaMVF --
#     a second, independent glaucoma dataset (n=170, kappa=0.37, MLE separates),
#     so it speaks to the same clinical question as the primary example
#     TH.data::GlaucomaM.  Manifest 8b currently records it as a single
#     configuration (naive 0.0034/0.0000 -> corrected 0.670/0.710), AMBER.
#
# (c) glaucoma_calibration.png -- observed event rate against mean predicted
#     probability by decile for the GlaucomaM ridge fit at lambda.1se.
#
# ridge_fit / pieces / mc_p / run are reused VERBATIM from glaucoma_deep.R.
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
Xr <- as.matrix(dfv[, setdiff(names(dfv), "Class")])
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
saveRDS(list(grid=resV, stab=stab), "T26_app_b_results.rds")

## ====================== (c) GlaucomaM calibration figure ===================
data("GlaucomaM", package="TH.data")
Xr2 <- as.matrix(GlaucomaM[, setdiff(names(GlaucomaM),"Class")])
y2  <- as.numeric(GlaucomaM$Class == "glaucoma")
X2  <- scale(Xr2); X2 <- X2[, apply(X2,2,function(c) all(is.finite(c)) & sd(c)>0), drop=FALSE]
n2 <- nrow(X2); p2 <- ncol(X2); X12 <- cbind(1,X2); D2 <- diag(c(0,rep(1,p2)))
set.seed(1); cv2 <- cv.glmnet(X2,y2,family="binomial",alpha=0)
l1se2 <- cv2$lambda.1se*n2
cat(sprintf("\n=== (c) GlaucomaM calibration at lambda.1se = %.1f (theory scale) ===\n", l1se2))
ft2 <- ridge_fit(X12,y2,l1se2,D2); ph2 <- ft2$pi
G <- 10
g2 <- pmin(ceiling(rank(ph2,ties.method="first")/(n2/G)), G)
idx <- split(seq_len(n2), g2)
obs  <- vapply(idx, function(I) mean(y2[I]), 0.0)
pred <- vapply(idx, function(I) mean(ph2[I]), 0.0)
ng   <- vapply(idx, length, 0L)
# Wilson 95% interval for the observed rate in each decile
wil <- function(k,n,z=1.96){ ph<-k/n; d<-1+z^2/n
  c((ph+z^2/(2*n)-z*sqrt(ph*(1-ph)/n+z^2/(4*n^2)))/d,
    (ph+z^2/(2*n)+z*sqrt(ph*(1-ph)/n+z^2/(4*n^2)))/d) }
CI <- t(mapply(function(I) wil(sum(y2[I]), length(I)), idx))
run2 <- mk_run(X12,y2,D2,n2); pv2 <- run2(l1se2,G,499,seed=110)
cat(sprintf("  naive dec=%.4f edge=%.4f | corrected dec=%.4f edge=%.4f\n",
    pv2["naive_dec"],pv2["naive_edge"],pv2["corr_dec"],pv2["corr_edge"]))
cal <- data.frame(decile=1:G, n=ng, pred_mean=round(pred,4),
                  obs_rate=round(obs,4), lo=round(CI[,1],4), hi=round(CI[,2],4))
print(cal, row.names=FALSE)
slope <- coef(glm(y2 ~ qlogis(pmin(pmax(ph2,1e-6),1-1e-6)), family=binomial))
cat(sprintf("  apparent calibration slope = %.3f ; CITL = %.3f ; AUC = %.3f\n",
    slope[2], slope[1],
    mean(outer(ph2[y2==1],ph2[y2==0],">")+0.5*outer(ph2[y2==1],ph2[y2==0],"=="))))

png("glaucoma_calibration.png", width=1500, height=1500, res=220)
par(mar=c(4.4,4.6,3.2,1.2))
plot(pred, obs, type="n", xlim=c(0,1), ylim=c(0,1), asp=1,
     xlab="Mean predicted probability within decile",
     ylab="Observed event rate within decile",
     main="GlaucomaM: ridge fit at the cross-validated penalty")
abline(0,1,col="grey35",lwd=2,lty=1)
grid(col="grey88", lty=1)
segments(pred, CI[,1], pred, CI[,2], col="grey55", lwd=2)
lines(pred, obs, col="#1f78b4", lwd=2)
points(pred, obs, pch=21, bg="#1f78b4", col="white", cex=1.5, lwd=1.4)
legend("topleft", bty="n", cex=0.82,
  legend=c(sprintf("n = %d, p = %d, kappa = %.2f", n2, p2, p2/n2),
           sprintf("lambda.1se = %.0f (theory scale)", l1se2),
           sprintf("G = 10 deciles, calibration slope %.2f", slope[2]),
           sprintf("naive p = %.4f (decile) / %.4f (EDGE)", pv2["naive_dec"], pv2["naive_edge"]),
           sprintf("corrected p = %.3f (decile) / %.3f (EDGE)", pv2["corr_dec"], pv2["corr_edge"])))
legend("bottomright", bty="n", cex=0.82, lty=c(1,1,1), lwd=c(2,2,2),
  col=c("grey35","#1f78b4","grey55"),
  legend=c("perfect calibration","observed by decile","95% Wilson interval"))
dev.off()
cat("  wrote glaucoma_calibration.png\n")
saveRDS(list(cal=cal, pvals=pv2, slope=slope, lam=l1se2), "T26_app_c_results.rds")
