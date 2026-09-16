# T3.2 -- BJ revision (co-author comment 1; the application's support): size of the corrected
# test on the application's OWN design. The GlaucomaM predictors (n = 196, p = 62, pairwise
# correlations up to 0.996) are held fixed and the response is simulated from a logistic model,
# so every rejection is false. lambda.1se is chosen by 10-fold cross-validation on each simulated
# response and then held fixed in the bootstrap, as in the application. Two truths:
#   T1  pi(beta.hat), the ridge fit to the real response at lambda.min
#   T2  pi(beta.tilde) at lambda.1se, the generator the application's own bootstrap drew from
# Reports the prepivoted corrected tests (NB = 499) and the uncorrected tests referred to the
# maximum likelihood covariance by Monte Carlo (4000 draws), decile and EDGE bases, G = 10.
# Output: T32_glaucoma_design_pvalues.csv and T32_glaucoma_design_results.rds.
suppressPackageStartupMessages({library(glmnet); library(TH.data); library(parallel)})
data("GlaucomaM", package="TH.data")
Xr <- as.matrix(GlaucomaM[, setdiff(names(GlaucomaM),"Class")])
yobs <- as.numeric(GlaucomaM$Class == "glaucoma")
X  <- scale(Xr); X <- X[, apply(X,2,function(c) all(is.finite(c)) & sd(c)>0), drop=FALSE]
n <- nrow(X); p <- ncol(X); X1 <- cbind(1,X); D <- diag(c(0,rep(1,p)))
G <- 10; NB <- 499; REPS <- 1000

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
       Om=diag(G)-U%*%solve(Fm,t(U)), Z=Z, bt=bt)
}
mc_p <- function(S,Om,Z=NULL,nd=4000){
  R<-chol(Om+diag(1e-9,nrow(Om))); W<-matrix(rnorm(nd*nrow(Om)),nd)%*%R
  if(is.null(Z)) return(mean(rowSums(W^2)>=S))
  A<-solve(crossprod(Z)); mean(apply(W,1,function(u){zr<-crossprod(Z,u); as.numeric(t(zr)%*%A%*%zr)})>=S)
}

# ---- the two truths, from the real data exactly as the application computed them ------------
set.seed(1); cv0 <- cv.glmnet(X, yobs, family="binomial", alpha=0)
lmin <- cv0$lambda.min*n; l1se <- cv0$lambda.1se*n
f_min <- ridge_fit(X1, yobs, lmin, D)
f_1se <- ridge_fit(X1, yobs, l1se, D); st0 <- pieces(X1, yobs, f_1se, l1se, D, G)
PI0 <- cbind(T1 = f_min$pi, T2 = as.numeric(1/(1+exp(-X1%*%st0$bt))))
cat(sprintf("real data: lambda.min %.1f, lambda.1se %.1f (theory scale)\n", lmin, l1se))
for (k in 1:2) cat(sprintf("  truth %s: probabilities in [%.2e, %.4f], sd %.3f, expected events %.1f\n",
    colnames(PI0)[k], min(PI0[,k]), max(PI0[,k]), sd(PI0[,k]), sum(PI0[,k])))

worker <- function(seed, truth) {
  set.seed(seed)
  y <- rbinom(n, 1, PI0[, truth])
  cv <- cv.glmnet(X, y, family="binomial", alpha=0, nfolds=10)
  lam <- cv$lambda.1se*n
  ft <- ridge_fit(X1,y,lam,D); st <- pieces(X1,y,ft,lam,D,G)
  gen <- pmin(pmax(as.numeric(1/(1+exp(-X1%*%st$bt))),1e-6),1-1e-6)
  Sd <- Se <- numeric(NB)
  for (b in 1:NB) { ys <- rbinom(n,1,gen); s2 <- pieces(X1,ys,ridge_fit(X1,ys,lam,D),lam,D,G)
    Sd[b] <- s2$Sc; Se[b] <- s2$Sce }
  c(seed=seed, lam=lam, events=sum(y),
    p_dec=(1+sum(Sd>=st$Sc))/(NB+1), p_edge=(1+sum(Se>=st$Sce,na.rm=TRUE))/(NB+1),
    pu_dec=mc_p(st$Sp,st$Om), pu_edge=if(is.null(st$Z)) NA_real_ else mc_p(st$Spe,st$Om,st$Z),
    pinned=mean(gen<=1e-6 | gen>=1-1e-6), gen_sd=sd(gen))
}

nc <- max(1, detectCores() - 2)
cl <- makeCluster(nc); on.exit(stopCluster(cl))
clusterEvalQ(cl, suppressPackageStartupMessages(library(glmnet)))
clusterExport(cl, c("G","NB","n","X","X1","D","PI0","ridge_fit","pieces","mc_p"))
clusterSetRNGStream(cl, 20260917)          # RNG kind as in the other studies; each replicate is set.seed()
cat("T3.2 glaucoma design on", nc, "workers, NB =", NB, "\n")

t0 <- Sys.time(); all <- NULL; summ <- NULL
for (tr in list(list(k="T1", seeds=400000L + seq_len(REPS)), list(k="T2", seeds=410000L + seq_len(REPS)))) {
  R <- t(parSapply(cl, tr$seeds, worker, truth=tr$k))
  d <- data.frame(truth=tr$k, R); all <- rbind(all, d)
  write.csv(all, "T32_glaucoma_design_pvalues.csv", row.names=FALSE)
  rate <- function(x) { x <- x[!is.na(x)]; c(round(mean(x < 0.05),4), round(sqrt(mean(x<0.05)*(1-mean(x<0.05))/length(x)),4)) }
  s <- data.frame(truth=tr$k, B=nrow(d), lam_med=round(median(d$lam),1),
    lam_q25=round(quantile(d$lam,.25),1), lam_q75=round(quantile(d$lam,.75),1),
    corr_dec=rate(d$p_dec)[1], se_corr_dec=rate(d$p_dec)[2],
    corr_edge=rate(d$p_edge)[1], se_corr_edge=rate(d$p_edge)[2],
    unc_dec=rate(d$pu_dec)[1], unc_edge=rate(d$pu_edge)[1],
    pinned_mean=round(mean(d$pinned),4), edge_na=sum(is.na(d$p_edge)))
  summ <- rbind(summ, s); print(s, row.names=FALSE)
  cat(sprintf("  %.0f min\n", as.numeric(difftime(Sys.time(),t0,units="mins"))))
  saveRDS(list(pvalues=all, summary=summ, lmin=lmin, l1se=l1se), "T32_glaucoma_design_results.rds")
}
cat("\n===== T3.2 SIZE ON THE GLAUCOMA DESIGN =====\n"); print(summ, row.names=FALSE)
cat("runtime:", round(as.numeric(difftime(Sys.time(),t0,units="mins")),1), "min\n")
