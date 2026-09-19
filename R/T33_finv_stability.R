# T3.3 -- how much do the glaucoma p-values of Table 5 depend
# on inverting a near-singular F? The whole 4 x 3 grid is re-run with F^{-1} replaced by an
# eigenvalue-floored inverse (eigenvalues below tau * largest raised to that floor), tau in
# {0 (plain solve, as in the paper), 1e-9, 1e-7, 1e-5}. Common random numbers: each cell draws
# one n x NB matrix of uniforms and every variant uses it, so only the inverse differs.
# The floored inverse is used everywhere F^{-1} enters: beta.tilde, and every bootstrap refit.
# Output: T33_finv_stability.csv and T33_finv_stability_results.rds.
suppressPackageStartupMessages({library(glmnet); library(TH.data)})
data("GlaucomaM", package="TH.data")
Xr <- as.matrix(GlaucomaM[, setdiff(names(GlaucomaM),"Class")])
y  <- as.numeric(GlaucomaM$Class == "glaucoma")
X  <- scale(Xr); X <- X[, apply(X,2,function(c) all(is.finite(c)) & sd(c)>0), drop=FALSE]
n <- nrow(X); p <- ncol(X); X1 <- cbind(1,X); D <- diag(c(0,rep(1,p)))
NB <- 499; taus <- c(0, 1e-9, 1e-7, 1e-5)

ridge_fit <- function(X1,y,lam,D,tol=1e-9,maxit=200){
  beta <- rep(0,ncol(X1))
  for(it in 1:maxit){ eta<-drop(X1%*%beta); pi<-1/(1+exp(-eta))
    w<-pmax(pi*(1-pi),1e-8); z<-eta+(y-pi)/w
    bn<-drop(solve(crossprod(X1,w*X1)+lam*D, crossprod(X1,w*z)))
    if(max(abs(bn-beta))<tol){beta<-bn;break}; beta<-bn }
  list(beta=beta, pi=1/(1+exp(-drop(X1%*%beta))))
}
# F^{-1} b, with the eigenvalues of F floored at tau * (largest eigenvalue); tau = 0 is solve()
finv_b <- function(Fm, b, tau) {
  if (tau == 0) return(drop(solve(Fm, b)))
  e <- eigen(Fm, symmetric=TRUE)
  d <- pmax(e$values, tau*e$values[1])
  drop(e$vectors %*% (crossprod(e$vectors, b) / d))
}
pieces <- function(X1,y,fit,lam,D,G,tau){
  pi<-pmin(pmax(fit$pi,1e-6),1-1e-6); w<-pmax(pi*(1-pi),1e-8)
  g<-pmin(ceiling(rank(pi,ties.method="first")/(length(y)/G)),G)
  idx<-split(seq_along(y),g); Vg<-vapply(idx,function(I) sum(w[I]),0.0)
  r<-(vapply(idx,function(I) sum(y[I]),0.0)-vapply(idx,function(I) sum(pi[I]),0.0))/sqrt(Vg)
  U<-t(vapply(idx,function(I) colSums(w[I]*X1[I,,drop=FALSE]),numeric(ncol(X1))))/sqrt(Vg)
  Fm<-crossprod(X1,w*X1); K<-lam*D
  bt<-fit$beta+finv_b(Fm,drop(K%*%fit$beta),tau); mu<-drop(U%*%solve(Fm+K,K%*%bt))
  pbar<-vapply(idx,function(I) mean(pi[I]),0.0); v<-r-mu
  Z<-tryCatch(as.matrix(stats::poly(pbar,3)),error=function(e) NULL)
  qf<-function(u,Z){zr<-crossprod(Z,u); as.numeric(t(zr)%*%solve(crossprod(Z))%*%zr)}
  list(Sc=sum(v^2), Sce=if(is.null(Z)) NA else qf(v,Z), bt=bt, Fm=Fm)
}

set.seed(1); cv <- cv.glmnet(X,y,family="binomial",alpha=0)
lmin <- cv$lambda.min*n; l1se <- cv$lambda.1se*n
lams <- c("lambda.min"=lmin, "0.5*1se"=l1se/2, "lambda.1se"=l1se, "2*1se"=2*l1se)
cat(sprintf("GlaucomaM n=%d p=%d; lambda.min %.1f, lambda.1se %.1f (theory scale)\n", n, p, lmin, l1se))

out <- NULL; t0 <- Sys.time()
for (ln in names(lams)) {
  lam <- lams[[ln]]; ft <- ridge_fit(X1,y,lam,D)
  pi <- pmin(pmax(ft$pi,1e-6),1-1e-6); w <- pmax(pi*(1-pi),1e-8)
  kF <- kappa(crossprod(X1,w*X1), exact=TRUE)
  gen0 <- NULL
  for (G in c(5,10,20)) {
    set.seed(7000 + G*10 + match(ln, names(lams)))
    Umat <- matrix(runif(n*NB), n, NB)
    for (tau in taus) {
      st <- pieces(X1,y,ft,lam,D,G,tau)
      gen <- pmin(pmax(as.numeric(1/(1+exp(-X1%*%st$bt))),1e-6),1-1e-6)
      if (tau == 0) gen0 <- gen
      Sd <- Se <- numeric(NB)
      for (b in 1:NB) { ys <- as.numeric(Umat[,b] < gen)
        s2 <- pieces(X1,ys,ridge_fit(X1,ys,lam,D),lam,D,G,tau); Sd[b] <- s2$Sc; Se[b] <- s2$Sce }
      out <- rbind(out, data.frame(lambda=ln, lam=round(lam,1), G=G, tau=tau, kappa_F=signif(kF,3),
        p_dec=(1+sum(Sd>=st$Sc))/(NB+1), p_edge=(1+sum(Se>=st$Sce,na.rm=TRUE))/(NB+1),
        norm_ratio=round(sqrt(sum(st$bt[-1]^2))/sqrt(sum(ft$beta[-1]^2)),1),
        gen_min=signif(min(gen),3), gen_max=signif(max(gen),4), gen_sd=round(sd(gen),4),
        pinned=sum(gen<=1e-6 | gen>=1-1e-6), events_gen=round(sum(gen),1), events_obs=sum(y),
        max_gen_shift=signif(max(abs(gen-gen0)),3)))
      cat(sprintf("  %-10s G=%2d tau=%-6g | p dec %.3f edge %.3f | ||bt||/||bh|| %7.1f | gen [%.2e, %.4f] | shift %.2e | %.1f min\n",
          ln, G, tau, tail(out$p_dec,1), tail(out$p_edge,1), tail(out$norm_ratio,1), min(gen), max(gen),
          tail(out$max_gen_shift,1), as.numeric(difftime(Sys.time(),t0,units="mins"))))
    }
  }
  write.csv(out, "T33_finv_stability.csv", row.names=FALSE)
}
base <- out[out$tau==0, c("lambda","G","p_dec","p_edge")]
cmp <- merge(out[out$tau>0,], base, by=c("lambda","G"), suffixes=c("","_plain"))
cmp$d_dec <- abs(cmp$p_dec - cmp$p_dec_plain); cmp$d_edge <- abs(cmp$p_edge - cmp$p_edge_plain)
cat("\n===== largest change in the p-value against the plain inverse, by floor =====\n")
print(aggregate(cbind(d_dec, d_edge, max_gen_shift) ~ tau, data=cmp, FUN=max))
saveRDS(list(grid=out, comparison=cmp, lmin=lmin, l1se=l1se), "T33_finv_stability_results.rds")
cat("runtime:", round(as.numeric(difftime(Sys.time(),t0,units="mins")),1), "min\n")
