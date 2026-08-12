# Does the blown-up beta.tilde make the bootstrap GENERATOR degenerate on GlaucomaM?
# If pi(beta.tilde) is pinned at 0/1, the bootstrap data are near-deterministic and the
# corrected p-values in the application cannot be trusted.
suppressPackageStartupMessages({library(glmnet); library(TH.data)})
data("GlaucomaM", package="TH.data")
X <- scale(as.matrix(GlaucomaM[, setdiff(names(GlaucomaM),"Class")]))
X <- X[, apply(X,2,function(c) all(is.finite(c)) & sd(c)>0), drop=FALSE]
y <- as.numeric(GlaucomaM$Class=="glaucoma"); n <- nrow(X); p <- ncol(X)
X1 <- cbind(1,X); D <- diag(c(0,rep(1,p)))
ridge_fit <- function(X1,y,lam,D,tol=1e-9,maxit=100){
  beta<-rep(0,ncol(X1))
  for(it in 1:maxit){eta<-drop(X1%*%beta);pi<-1/(1+exp(-eta))
    w<-pmax(pi*(1-pi),1e-8);z<-eta+(y-pi)/w
    bn<-drop(solve(crossprod(X1,w*X1)+lam*D,crossprod(X1,w*z)))
    if(max(abs(bn-beta))<tol){beta<-bn;break};beta<-bn}
  list(beta=beta,pi=1/(1+exp(-drop(X1%*%beta))))}
set.seed(1); cv <- cv.glmnet(X,y,family="binomial",alpha=0)
for (lam in c(cv$lambda.min*n, cv$lambda.1se*n)) {
  ft <- ridge_fit(X1,y,lam,D)
  pi <- pmin(pmax(ft$pi,1e-8),1-1e-8); w <- pmax(pi*(1-pi),1e-10)
  F <- crossprod(X1,w*X1)
  bt <- ft$beta + drop(solve(F, lam*D %*% ft$beta))
  gen <- as.numeric(1/(1+exp(-drop(X1 %*% bt))))
  cat(sprintf("\n=== lambda = %.1f ===\n", lam))
  cat(sprintf("  pi.hat  (fit)      : range [%.4f, %.4f]  sd %.4f\n", min(ft$pi), max(ft$pi), sd(ft$pi)))
  cat(sprintf("  pi(beta.tilde) GEN : range [%.2e, %.6f]  sd %.4f\n", min(gen), max(gen), sd(gen)))
  cat(sprintf("  GENERATOR pinned  : %d of %d below 1e-6 or above 1-1e-6  (%.0f%%)\n",
      sum(gen < 1e-6 | gen > 1-1e-6), n, 100*mean(gen < 1e-6 | gen > 1-1e-6)))
  cat(sprintf("  expected events under generator: %.1f  (observed %d)\n", sum(gen), sum(y)))
}
