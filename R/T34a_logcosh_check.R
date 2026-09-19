# T3.4a -- a smooth penalty that is not quadratic.
# IMPLEMENTATION CHECK, run before any size study is trusted.
#
# Penalty on each penalized coefficient (the intercept is never penalized):
#   pen(b)   = lambda * delta^2 * log cosh(b / delta)
#   a(b)     = pen'(b)  = lambda * delta * tanh(b / delta)     -> lambda * b as delta -> Inf (ridge)
#   K(b)     = pen''(b) = lambda / cosh(b / delta)^2           (0 < K <= lambda: convex)
# Near zero it is ridge; for |b| >> delta it is a lasso with weight lambda * delta.
# The correction follows the paper's argument with a(.) in place of K b:
#   beta.tilde = beta.hat + F^{-1} a(beta.hat),  mu.hat = U M^{-1} a(beta.tilde),  M = F + diag(K(beta.hat)).
#
# Check 1: with delta large the fit and both corrected statistics equal the ridge ones.
# Check 2: first-order mean and variance of r - mu.hat under a correct model, no bootstrap: the
#          mean should be near 0, as the paper checks for ridge.
# Check 3: choice of delta in design B, so that the penalty works in its non-quadratic range.
#
# v2: log cosh(x) is computed as log1p(2 sinh(x/2)^2) for |x| < 1 (the v1 formula lost precision
# to cancellation there), and Newton stops on the size of the Newton step, not the damped change.
suppressPackageStartupMessages(library(parallel))
G <- 10

lcosh <- function(x) { a <- abs(x); ifelse(a < 1, log1p(2*sinh(a/2)^2), a + log1p(exp(-2*a)) - log(2)) }
lc_fit <- function(X1, y, lam, delta, pen, tol=1e-9, maxit=200) {
  obj <- function(b) { eta <- drop(X1%*%b)
    sum(y*eta - (pmax(eta,0) + log1p(exp(-abs(eta))))) - lam*delta^2*sum(pen*lcosh(b/delta)) }
  beta <- rep(0, ncol(X1))
  for (it in 1:maxit) {
    eta <- drop(X1%*%beta); pi <- 1/(1+exp(-eta)); w <- pmax(pi*(1-pi), 1e-8)
    g <- drop(crossprod(X1, y-pi)) - pen*lam*delta*tanh(beta/delta)
    H <- crossprod(X1, w*X1) + diag(pen*lam/cosh(beta/delta)^2)
    step <- drop(solve(H, g))
    if (max(abs(step)) < tol) break
    f0 <- obj(beta); t <- 1
    for (h in 1:30) { if (obj(beta + t*step) >= f0 - 1e-10*max(1, abs(f0))) break; t <- t/2 }
    beta <- beta + t*step
  }
  list(beta=beta, pi=1/(1+exp(-drop(X1%*%beta))), iter=it)
}
ridge_fit <- function(X1,y,lam,D,tol=1e-10,maxit=200){
  beta <- rep(0,ncol(X1))
  for(it in 1:maxit){ eta<-drop(X1%*%beta); pi<-1/(1+exp(-eta))
    w<-pmax(pi*(1-pi),1e-8); z<-eta+(y-pi)/w
    bn<-drop(solve(crossprod(X1,w*X1)+lam*D, crossprod(X1,w*z)))
    if(max(abs(bn-beta))<tol){beta<-bn;break}; beta<-bn }
  list(beta=beta, pi=1/(1+exp(-drop(X1%*%beta))))
}
pieces_gen <- function(X1, y, fit, a_fun, K_fun) {
  pi<-pmin(pmax(fit$pi,1e-6),1-1e-6); w<-pmax(pi*(1-pi),1e-8)
  g<-pmin(ceiling(rank(pi,ties.method="first")/(length(y)/G)),G)
  idx<-split(seq_along(y),g); Vg<-vapply(idx,function(I) sum(w[I]),0.0)
  r<-(vapply(idx,function(I) sum(y[I]),0.0)-vapply(idx,function(I) sum(pi[I]),0.0))/sqrt(Vg)
  U<-t(vapply(idx,function(I) colSums(w[I]*X1[I,,drop=FALSE]),numeric(ncol(X1))))/sqrt(Vg)
  Fm<-crossprod(X1,w*X1)
  bt<-fit$beta+drop(solve(Fm, a_fun(fit$beta)))
  mu<-drop(U%*%solve(Fm+K_fun(fit$beta), a_fun(bt)))
  pbar<-vapply(idx,function(I) mean(pi[I]),0.0); v<-r-mu
  Z<-tryCatch(as.matrix(stats::poly(pbar,3)),error=function(e) NULL)
  ed<-if(is.null(Z)) NA_real_ else {zr<-crossprod(Z,v); as.numeric(t(zr)%*%solve(crossprod(Z))%*%zr)}
  list(dec=sum(v^2), edge=ed, r=r, v=v, bt=bt, trOm=G-sum(diag(U%*%solve(Fm,t(U)))))
}
designA <- function(seed) { set.seed(seed); n<-500; p<-5; b0<-c(0.5,-0.4,0.3,-0.3,0.2)
  X<-matrix(rnorm(n*p),n,p); list(X1=cbind(1,X), y=rbinom(n,1,1/(1+exp(-drop(X%*%b0)))), p=p) }
designB <- function(seed) { set.seed(seed); n<-400; p<-100; Ch<-chol(0.7^abs(outer(1:p,1:p,"-")))
  b0<-rep(c(.35,-.30,.25,-.20,.15),length.out=p)*0.8887
  X<-matrix(rnorm(n*p),n,p)%*%Ch; list(X1=cbind(1,X), y=rbinom(n,1,1/(1+exp(-drop(X%*%b0)))), p=p) }
lc_pieces <- function(d, lam, delta) { pen <- c(0, rep(1, d$p)); f <- lc_fit(d$X1, d$y, lam, delta, pen)
  list(f=f, s=pieces_gen(d$X1, d$y, f, function(b) pen*lam*delta*tanh(b/delta), function(b) diag(pen*lam/cosh(b/delta)^2))) }

cat("=== Check 1: log-cosh with large delta against ridge (same data, same lambda) ===\n")
for (des in list(list(f=designA, lam=137, nm="A"), list(f=designB, lam=416, nm="B"))) for (delta in c(50, 1e4)) {
  d <- des$f(1); D <- diag(c(0, rep(1, d$p)))
  L <- lc_pieces(d, des$lam, delta); fr <- ridge_fit(d$X1, d$y, des$lam, D)
  sr <- pieces_gen(d$X1, d$y, fr, function(b) drop(des$lam*D%*%b), function(b) des$lam*D)
  cat(sprintf("  design %s delta=%-6g: max|beta diff| %.2e | dec %.8f vs %.8f | edge %.8f vs %.8f | Newton iters %d\n",
      des$nm, delta, max(abs(L$f$beta-fr$beta)), L$s$dec, sr$dec, L$s$edge, sr$edge, L$f$iter))
}

cat("\n=== Check 3: design B, lambda = 416, choice of delta ===\n")
for (delta in c(0.02, 0.01, 0.005)) {
  L <- lc_pieces(designB(1), 416, delta); fr <- ridge_fit(designB(1)$X1, designB(1)$y, 416, diag(c(0, rep(1, 100))))
  cat(sprintf("  delta=%-6g | median |beta.hat|/delta %.2f | ||beta.hat|| log-cosh %.3f vs ridge %.3f | Newton iters %d\n",
      delta, median(abs(L$f$beta[-1]))/delta, sqrt(sum(L$f$beta[-1]^2)), sqrt(sum(fr$beta[-1]^2)), L$f$iter))
}

cat("\n=== Check 2: mean and covariance trace of r and r - mu.hat under a correct model ===\n")
check2 <- function(seed, des, lam, delta, kind) {
  d <- if (des == "A") designA(seed) else designB(seed); D <- diag(c(0, rep(1, d$p)))
  if (kind == "ridge") { f <- ridge_fit(d$X1, d$y, lam, D)
    s <- pieces_gen(d$X1, d$y, f, function(b) drop(lam*D%*%b), function(b) lam*D)
  } else { L <- lc_pieces(d, lam, delta); f <- L$f; s <- L$s }
  c(s$r, s$v, s$trOm, median(abs(f$beta[-1]))/delta)
}
nc <- max(1, detectCores()-2); cl <- makeCluster(nc)
clusterExport(cl, c("G","lcosh","lc_fit","ridge_fit","pieces_gen","designA","designB","lc_pieces"))
for (cfg in list(list(des="A", kind="ridge", lam=137, delta=1, N=2000), list(des="A", kind="logcosh", lam=137, delta=0.1, N=2000),
                 list(des="B", kind="ridge", lam=416, delta=1, N=1000), list(des="B", kind="logcosh", lam=416, delta=0.01, N=1000),
                 list(des="B", kind="logcosh", lam=416, delta=0.005, N=1000))) {
  M <- t(parSapply(cl, 20000 + seq_len(cfg$N), check2, des=cfg$des, lam=cfg$lam, delta=cfg$delta, kind=cfg$kind))
  R <- M[,1:G]; V <- M[,(G+1):(2*G)]
  cat(sprintf("  design %s %-7s lambda=%d delta=%-6g | ||mean r|| %.3f  ||mean(r-mu)|| %.3f | tr cov(r-mu) %.2f vs mean tr(I-UF^-1U') %.2f | median |beta|/delta %.2f\n",
      cfg$des, cfg$kind, cfg$lam, cfg$delta, sqrt(sum(colMeans(R)^2)), sqrt(sum(colMeans(V)^2)),
      sum(diag(cov(V))), mean(M[,2*G+1]), median(M[,2*G+2])))
}
stopCluster(cl)
