# E13, DECISIVE TEST: is the generator overshoot the CAUSE of EDGE's conservatism at
# design B, lambda = 50, or just a correlate?
#
# 8h(b) measured a 26% overshoot in ||beta.tilde|| in that cell and <=5% in all five
# others, and that cell is the only conservative one. But E[S*]/E[S] = 0.988 there, so
# the mean does not carry the effect. Two things settle it:
#   (1) var(S*)/var(S)  -- does the bootstrap world have a FATTER right tail?
#   (2) INTERVENTION: rerun the same cell with the debiasing SHRUNK back so the
#       generator has the right scale, holding everything else fixed. If the
#       conservatism disappears, the overshoot caused it. If it survives, it did not.
# The intervention is the real test: it changes the suspected cause and nothing else.
suppressPackageStartupMessages({library(parallel)})
G <- 10; NB <- 149; REPS <- 400

design <- function(panel) {
  if (panel == "A") {
    n <- 500; p <- 5; b0 <- c(0.5,-0.4,0.3,-0.3,0.2); X <- matrix(rnorm(n*p), n, p)
  } else {
    n <- 400; p <- 100; S <- 0.7^abs(outer(1:p,1:p,"-")); Ch <- chol(S)
    b0 <- rep(c(.35,-.30,.25,-.20,.15), length.out = p) * 0.8887
    X <- matrix(rnorm(n*p), n, p) %*% Ch
  }
  list(X1 = cbind(1,X), y = rbinom(n,1,1/(1+exp(-drop(X %*% b0)))),
       D = diag(c(0,rep(1,p))), b0 = c(0,b0), n = n, p = p)
}
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

# shrink: 1 = the procedure as published (beta.tilde). s < 1 pulls the GENERATOR back
# toward beta.hat: gen.beta = beta.hat + s*(beta.tilde - beta.hat). The statistic and
# mu.hat are UNCHANGED -- only the bootstrap generator moves.
worker <- function(seed, panel, lam, shrink) {
  set.seed(seed)
  d  <- design(panel)
  ft <- ridge_fit(d$X1, d$y, lam, d$D)
  st <- pieces(d$X1, d$y, ft, lam, d$D)
  gb  <- ft$beta + shrink * (st$bt - ft$beta)
  gen <- pmin(pmax(as.numeric(1/(1+exp(-d$X1 %*% gb))), 1e-6), 1-1e-6)
  Sd <- Se <- numeric(NB)
  for (b in 1:NB) {
    ys <- rbinom(d$n, 1, gen)
    s2 <- pieces(d$X1, ys, ridge_fit(d$X1, ys, lam, d$D), lam, d$D)
    Sd[b] <- s2$dec; Se[b] <- s2$edge
  }
  c(p_dec  = (1+sum(Sd >= st$dec))/(NB+1),
    p_edge = (1+sum(Se >= st$edge, na.rm=TRUE))/(NB+1),
    vr_dec = var(Sd)/1, vr_edge = var(Se, na.rm=TRUE),
    S_dec = st$dec, S_edge = st$edge,
    mS_dec = mean(Sd), mS_edge = mean(Se, na.rm=TRUE),
    ngen = sqrt(sum(gb^2)), nb0 = sqrt(sum(d$b0^2)))
}

nc <- max(1, detectCores() - 2)
cl <- makeCluster(nc); on.exit(stopCluster(cl))
clusterExport(cl, c("G","NB","design","ridge_fit","pieces"))
clusterSetRNGStream(cl, 20260809)
cat("E13 intervention test on", nc, "workers; REPS =", REPS, "NB =", NB, "\n\n")

run <- function(panel, lam, shrink, tag) {
  t0 <- Sys.time()
  M <- simplify2array(parLapply(cl, seq_len(REPS), worker,
                                panel=panel, lam=lam, shrink=shrink))
  rd <- mean(M["p_dec",] <= 0.05); re <- mean(M["p_edge",] <= 0.05)
  se <- sqrt(re*(1-re)/REPS)
  # variance of the bootstrap statistic vs variance of the observed statistic across reps
  vr_e <- mean(M["vr_edge",]) / var(M["S_edge",], na.rm=TRUE)
  vr_d <- mean(M["vr_dec",])  / var(M["S_dec",],  na.rm=TRUE)
  cat(sprintf("%-34s | size dec %.4f  edge %.4f (SE %.4f) | var(S*)/var(S): dec %.2f edge %.2f | gen/truth %.2f  (%.0fs)\n",
      tag, rd, re, se, vr_d, vr_e,
      mean(M["ngen",])/mean(M["nb0",]), as.numeric(difftime(Sys.time(),t0,units="secs"))))
  c(size_dec=rd, size_edge=re, se_edge=se, vr_dec=vr_d, vr_edge=vr_e,
    gen_ratio=mean(M["ngen",])/mean(M["nb0",]))
}

out <- list()
out$B50_published  <- run("B",   50, 1.00, "B lam=50   PUBLISHED (shrink=1)")
## shrink chosen so the generator lands on the TRUTH's scale: with ||bhat||=0.84,
## ||btilde||=2.90, ||b0||=2.31, collinearity gives 0.84 + s(2.90-0.84) = 2.31 => s ~ 0.71.
## The printed gen/truth column verifies it landed there.
out$B50_shrunk     <- run("B",   50, 0.71, "B lam=50   INTERVENED (shrink=0.71)")
out$B1000_published<- run("B", 1000, 1.00, "B lam=1000 PUBLISHED (control)")
saveRDS(out, "overshoot_test_results.rds")
cat("\nsaved overshoot_test_results.rds\n")
