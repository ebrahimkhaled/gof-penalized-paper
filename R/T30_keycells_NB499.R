# T3.0 -- BJ revision, co-author comment 3: are the size results sensitive to the number of
# bootstrap replicates? The decisive cells of Table 2 are re-run with NB = 499 on EXACTLY the
# datasets of the original runs: same seeds, same code, same RNG kind (clusterSetRNGStream sets
# L'Ecuyer-CMRG on the workers before each set.seed, as in T23_rerun.R and verify_E13.R). The
# first 149 bootstrap draws are therefore the original ones, and the NB = 149 p-values must
# reproduce Table 2 and Table S6 -- a built-in check that the archived code gives the paper.
#   design A, lambda = 137 : seeds 1001..3000        (T23_rerun.R)
#   design B, lambda = 416 : seeds 1001..2000        (T23_rerun.R)
#   design B, lambda = 50  : seeds 1001..2000        (T23_rerun.R, block 1)
#                            seeds 1..1500           (verify_E13.R, block 2)
#                            seeds 100001..101500    (verify_E13.R, block 3)
# Rules: "p < 0.05" is the one Table 2 used. At NB = 499, alpha * (NB + 1) = 25 is an integer, so
# "p <= 0.05" has exact level 0.05; at NB = 149 the two rules coincide (attainable 7/150).
# Output: T30_keycells_pvalues.csv (per replicate, p at NB = 149 and 499), T30_keycells_results.rds.
suppressPackageStartupMessages({library(glmnet); library(parallel)})
G <- 10; NB <- 499

worker <- function(seed, panel, lam) {
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
  if (panel=="A") {
    n<-500; p<-5; b0<-c(0.5,-0.4,0.3,-0.3,0.2); D<-diag(c(0,rep(1,p)))
    X<-matrix(rnorm(n*p),n,p); X1<-cbind(1,X); y<-rbinom(n,1,1/(1+exp(-drop(X%*%b0))))
  } else {
    n<-400; p<-100; S<-0.7^abs(outer(1:p,1:p,"-")); Ch<-chol(S)
    b0<-rep(c(.35,-.30,.25,-.20,.15),length.out=p)*0.8887; D<-diag(c(0,rep(1,p)))
    X<-matrix(rnorm(n*p),n,p)%*%Ch; X1<-cbind(1,X); y<-rbinom(n,1,1/(1+exp(-drop(X%*%b0))))
  }
  ft<-ridge_fit(X1,y,lam,D); st<-pieces(X1,y,ft,lam,D)
  gen<-pmin(pmax(as.numeric(1/(1+exp(-X1%*%st$bt))),1e-6),1-1e-6)
  Sd<-Se<-numeric(NB)
  for(bb in 1:NB){ ys<-rbinom(n,1,gen); s2<-pieces(X1,ys,ridge_fit(X1,ys,lam,D),lam,D)
    Sd[bb]<-s2$dec; Se[bb]<-s2$edge }
  c(seed=seed,
    dec149=(1+sum(Sd[1:149]>=st$dec))/150, edge149=(1+sum(Se[1:149]>=st$edge,na.rm=TRUE))/150,
    dec499=(1+sum(Sd>=st$dec))/(NB+1), edge499=(1+sum(Se>=st$edge,na.rm=TRUE))/(NB+1))
}

nc <- max(1, detectCores() - 2)
cl <- makeCluster(nc); on.exit(stopCluster(cl))
clusterEvalQ(cl, suppressPackageStartupMessages(library(glmnet)))
clusterExport(cl, c("G","NB"))
clusterSetRNGStream(cl, 20260807)          # the same stream set-up as T23_rerun.R and verify_E13.R
cat("T3.0 key cells at NB =", NB, "on", nc, "workers\n")

cells <- list(
  list(p="A", lam=137, seeds=1001:3000,      block="T23",   tab2=c(dec=0.0475, edge=0.0505)),
  list(p="B", lam=416, seeds=1001:2000,      block="T23",   tab2=c(dec=0.0710, edge=0.0540)),
  list(p="B", lam=50,  seeds=1001:2000,      block="blk1",  tab2=c(dec=0.0630, edge=0.0270)),
  list(p="B", lam=50,  seeds=1:1500,         block="blk2",  tab2=c(dec=0.0587, edge=0.0440)),
  list(p="B", lam=50,  seeds=100001:101500,  block="blk3",  tab2=c(dec=0.0580, edge=0.0367)))
all <- NULL; summ <- NULL; t0 <- Sys.time()
rate <- function(x, rule) { x <- x[!is.na(x)]; m <- if (rule=="lt") mean(x < 0.05) else mean(x <= 0.05)
  c(rate=round(m,4), se=round(sqrt(m*(1-m)/length(x)),4)) }
for (ce in cells) {
  R <- t(parSapply(cl, ce$seeds, worker, panel=ce$p, lam=ce$lam))
  d <- data.frame(panel=ce$p, lambda=ce$lam, block=ce$block, R); all <- rbind(all, d)
  write.csv(all, "T30_keycells_pvalues.csv", row.names=FALSE)
  s <- data.frame(panel=ce$p, lambda=ce$lam, block=ce$block, B=nrow(d),
    dec149=rate(d$dec149,"lt")[["rate"]], edge149=rate(d$edge149,"lt")[["rate"]],
    tab_dec=ce$tab2[["dec"]], tab_edge=ce$tab2[["edge"]],
    dec499_lt=rate(d$dec499,"lt")[["rate"]], se_dec499=rate(d$dec499,"lt")[["se"]],
    edge499_lt=rate(d$edge499,"lt")[["rate"]], se_edge499=rate(d$edge499,"lt")[["se"]],
    dec499_le=rate(d$dec499,"le")[["rate"]], edge499_le=rate(d$edge499,"le")[["rate"]],
    flip_dec=sum((d$dec149<0.05) != (d$dec499<0.05)), flip_edge=sum((d$edge149<0.05) != (d$edge499<0.05), na.rm=TRUE))
  summ <- rbind(summ, s); print(s, row.names=FALSE)
  cat(sprintf("  reproduces Table 2 at NB = 149: %s | %.0f min\n",
      isTRUE(all.equal(c(s$dec149, s$edge149), c(s$tab_dec, s$tab_edge), tolerance=1e-8)),
      as.numeric(difftime(Sys.time(),t0,units="mins"))))
  saveRDS(list(pvalues=all, summary=summ), "T30_keycells_results.rds")
}
b50 <- all[all$panel=="B" & all$lambda==50, ]
pooled <- data.frame(panel="B", lambda=50, block="pooled", B=nrow(b50),
  dec149=rate(b50$dec149,"lt")[["rate"]], edge149=rate(b50$edge149,"lt")[["rate"]],
  dec499_lt=rate(b50$dec499,"lt")[["rate"]], se_dec499=rate(b50$dec499,"lt")[["se"]],
  edge499_lt=rate(b50$edge499,"lt")[["rate"]], se_edge499=rate(b50$edge499,"lt")[["se"]],
  dec499_le=rate(b50$dec499,"le")[["rate"]], edge499_le=rate(b50$edge499,"le")[["rate"]])
cat("\n===== T3.0 KEY CELLS: NB = 149 (Table 2) against NB = 499 =====\n"); print(summ, row.names=FALSE)
cat("\npooled lambda = 50 cell (Table 2 reports decile 0.0595, EDGE 0.0370):\n"); print(pooled, row.names=FALSE)
saveRDS(list(pvalues=all, summary=summ, pooled50=pooled), "T30_keycells_results.rds")
cat("runtime:", round(as.numeric(difftime(Sys.time(),t0,units="mins")),1), "min\n")
