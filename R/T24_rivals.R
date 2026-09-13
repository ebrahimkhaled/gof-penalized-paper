# =====================================================================
# T2.4 — HEAD-TO-HEAD IN DGP B (n=400, p=100, AR(1) 0.7, lambda=416)
#        our prepivoted corrected test (decile + EDGE bases)
#        vs GRPtests::GRPtest  (Jankova, Shah, Buehlmann, Samworth)
#        vs RPtests::RPtest    (Shah & Buehlmann; Gaussian linear model - off-label here)
#
# PAIRED: every method sees the SAME simulated dataset in every replicate.
#
# Package notes (reported honestly in the write-up):
#  * The CRAN package is RPtests (plural), not "RPtest"; v0.1.5 installs cleanly.
#    Its RPtest() targets the GAUSSIAN linear model. Applying it to a binary
#    response is off-label; we run it anyway and let its SIZE row decide whether
#    it is usable as a comparator at all.
#  * GRPtests was REMOVED from CRAN (archived 2022-05-08). Installed from source
#    tarball GRPtests_0.1.2.tar.gz from the CRAN archive. It is the GLM
#    (binomial) descendant of RPtests and is the genuine rival.
#  * GRPtest is run at its package DEFAULT nsplits = 5 (the recommended,
#    aggregated version) AND at nsplits = 1 (single split), so the rival is not
#    handicapped by a bad tuning choice.
#
# REUSED VERBATIM from edge_basis_study.R: ridge_fit, pieces, p_prepivot2.
# =====================================================================
suppressPackageStartupMessages({library(glmnet); library(parallel)})
G <- 10; NB <- 149; LAM <- 416
B <- as.integer(Sys.getenv("T24_B", "200"))
NW <- as.integer(Sys.getenv("T24_NW", "9"))

## ---------- verbatim from edge_basis_study.R ----------
ridge_fit <- function(X1, y, lam, D, tol=1e-9, maxit=60) {
  beta <- rep(0, ncol(X1))
  for (it in 1:maxit) {
    eta <- drop(X1 %*% beta); pi <- 1/(1+exp(-eta))
    w <- pmax(pi*(1-pi), 1e-8); z <- eta + (y-pi)/w
    bn <- drop(solve(crossprod(X1, w*X1) + lam*D, crossprod(X1, w*z)))
    if (max(abs(bn-beta)) < tol) { beta <- bn; break }
    beta <- bn
  }
  eta <- drop(X1 %*% beta); list(beta=beta, pi=1/(1+exp(-eta)))
}
pieces <- function(X1, y, fit, lam, D) {
  pi <- pmin(pmax(fit$pi, 1e-6), 1-1e-6); w <- pmax(pi*(1-pi),1e-8)
  g  <- pmin(ceiling(rank(pi, ties.method="first")/(length(y)/G)), G)
  idx <- split(seq_along(y), g)
  Vg  <- vapply(idx, function(I) sum(w[I]), 0.0)
  og  <- vapply(idx, function(I) sum(y[I]), 0.0)
  eg  <- vapply(idx, function(I) sum(pi[I]), 0.0)
  pbar<- vapply(idx, function(I) mean(pi[I]), 0.0)
  r   <- (og - eg)/sqrt(Vg)
  U   <- t(vapply(idx, function(I) colSums(w[I]*X1[I,,drop=FALSE]), numeric(ncol(X1))))/sqrt(Vg)
  Fm  <- crossprod(X1, w*X1); K <- lam*D
  bt  <- fit$beta + drop(solve(Fm, K %*% fit$beta))
  mu  <- drop(U %*% solve(Fm+K, K %*% bt))
  rc  <- r - mu
  Z <- tryCatch(as.matrix(stats::poly(pbar, 3)), error=function(e) NULL)
  qf <- function(v, Z) { Zr <- crossprod(Z, v); as.numeric(t(Zr) %*% solve(crossprod(Z)) %*% Zr) }
  list(dec_c = sum(rc^2), edge_c = if (is.null(Z)) NA_real_ else qf(rc, Z),
       dec_p = sum(r^2),  edge_p = if (is.null(Z)) NA_real_ else qf(r, Z),
       Z = Z, bt = bt, Om = diag(G) - U %*% solve(Fm, t(U)))
}
p_prepivot2 <- function(X1, y, lam, D) {
  ft <- ridge_fit(X1,y,lam,D); st <- pieces(X1,y,ft,lam,D)
  gen <- pmin(pmax(as.numeric(1/(1+exp(-X1 %*% st$bt))),1e-6),1-1e-6)
  Sd <- Se <- numeric(NB)
  for (bb in 1:NB) {
    ys <- rbinom(nrow(X1),1,gen)
    s2 <- pieces(X1,ys,ridge_fit(X1,ys,lam,D),lam,D)
    Sd[bb] <- s2$dec_c; Se[bb] <- s2$edge_c
  }
  c(dec = (1+sum(Sd>=st$dec_c))/(NB+1),
    edge= (1+sum(Se>=st$edge_c, na.rm=TRUE))/(NB+1))
}
## ------------------------------------------------------

## DGP B, frozen (identical to oracle_knownnull.R / edge_basis_study.R)
nB<-400; pB<-100; SigB<-0.7^abs(outer(1:pB,1:pB,"-")); ChB<-chol(SigB)
bB<-rep(c(.35,-.30,.25,-.20,.15),length.out=pB)*0.8887; DB<-diag(c(0,rep(1,pB)))
vB<-as.numeric(t(bB)%*%SigB%*%bB)
mkB <- function(g,k){ X<-matrix(rnorm(nB*pB),nB,pB)%*%ChB; e0<-drop(X%*%bB)
  eta <- if(k=="cubic") e0+g*(e0^3-3*vB*e0)/sqrt(6*vB^3) else if(k=="quad") e0+g*(X[,1]^2-1) else e0
  list(X=X, X1=cbind(1,X), y=rbinom(nB,1,1/(1+exp(-eta))), eta0=e0) }

one_rep <- function(b, gam, kind) {
  suppressPackageStartupMessages({library(GRPtests); library(RPtests)})
  d <- mkB(gam, kind)
  t0 <- proc.time()[3]; pp <- p_prepivot2(d$X1, d$y, LAM, DB); t_ours <- proc.time()[3]-t0
  t0 <- proc.time()[3]
  g5 <- tryCatch(as.numeric(GRPtest(d$X, d$y, fam="binomial", nsplits=5L)), error=function(e) NA_real_)
  t_g5 <- proc.time()[3]-t0
  t0 <- proc.time()[3]
  g1 <- tryCatch(as.numeric(GRPtest(d$X, d$y, fam="binomial", nsplits=1L)), error=function(e) NA_real_)
  t_g1 <- proc.time()[3]-t0
  t0 <- proc.time()[3]
  rp <- tryCatch(as.numeric(RPtests::RPtest(d$X, d$y, test="nonlin", B=49L)), error=function(e) NA_real_)
  t_rp <- proc.time()[3]-t0
  c(rep=b, gamma=gam, ours_dec=pp[["dec"]], ours_edge=pp[["edge"]],
    grp5=g5, grp1=g1, rptest=rp,
    t_ours=t_ours, t_grp5=t_g5, t_grp1=t_g1, t_rptest=t_rp)
}

cells <- list(
  list(k="null",  g=0.0),
  list(k="quad",  g=1.0),
  list(k="quad",  g=1.5),
  list(k="quad",  g=0.5),
  list(k="cubic", g=1.5),
  list(k="cubic", g=1.0)
)

cl <- makeCluster(NW)
clusterSetRNGStream(cl, 2026)
clusterExport(cl, c("G","NB","LAM","ridge_fit","pieces","p_prepivot2",
                    "nB","pB","ChB","bB","DB","vB","mkB","one_rep"))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(glmnet))))

allp <- NULL; summ <- NULL
T0 <- Sys.time()
for (ce in cells) {
  tc <- Sys.time()
  L <- parLapply(cl, 1:B, one_rep, gam=ce$g, kind=ce$k)
  M <- as.data.frame(do.call(rbind, L)); M$alt <- ce$k
  allp <- rbind(allp, M)
  rj <- function(v) mean(v < 0.05, na.rm=TRUE)
  se <- function(v) sqrt(rj(v)*(1-rj(v))/sum(!is.na(v)))
  s <- data.frame(alt=ce$k, gamma=ce$g, B=B,
    ours_dec=rj(M$ours_dec),  se_dec =se(M$ours_dec),
    ours_edge=rj(M$ours_edge),se_edge=se(M$ours_edge),
    grp5=rj(M$grp5),          se_grp5=se(M$grp5),
    grp1=rj(M$grp1),          se_grp1=se(M$grp1),
    rptest=rj(M$rptest),      se_rp  =se(M$rptest),
    na_grp5=sum(is.na(M$grp5)), na_grp1=sum(is.na(M$grp1)), na_rp=sum(is.na(M$rptest)),
    t_ours=mean(M$t_ours), t_grp5=mean(M$t_grp5), t_grp1=mean(M$t_grp1), t_rp=mean(M$t_rptest))
  summ <- rbind(summ, s)
  cat(sprintf("[%5.1f min] %-5s g=%.1f | OURS dec=%.3f edge=%.3f | GRP5=%.3f GRP1=%.3f | RPtest=%.3f\n",
      as.numeric(difftime(Sys.time(),T0,units="mins")), ce$k, ce$g,
      s$ours_dec, s$ours_edge, s$grp5, s$grp1, s$rptest)); flush.console()
  saveRDS(list(summary=summ, pvalues=allp, B=B, NB=NB, lambda=LAM, workers=NW),
          "T24_rivals_results.rds")
  write.csv(allp, "T24_rivals_pvalues.csv", row.names=FALSE)
}
stopCluster(cl)

cat("\n===== T2.4 HEAD-TO-HEAD, DGP B (n=400, p=100, lambda=416), B =", B, " =====\n")
print(summ, row.names=FALSE, digits=3)
cat("\nTotal wall clock:", round(as.numeric(difftime(Sys.time(),T0,units="mins")),1), "min\n")
