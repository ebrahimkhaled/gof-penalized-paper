# T2.5 — POWER AT HIGH PRECISION. Supersedes manifest 6.1 and 6.2 (both B = 200,
# MC SE <= 0.035). Same designs, same lambdas, same departures, same estimator and
# basis definitions as edge_basis_study.R / power_study.R -- ONLY B changes, so the
# numbers are directly comparable.
#
#   Panel A: n=500,  p=5,   lambda=100, departures gamma*(x1^2-1) AND gamma*x1*x2
#   Panel B: n=400,  p=100, AR(1) 0.7, lambda=416, departure gamma*(x1^2-1)
#
# Adds three things the B=200 runs could not support:
#   (1) per-replication p-values written to CSV (not just rejection rates), so any
#       cell can be re-derived and any MC SE recomputed by a reader;
#   (2) SIZE-ADJUSTED power -- the gamma=0 cell of each panel supplies the null
#       p-value distribution, and the alternative cells are thresholded at its 5th
#       percentile. E13 requires this for panel B, where EDGE under-rejects.
#   (3) the interaction departure in panel A, which is half of the eight-alternative
#       paired study behind the "-0.010 power cost" claim.
suppressPackageStartupMessages(library(parallel))
G <- 10; NB <- 149; ND <- 4000; B <- 2000

ridge_fit <- function(X1, y, lam, D, tol=1e-9, maxit=60) {
  beta <- rep(0, ncol(X1))
  for (it in 1:maxit) {
    eta <- drop(X1 %*% beta); pi <- 1/(1+exp(-eta))
    w <- pmax(pi*(1-pi), 1e-8); z <- eta + (y-pi)/w
    bn <- drop(solve(crossprod(X1, w*X1) + lam*D, crossprod(X1, w*z)))
    if (max(abs(bn-beta)) < tol) { beta <- bn; break }
    beta <- bn
  }
  list(beta=beta, pi=1/(1+exp(-drop(X1 %*% beta))))
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
p_mle2 <- function(X1, y, D) {
  ft <- ridge_fit(X1,y,0,D); st <- pieces(X1,y,ft,0,D)
  R <- chol(st$Om + diag(1e-9,G)); Zr <- matrix(rnorm(ND*G),ND) %*% R
  pd <- mean(rowSums(Zr^2) >= st$dec_p)
  Zb <- st$Z; A <- solve(crossprod(Zb))
  se <- apply(Zr, 1, function(v){ zz <- crossprod(Zb, v); as.numeric(t(zz)%*%A%*%zz) })
  c(dec = pd, edge = mean(se >= st$edge_p))
}

# one replication -> four p-values
worker <- function(seed, panel, alt, gam, lam, do_mle) {
  set.seed(seed)
  if (panel == "A") {
    n <- 500; p <- 5; b0 <- c(0.5,-0.4,0.3,-0.3,0.2); D <- diag(c(0,rep(1,p)))
    X <- matrix(rnorm(n*p), n, p)
    dep <- if (alt == "quad") gam*(X[,1]^2 - 1) else gam*X[,1]*X[,2]
  } else {
    n <- 400; p <- 100; S <- 0.7^abs(outer(1:p,1:p,"-")); Ch <- chol(S)
    b0 <- rep(c(.35,-.30,.25,-.20,.15), length.out=p)*0.8887; D <- diag(c(0,rep(1,p)))
    X <- matrix(rnorm(n*p), n, p) %*% Ch
    dep <- gam*(X[,1]^2 - 1)
  }
  X1 <- cbind(1, X); y <- rbinom(n, 1, 1/(1+exp(-(drop(X %*% b0) + dep))))
  pp <- p_prepivot2(X1, y, lam, D)
  pm <- if (do_mle) p_mle2(X1, y, D) else c(dec=NA_real_, edge=NA_real_)
  c(p_pre_dec=pp["dec"], p_pre_edge=pp["edge"], p_mle_dec=pm["dec"], p_mle_edge=pm["edge"])
}

nc <- max(1, detectCores()-2)
cl <- makeCluster(nc); on.exit(stopCluster(cl))
clusterExport(cl, c("G","NB","ND","ridge_fit","pieces","p_prepivot2","p_mle2"))
clusterSetRNGStream(cl, 20260810)
cat("T2.5 power at B =", B, "on", nc, "workers\n\n")

cells <- c(
  lapply(c(0,0.25,0.5,0.75,1.0), function(g) list(p="A", alt="quad",  gam=g, lam=100, mle=TRUE)),
  lapply(c(0,0.25,0.5,0.75,1.0), function(g) list(p="A", alt="inter", gam=g, lam=100, mle=TRUE)),
  lapply(c(0,0.5,1.0,1.5),       function(g) list(p="B", alt="quad",  gam=g, lam=416, mle=FALSE)))

allp <- NULL; t00 <- Sys.time()
for (ce in cells) {
  t0 <- Sys.time()
  M <- simplify2array(parLapply(cl, seq_len(B), worker, panel=ce$p, alt=ce$alt,
                                gam=ce$gam, lam=ce$lam, do_mle=ce$mle))
  df <- data.frame(panel=ce$p, alt=ce$alt, gamma=ce$gam, lambda=ce$lam,
                   rep=seq_len(B), t(M))
  names(df)[6:9] <- c("p_pre_dec","p_pre_edge","p_mle_dec","p_mle_edge")
  allp <- rbind(allp, df)
  write.csv(allp, "T25_power_pvalues.csv", row.names=FALSE)   # incremental, crash-safe
  rd <- mean(df$p_pre_dec<=0.05); re <- mean(df$p_pre_edge<=0.05, na.rm=TRUE)
  md <- mean(df$p_mle_dec<=0.05); me <- mean(df$p_mle_edge<=0.05)
  cat(sprintf("%s %-5s gam=%.2f | prepivot dec %.4f edge %.4f | MLE dec %s edge %s | %.1f min (tot %.1f)\n",
      ce$p, ce$alt, ce$gam, rd, re,
      ifelse(ce$mle, sprintf("%.4f", md), "  --  "),
      ifelse(ce$mle, sprintf("%.4f", me), "  --  "),
      as.numeric(difftime(Sys.time(),t0,units="mins")),
      as.numeric(difftime(Sys.time(),t00,units="mins"))))
}
saveRDS(allp, "T25_power_results.rds")

## ---- summaries -------------------------------------------------------------
cat("\n===== NOMINAL vs SIZE-ADJUSTED POWER (alpha = 0.05) =====\n")
summ <- NULL
for (pn in unique(allp$panel)) for (al in unique(allp$alt[allp$panel==pn])) {
  sub <- allp[allp$panel==pn & allp$alt==al, ]
  if (!nrow(sub)) next
  nul <- sub[sub$gamma==0, ]
  # size-adjusted threshold: the 5th percentile of the NULL p-value distribution
  td <- quantile(nul$p_pre_dec,  0.05, type=1, na.rm=TRUE)
  te <- quantile(nul$p_pre_edge, 0.05, type=1, na.rm=TRUE)
  for (g in sort(unique(sub$gamma))) {
    s <- sub[sub$gamma==g, ]
    summ <- rbind(summ, data.frame(panel=pn, alt=al, gamma=g,
      dec=mean(s$p_pre_dec<=0.05), edge=mean(s$p_pre_edge<=0.05, na.rm=TRUE),
      dec_adj=mean(s$p_pre_dec<=td), edge_adj=mean(s$p_pre_edge<=te, na.rm=TRUE),
      mle_dec=mean(s$p_mle_dec<=0.05), mle_edge=mean(s$p_mle_edge<=0.05),
      se=round(sqrt(0.25/B),4)))
  }
  cat(sprintf("  [%s %s] size-adjusted thresholds: decile %.4f  EDGE %.4f\n", pn, al, td, te))
}
print(summ, row.names=FALSE, digits=4)
write.csv(summ, "T25_power_summary.csv", row.names=FALSE)

## the paired "cost of the correction" claim, panel A, decile basis, 8 alternatives
pa <- summ[summ$panel=="A" & summ$gamma>0, ]
d  <- pa$dec - pa$mle_dec
tt <- t.test(d)
cat(sprintf("\nPAIRED COST OF THE CORRECTION (panel A, decile, %d alternatives):\n", length(d)))
cat(sprintf("  mean difference %.4f, SE %.4f, t = %.2f on %d df, p = %.3f\n",
    mean(d), sd(d)/sqrt(length(d)), tt$statistic, tt$parameter, tt$p.value))
cat(sprintf("  (B=200 study gave -0.0100, SE 0.0050, t = -1.87, p = 0.10)\n"))
cat("\nwrote T25_power_pvalues.csv, T25_power_summary.csv, T25_power_results.rds\n")
