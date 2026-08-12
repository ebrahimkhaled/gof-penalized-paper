# Which "uncorrected test"? Two different things have been called that.
#
#   HL   = textbook Hosmer-Lemeshow: S = sum (O-E)^2 / {m pbar(1-pbar)}  ~ chi^2_{G-2}
#          This is what section 1 of the paper describes and what practice applies.
#   OmML = the same grouped residuals, standardised by Vg = sum w_i, with NO mu-hat
#          correction, referred by Monte Carlo to Omega_MLE = I - U F^-1 U'.
#          This is what glaucoma_deep.R computed and what the manifest's
#          "naive 0.0000-0.0060" refers to.
#
# They are not the same test and they do not agree. This script reports both for all
# twelve lambda x G configurations so the paper can state exactly what each one does.
suppressPackageStartupMessages({library(glmnet); library(TH.data)})
data("GlaucomaM", package = "TH.data")
X <- scale(as.matrix(GlaucomaM[, setdiff(names(GlaucomaM), "Class")]))
X <- X[, apply(X, 2, function(c) all(is.finite(c)) & sd(c) > 0), drop = FALSE]
y <- as.numeric(GlaucomaM$Class == "glaucoma"); n <- nrow(X); p <- ncol(X)
X1 <- cbind(1, X); D <- diag(c(0, rep(1, p)))
set.seed(1); cv <- cv.glmnet(X, y, family = "binomial", alpha = 0)
lmin <- cv$lambda.min * n; l1se <- cv$lambda.1se * n

ridge_fit <- function(X1,y,lam,D,tol=1e-9,maxit=100){
  beta <- rep(0,ncol(X1))
  for(it in 1:maxit){ eta<-drop(X1%*%beta); pi<-1/(1+exp(-eta))
    w<-pmax(pi*(1-pi),1e-8); z<-eta+(y-pi)/w
    bn<-drop(solve(crossprod(X1,w*X1)+lam*D, crossprod(X1,w*z)))
    if(max(abs(bn-beta))<tol){beta<-bn;break}; beta<-bn }
  list(beta=beta, pi=1/(1+exp(-drop(X1%*%beta))))
}

both <- function(lam, G, nd = 200000) {
  ft <- ridge_fit(X1, y, lam, D)
  pi <- pmin(pmax(ft$pi,1e-8),1-1e-8); w <- pmax(pi*(1-pi),1e-10)
  g <- pmin(ceiling(rank(pi, ties.method="first")/(n/G)), G)
  idx <- split(seq_len(n), g)
  O <- vapply(idx, function(I) sum(y[I]), 0.0)
  E <- vapply(idx, function(I) sum(pi[I]), 0.0)
  m <- vapply(idx, function(I) length(I), 0.0)
  pb <- E/m
  # (1) textbook Hosmer-Lemeshow
  S_hl <- sum((O-E)^2/(m*pb*(1-pb)));  p_hl <- pchisq(S_hl, G-2, lower.tail=FALSE)
  # (2) uncorrected statistic vs Omega_MLE, Monte Carlo reference
  Vg <- vapply(idx, function(I) sum(w[I]), 0.0)
  r  <- (O-E)/sqrt(Vg)
  U  <- t(vapply(idx, function(I) colSums(w[I]*X1[I,,drop=FALSE]), numeric(ncol(X1))))/sqrt(Vg)
  Fm <- crossprod(X1, w*X1)
  Om <- diag(G) - U %*% solve(Fm, t(U))
  R  <- chol(Om + diag(1e-9, G))
  Z  <- matrix(rnorm(nd*G), nd) %*% R
  S_om <- sum(r^2); p_om <- mean(rowSums(Z^2) >= S_om)
  c(S_hl=S_hl, p_hl=p_hl, S_om=S_om, p_om=p_om, dfOm=sum(diag(Om)))
}

set.seed(20260811)
lams <- c(lambda.min=lmin, half.1se=l1se/2, lambda.1se=l1se, twice.1se=2*l1se)
res <- NULL
for (nm in names(lams)) for (G in c(5,10,20)) {
  b <- both(lams[[nm]], G)
  res <- rbind(res, data.frame(lambda=nm, lam=round(lams[[nm]],1), G=G,
    HL_S=round(b["S_hl"],2), HL_p=b["p_hl"],
    Om_S=round(b["S_om"],2), Om_p=b["p_om"], Om_df=round(b["dfOm"],2)))
}
rownames(res) <- NULL
fmt <- function(x) ifelse(x < 1e-4, sprintf("%.1e", x), sprintf("%.4f", x))
res$HL_p_fmt <- fmt(res$HL_p); res$Om_p_fmt <- fmt(res$Om_p)
print(res[,c("lambda","G","HL_S","HL_p_fmt","Om_S","Om_p_fmt","Om_df")], row.names=FALSE)
cat(sprintf("\nHosmer-Lemeshow (chi2_{G-2}) range : [%.2e, %.4f]   cells with p>=0.05: %d/12\n",
            min(res$HL_p), max(res$HL_p), sum(res$HL_p >= 0.05)))
cat(sprintf("Uncorrected vs Omega_MLE     range : [%.2e, %.4f]   cells with p>=0.05: %d/12\n",
            min(res$Om_p), max(res$Om_p), sum(res$Om_p >= 0.05)))
cat(sprintf("\nHEADLINE (lambda.1se, G=10): HL p = %.2e ; Omega_MLE p = %.2e\n",
   res$HL_p[res$lambda=="lambda.1se" & res$G==10],
   res$Om_p[res$lambda=="lambda.1se" & res$G==10]))
write.csv(res, "naive_glaucoma_both.csv", row.names=FALSE)
cat("wrote naive_glaucoma_both.csv\n")
