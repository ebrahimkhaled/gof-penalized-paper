# Does Omega_K >= Omega_MLE or <= ?  The paper claims the penalized covariance is SMALLER.
# The referee claims the identity  F^-1 - M^-1(F+2K)M^-1 = M^-1 K F^-1 K M^-1 >= 0
# makes it LARGER. Settle numerically on random well-conditioned F, K.
set.seed(11)
worst <- Inf; worstd <- Inf
for (rep in 1:500) {
  p <- sample(3:8, 1)
  A <- matrix(rnorm(p*p), p); F <- crossprod(A) + diag(p)*0.5      # F  > 0
  K <- diag(c(0, abs(rnorm(p-1))*runif(1,0,50)))                   # K >= 0, intercept unpenalised
  M <- F + K
  lhs <- solve(F) - solve(M) %*% (F + 2*K) %*% solve(M)
  rhs <- solve(M) %*% K %*% solve(F) %*% K %*% solve(M)
  stopifnot(max(abs(lhs - rhs)) < 1e-8)                            # the IDENTITY
  worst <- min(worst, min(eigen(lhs, symmetric=TRUE, only.values=TRUE)$values))
  # and the consequence for Omega, with a random U
  G <- 10; U <- matrix(rnorm(G*p), G, p)
  OmK  <- diag(G) - U %*% solve(M) %*% (F + 2*K) %*% solve(M) %*% t(U)
  OmML <- diag(G) - U %*% solve(F) %*% t(U)
  worstd <- min(worstd, min(eigen(OmK - OmML, symmetric=TRUE, only.values=TRUE)$values))
}
cat(sprintf("identity F^-1 - M^-1(F+2K)M^-1 = M^-1 K F^-1 K M^-1 : HOLDS in all 500 draws\n"))
cat(sprintf("min eigenvalue of  F^-1 - M^-1(F+2K)M^-1  over 500 draws : %.3e  (>=0 means M^-1(F+2K)M^-1 <= F^-1)\n", worst))
cat(sprintf("min eigenvalue of  Omega_K - Omega_MLE     over 500 draws : %.3e\n", worstd))
cat(sprintf("\nVERDICT: Omega_K %s Omega_MLE\n", ifelse(worstd >= -1e-9, ">= (penalized covariance is LARGER)", "<= ")))
