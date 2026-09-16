# Every number in Table S8 (size with a cross-validated penalty) and in the like-for-like
# comparison of Section 5.3, computed from the per-replication p-values the studies archived.
# Run from R/. Reads ../data/T31_cv_lambda_pvalues.csv, ../data/T32_glaucoma_design_pvalues.csv,
# ../data/T37_glaucoma_T1_pvalues.csv and ../data/T30_keycells_pvalues.csv.
# Exact binomial intervals are printed with their bounds rounded outwards, as in the table.
# Ends in stopifnot() guards on the values the paper prints.

rate  <- function(p) mean(p < 0.05)
ci    <- function(p) { k <- sum(p < 0.05); binom.test(k, length(p))$conf.int }
outw  <- function(x) c(floor(x[1]*1000)/1000, ceiling(x[2]*1000)/1000)
lamq  <- function(l) round(quantile(l, c(.5, .25, .75), na.rm = TRUE), 1)
row   <- function(lab, lam, pd, pe) {
  q <- lamq(lam); a <- outw(ci(pd)); b <- outw(ci(pe))
  cat(sprintf("%-34s %5d  %6.1f (%6.1f, %6.1f)  %.4f [%.3f, %.3f]  %.4f [%.3f, %.3f]\n",
              lab, length(pd), q[1], q[2], q[3], rate(pd), a[1], a[2], rate(pe), b[1], b[2]))
  invisible(list(n = length(pd), q = unname(q), dec = rate(pd), edge = rate(pe),
                 ci_dec = a, ci_edge = b))
}

t31 <- read.csv("../data/T31_cv_lambda_pvalues.csv")
t32 <- read.csv("../data/T32_glaucoma_design_pvalues.csv")
t37 <- read.csv("../data/T37_glaucoma_T1_pvalues.csv")
t30 <- read.csv("../data/T30_keycells_pvalues.csv")

cat("Table S8: size of the prepivoted corrected test, penalty chosen by 10-fold CV in every replicate\n")
cat(sprintf("%-34s %5s  %-22s  %-22s  %-22s\n", "design / rule or truth", "reps",
            "lambda median (IQR)", "decile [95%]", "EDGE [95%]"))
A <- t31[t31$panel == "A", ]; B <- t31[t31$panel == "B", ]
a1 <- row("A, 1se",  A$lam_1se, A$dec_1se, A$edge_1se)
a0 <- row("A, min",  A$lam_min, A$dec_min, A$edge_min)
b1 <- row("B, 1se",  B$lam_1se, B$dec_1se, B$edge_1se)
b0 <- row("B, min",  B$lam_min, B$dec_min, B$edge_min)
T1a <- t32[t32$truth == "T1", ]; T2 <- t32[t32$truth == "T2", ]
cat(sprintf("  glaucoma truth pi(beta.hat) at lambda.min, two independent blocks:\n"))
cat(sprintf("    seeds 400001+: decile %.3f, EDGE %.3f\n", rate(T1a$p_dec), rate(T1a$p_edge)))
cat(sprintf("    seeds 420001+: decile %.3f, EDGE %.3f\n", rate(t37$p_dec), rate(t37$p_edge)))
g1 <- row("glaucoma, pi(beta.hat), pooled", c(T1a$lam, t37$lam),
          c(T1a$p_dec, t37$p_dec), c(T1a$p_edge, t37$p_edge))
g2 <- row("glaucoma, pi(beta.tilde)", T2$lam, T2$p_dec, T2$p_edge)

cat("\nSection 5.3: cross-validated rates against the fixed-penalty cells at the same bootstrap size\n")
f416 <- t30[t30$panel == "B" & t30$lambda == 416, ]
f50  <- t30[t30$panel == "B" & t30$lambda == 50, ]
cmp <- function(lab, cv, fx) {
  d  <- rate(cv) - rate(fx)
  se <- sqrt(rate(cv)*(1 - rate(cv))/length(cv) + rate(fx)*(1 - rate(fx))/length(fx))
  cat(sprintf("  %-22s CV %.4f  fixed %.4f  difference %+.4f  = %.2f SE of the difference\n",
              lab, rate(cv), rate(fx), d, d/se))
  d/se
}
z416d <- cmp("lambda 416, decile", B$dec_1se,  f416$dec499)
z416e <- cmp("lambda 416, EDGE",   B$edge_1se, f416$edge499)
z50d  <- cmp("lambda 50, decile",  B$dec_min,  f50$dec499)
z50e  <- cmp("lambda 50, EDGE",    B$edge_min, f50$edge499)

# ---- guards: the values the paper and supplement print ----------------------------------------
near <- function(x, y, tol = 5e-5) abs(x - y) <= tol + 1e-12   # printed to 4 dp: within half a unit
stopifnot(
  near(a1$dec, .0580), near(a1$edge, .0595), near(a0$dec, .0520), near(a0$edge, .0560),
  near(b1$dec, .0710), near(b1$edge, .0530), near(b0$dec, .0780), near(b0$edge, .0510),
  all(near(a1$q, c(117.6, 93.9, 147.8), .051)), all(near(b1$q, c(427.0, 269.5, 616.5), .051)),
  all(near(a0$q, c(10.0, 6.8, 13.0), .051)),   all(near(b0$q, c(51.8, 30.9, 87.9), .051)),
  near(g1$dec, .0425), near(g1$edge, .0760), all(near(g1$ci_edge, c(.064, .089))),
  all(near(g1$q, c(202.6, 148.0, 272.5), .051)),
  near(g2$dec, .038), near(g2$edge, .048), all(near(g2$q, c(266.7, 171.7, 387.0), .051)),
  near(rate(f416$dec499), .0680), near(rate(f416$edge499), .0560),
  near(rate(f50$dec499), .0592),  near(rate(f50$edge499), .0408),
  abs(z50d - 2) < 0.1, abs(z50e - 4/3) < 0.1
)
cat("\nall printed values reproduced\n")
