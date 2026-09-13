# Power curves re-parameterised by the VISIBLE EFFECT SIZE (task T2.2 -> figure).
# Replaces the gamma-axis curves, which are non-monotone (manifest 8d) and compare
# the two designs on an axis that means different things in each.
#
# Uses the df-CORRECTED tau (manifest E12: the raw plug-in overstates by up to 90%
# because poly(eta.hat,3) spends 3 fitted degrees of freedom).
# All power cells: known-null Monte Carlo, B1 = 1000, MC SE <= 0.0158.
d <- data.frame(
  panel = rep(c("A","A","B","B"), each = 5),
  alt   = rep(c("quad","cubic","quad","cubic"), each = 5),
  gamma = c(0.340,0.680,1.020,1.531,2.126,  0.123,0.247,0.370,0.555,0.770,
            1.117,2.235,3.352,5.028,6.000,  0.525,1.050,1.575,2.362,3.281),
  tau   = c(0.436,0.777,1.130,1.775,2.828,  0.280,0.585,0.933,1.479,2.247,
            0.195,0.777,2.088,5.370,7.422,  0.412,0.835,1.191,1.297,1.673),
  nu    = c(0.172,0.284,0.384,0.536,0.747,  0.111,0.222,0.333,0.480,0.647,
            0.161,0.413,0.789,1.552,1.983,  0.188,0.337,0.417,0.408,0.481),
  dec   = c(0.101,0.217,0.294,0.374,0.458,  0.057,0.062,0.108,0.215,0.372,
            0.077,0.102,0.166,0.180,0.192,  0.062,0.096,0.101,0.146,0.148),
  edge  = c(0.179,0.368,0.474,0.516,0.568,  0.049,0.091,0.166,0.362,0.583,
            0.071,0.078,0.131,0.155,0.149,  0.049,0.056,0.058,0.081,0.076))
d$nsqrtn <- d$nu * ifelse(d$panel == "A", sqrt(500), sqrt(400))

acc <- "#1a5276"; wrn <- "#b03a2e"; grn <- "#1e8449"; gry <- "#7f8c8d"
png("power_tau.png", width = 1620, height = 560, res = 130)
par(mfrow = c(1,3), mar = c(4.5,4.5,3.6,1), cex.lab = 1.02)

## (a) THE SUCCESS: index-aligned departure, decile basis, matched on nu*sqrt(n)
a <- d[d$alt=="cubic" & d$panel=="A",]; b <- d[d$alt=="cubic" & d$panel=="B",]
plot(NA, xlim=c(0,15), ylim=c(0,0.42), xlab=expression(paste(nu, sqrt(n))),
     ylab="power (rejection rate)",
     main="A. index-aligned departure\ncurves collapse (gap 0.016)")
abline(h=0.05, lty=2, col=gry)
lines(a$nsqrtn, a$dec, col=acc, lwd=2.6); points(a$nsqrtn, a$dec, col=acc, pch=19, cex=1.2)
lines(b$nsqrtn, b$dec, col=wrn, lwd=2.6, lty=2); points(b$nsqrtn, b$dec, col=wrn, pch=17, cex=1.2)
legend("topleft", bty="n", cex=0.82, legend=c("design A  (p=5)","design B  (p=100)"),
       col=c(acc,wrn), lwd=2.6, lty=c(1,2), pch=c(19,17))
text(8, 0.40, "decile basis", cex=0.85, col=gry)

## (b) THE CEILING: coordinate-wise departure does NOT collapse
a <- d[d$alt=="quad" & d$panel=="A",]; b <- d[d$alt=="quad" & d$panel=="B",]
plot(NA, xlim=c(0,7.6), ylim=c(0,0.62), xlab=expression(paste("visible effect size  ", tau)),
     ylab="", main="B. coordinate-wise departure\ndesign B ceilings at 0.19")
abline(h=0.05, lty=2, col=gry)
lines(a$tau, a$dec, col=acc, lwd=2.6); points(a$tau, a$dec, col=acc, pch=19, cex=1.2)
lines(b$tau, b$dec, col=wrn, lwd=2.6, lty=2); points(b$tau, b$dec, col=wrn, pch=17, cex=1.2)
abline(h=0.192, lty=3, col=wrn)
text(4.6, 0.225, "ceiling 0.19", col=wrn, cex=0.8)
arrows(2.83, 0.458, 2.83, 0.55, length=0.07, col=acc); text(2.83, 0.585, "still rising", col=acc, cex=0.8)
legend("topleft", bty="n", cex=0.82, legend=c("design A  (p=5)","design B  (p=100)"),
       col=c(acc,wrn), lwd=2.6, lty=c(1,2), pch=c(19,17))

## (c) THE REVERSAL: EDGE beats decile in A, loses in B
a <- d[d$alt=="cubic" & d$panel=="A",]; b <- d[d$alt=="cubic" & d$panel=="B",]
plot(NA, xlim=c(0,2.4), ylim=c(0,0.62), xlab=expression(paste("visible effect size  ", tau)),
     ylab="", main="C. the EDGE advantage reverses\nwith dimension")
abline(h=0.05, lty=2, col=gry)
lines(a$tau, a$edge, col=grn, lwd=2.6); points(a$tau, a$edge, col=grn, pch=19, cex=1.2)
lines(a$tau, a$dec,  col=acc, lwd=2.2, lty=3); points(a$tau, a$dec, col=acc, pch=1, cex=1.1)
lines(b$tau, b$edge, col=grn, lwd=2.6, lty=2); points(b$tau, b$edge, col=grn, pch=15, cex=1.1)
lines(b$tau, b$dec,  col=wrn, lwd=2.2, lty=3); points(b$tau, b$dec, col=wrn, pch=0, cex=1.1)
legend("topleft", bty="n", cex=0.72,
  legend=c("A: EDGE","A: decile","B: EDGE","B: decile"),
  col=c(grn,acc,grn,wrn), lwd=c(2.6,2.2,2.6,2.2), lty=c(1,3,2,3), pch=c(19,1,15,0))
text(1.6, 0.13, "EDGE below decile\nthroughout design B", cex=0.72, col=wrn)
dev.off()
cat("wrote power_tau.png\n")

## the collapse arithmetic, printed for the caption
ia <- approx(d$nsqrtn[d$alt=="cubic"&d$panel=="A"], d$dec[d$alt=="cubic"&d$panel=="A"],
             xout=d$nsqrtn[d$alt=="cubic"&d$panel=="B"], rule=1)$y
gap <- abs(ia - d$dec[d$alt=="cubic"&d$panel=="B"])
cat("cubic/decile matched on nu*sqrt(n): mean |gap| =", round(mean(gap,na.rm=TRUE),4),
    " (MC SE of a single difference ~0.023)\n")
saveRDS(d,"tau_curve_data.rds")
