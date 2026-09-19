# Figure 1 for the paper — THE PENALTY PATH. EK house style.
#
# Replaces Table 1. Its job is to render, as a SHAPE, the fact a table cannot show:
# ridge shrinkage improves the fitted probabilities and destroys the goodness-of-fit test
# over the SAME stretch of the penalty path, and cross-validation lands beyond the point
# at which the test has already failed.
#
# Every number is registered. Manifest section given at each vector.
## Run from this directory:  cd paper_seriesC/R && Rscript fig1_penaltypath_ek.R
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(scales)
                                library(ggrepel)})
source("_ek_theme.R")

## ---- registered data ------------------------------------------------------------
## manifest 2: naive size under H0, B = 500, MC SE <= 0.022
lgA  <- c(0, 0.02, 0.05, 0.10, 0.20, 0.274)
hlA  <- c(0.068, 0.068, 0.118, 0.362, 0.916, 0.988)
edA  <- c(0.034, 0.046, 0.170, 0.888, 1.000, 1.000)
lgB  <- c(0, 0.01, 0.05, 0.124, 0.513, 1.5)
hlB  <- c(0.178, 0.432, 0.990, 0.998, 1.000, 1.000)   # 0.124 column measured by T38 (was typed 0.994)
edB  <- c(0.284, 0.994, 1.000, 1.000, 1.000, 1.000)
maeB <- c(0.1859, 0.1469, 0.1330, 0.1388, 0.1606, 0.1840)   # manifest 2, DGP B; 0.124 from T38

## manifest 2: cross-validated penalties (medians; IQR registered for A only)
cvA_min <- 0.0197; cvA_1se <- 0.2331; cvA_lo <- 0.185; cvA_hi <- 0.295
cvB_min <- 0.1245; cvB_1se <- 1.0404

## manifest 8f + 8i: the PREPIVOTED corrected test, decile basis. A separate experiment
## at a different replication count -- plotted as ISOLATED points, never joined by a line.
## Design B lambda=50 uses the POOLED N=4000 value (8i), NOT the superseded 8f 0.0630.
corA_x <- c(100, 137, 200)/500;  corA_y <- c(0.0435, 0.0475, 0.0465)
corA_s <- c(0.0046, 0.0048, 0.0047)
corB_x <- c(50, 416, 1000)/400;  corB_y <- c(0.0595, 0.0710, 0.0520)
corB_s <- c(0.0037, 0.0081, 0.0070)

GREY <- "grey45"; GRN <- "#009E73"; BLU <- "#0072B2"; ORG <- "#D55E00"

## pseudo-log so the lambda_g = 0 column (the maximum-likelihood fit, which carries the
## Sur-Candes point) is not silently dropped by a plain log axis.
psl <- pseudo_log_trans(sigma = 0.01, base = 10)

panel_size <- function(lg, hl, ed, brks, ttl, sub, cx, cy, cs, vmin, v1se,
                       band = NULL, ylab = NULL, note = NULL, note_x = NULL) {
  d <- rbind(data.frame(lg, r = hl, basis = "decile"),
             data.frame(lg, r = ed, basis = "EDGE"))
  d <- d[order(d$basis, d$lg), ]
  p <- ggplot(d, aes(lg, r, colour = basis, linetype = basis, shape = basis))
  if (!is.null(band))
    p <- p + annotate("rect", xmin = band[1], xmax = band[2], ymin = -Inf, ymax = Inf,
                      fill = "grey70", alpha = 0.12)
  p +
    ek_nominal() +
    geom_vline(xintercept = vmin, linetype = "13", colour = "grey55", linewidth = 0.35) +
    geom_vline(xintercept = v1se, linetype = "solid", colour = "grey40", linewidth = 0.4) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.9) +
    ## the corrected test: isolated points + 2 s.e., deliberately NOT joined
    annotate("errorbar", x = cx, ymin = pmax(cy - 2*cs, 0), ymax = cy + 2*cs,
             width = 0.06, colour = BLU, linewidth = 0.35) +
    annotate("point", x = cx, y = cy, shape = 22, size = 2.1, colour = BLU, fill = "white") +
    scale_colour_manual(values = c(decile = GREY, EDGE = GRN), name = NULL) +
    scale_linetype_manual(values = c(decile = "solid", EDGE = "22"), name = NULL) +
    scale_shape_manual(values = c(decile = 19, EDGE = 17), name = NULL) +
    scale_x_continuous(trans = psl, breaks = brks, labels = brks) +
    scale_y_continuous(limits = c(0, 1.06), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
    labs(title = ttl, x = expression(lambda[g]), y = ylab) +
    theme_ek(base_size = 9) +
    theme(legend.position = "none",
          plot.title.position = "panel",
          plot.title = element_text(hjust = 0, face = "bold", size = rel(0.95)),
          plot.subtitle = element_blank())
}

pA <- panel_size(lgA, hlA, edA, c(0, 0.05, 0.1, 0.274), "A.  Design A  (p = 5)",
                 "fails before CV stops",
                 corA_x, corA_y, corA_s, cvA_min, cvA_1se,
                 band = c(cvA_lo, cvA_hi), ylab = "rejection rate of a CORRECT model") +
      annotate("text", x = 0.023, y = 0.92, hjust = 0, size = 2.2, colour = "grey35",
           lineheight = 0.95, label = "cross-validation\nlands in\nthe band") +
  annotate("text", x = 0.30, y = 0.135, hjust = 0.5, size = 2.2, colour = BLU,
           label = "corrected")

pB <- panel_size(lgB, hlB, edB, c(0, 0.124, 0.513, 1.5), "B.  Design B  (p = 100)",
                 "invalid even unpenalized",
                 corB_x, corB_y, corB_s, cvB_min, cvB_1se) +
  annotate("text", x = 0.055, y = 0.30, hjust = 0, size = 2.2, colour = "grey25",
           lineheight = 0.95, label = "unpenalized MLE:\nalready 0.178") +
  annotate("text", x = 0.30, y = 0.16, hjust = 0, size = 2.2, colour = BLU, label = "corrected")

## ---- panel C: the trade ----------------------------------------------------------
dC <- data.frame(mae = maeB, rej = hlB, lg = lgB)
pC <- ggplot(dC, aes(mae, rej)) +
  ek_nominal() +
  geom_path(arrow = arrow(length = unit(0.05, "in"), type = "closed"),
            colour = ORG, linewidth = 0.6) +
  geom_point(size = 2, colour = ORG) +
  annotate("point", x = 0.1330, y = 0.990, shape = 23, size = 3.1,
           colour = ORG, fill = ORG) +
  geom_text_repel(aes(label = sprintf("%.3g", lg)), size = 2.1, colour = "grey30",
                  min.segment.length = 0.1, segment.size = 0.2, segment.colour = "grey60",
                  box.padding = 0.28, seed = 7) +
  annotate("text", x = 0.152, y = 0.66, hjust = 0.5, size = 2.15, colour = "grey20",
           lineheight = 0.95,
           label = "the best probabilities this\nmodel ever produces,\nrejected 99% of the time") +
  scale_x_reverse(limits = c(0.196, 0.122)) +
  scale_y_continuous(limits = c(0, 1.06), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(title = "C.  The trade, design B",
       x = "MAE  (better →)",
       y = NULL) +
  theme_ek(base_size = 9) +
  theme(plot.title.position = "panel",
        plot.title = element_text(hjust = 0, face = "bold", size = rel(0.95)),
        plot.subtitle = element_blank())

pl <- pA | pB | pC
BJ_W <- 423.72342/72.27
ggsave("../Fig/fig1_penaltypath.pdf", pl, width = BJ_W, height = BJ_W*0.40, device = cairo_pdf)
ggsave("../Fig/fig1_penaltypath.png", pl, width = BJ_W, height = BJ_W*0.40, dpi = EK_DPI)

## ---- guards: every claim the caption makes -----------------------------------------
stopifnot(lgB[which.min(maeB)] == 0.05,              # the MAE minimum is at lambda_g = 0.05
          min(maeB) == 0.1330, maeB[1] == 0.1859,    # and it is a 28% improvement
          hlB[which.min(maeB)] == 0.990)             # where the test rejects 99% of the time
stopifnot(all(hlA[lgA >= 0.20] >= 0.916))            # design A has failed by the CV band
stopifnot(hlB[1] == 0.178)                           # the unpenalized column is already invalid
stopifnot(cvA_1se > 0.05)                            # CV lands beyond the failure point
cat(sprintf("MAE falls %.1f%% (%.4f -> %.4f) while HL rises %.3f -> %.3f\n",
            100*(maeB[1]-min(maeB))/maeB[1], maeB[1], min(maeB), hlB[1], hlB[which.min(maeB)]))
cat("wrote ../Fig/fig1_penaltypath.pdf/.png\n")
