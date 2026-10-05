# Figure 5 for the paper -- THE VALIDITY MAP. EK house style.
#
# One question per panel: in this design, which tests reject correct models at the rate they claim, as p/n
# grows? The shaded band is the pre-declared range [0.03, 0.08]. The rate axis is logarithmic because the
# interesting failures run in both directions: a test that never rejects (GRPtests with five splits) is as
# useless as one that always does (HL, and PLStests once the model has an intercept). The classical test on
# the maximum likelihood fit is drawn only where that estimate exists, with the fraction of datasets in which
# it does printed beside each point. BAGofT, at about an hour per test, was run only in panel B at
# p/n = 0.05 and 0.25 (100 replicates per point).
#
# Data: gof-penalized-paper/data/T40_validity_map.csv (T40_validity_map_data.py; T40, T41, T42 and the
# Amendment 5/6 second blocks pooled, 1000-2000 replicates per point).
## Run from this directory:  cd paper_seriesC/R && Rscript fig5_validity_map_ek.R
suppressPackageStartupMessages({library(ggplot2)})
source("_ek_theme.R")

# runs from R/ of the archive (data in ../data) or of the paper folder (data in the archive next to it)
SRC <- c("../data/T40_validity_map.csv", "../../gof-penalized-paper/data/T40_validity_map.csv")
d <- read.csv(SRC[file.exists(SRC)][1])
d$design <- factor(d$design, levels = c("Dense, correlation 0.4", "Dense, correlation 0.7",
                                        "Dense, correlation 0.8", "Sparse, correlation 0.7",
                                        "Dense, 25% events"),
                   labels = c("A. Dense signal, correlation 0.4", "B. Dense signal, correlation 0.7",
                              "C. Dense signal, correlation 0.8", "D. Sparse signal, correlation 0.7",
                              "E. Dense signal, 25% events"))
LEV <- c("SC.EDGE", "SC.HL", "HL", "HL on MLE", "GRPtests (5 splits)", "GRPtests (1 split)", "PLStests", "BAGofT")
d$method <- factor(d$method, levels = LEV)
d <- d[!is.na(d$rate), ]
FLOOR <- 0.001                                   # a rate of 0 is drawn on the floor of the log axis
d$y <- pmax(d$rate, FLOOR)

pal <- c("SC.EDGE" = "#009E73", "SC.HL" = "#7FC9B0", "HL" = "#56B4E9", "HL on MLE" = "#000000",
         "GRPtests (5 splits)" = "#E69F00", "GRPtests (1 split)" = "#D55E00", "PLStests" = "#CC79A7",
         "BAGofT" = "#8C6D31")
lty <- c("SC.EDGE" = "solid", "SC.HL" = "solid", "HL" = "22", "HL on MLE" = "13",
         "GRPtests (5 splits)" = "42", "GRPtests (1 split)" = "42", "PLStests" = "42", "BAGofT" = "solid")
shp <- c("SC.EDGE" = 19, "SC.HL" = 1, "HL" = 2, "HL on MLE" = 15,
         "GRPtests (5 splits)" = 0, "GRPtests (1 split)" = 5, "PLStests" = 6, "BAGofT" = 8)
lwd <- c("SC.EDGE" = 0.9, "SC.HL" = 0.6, "HL" = 0.45, "HL on MLE" = 0.45,
         "GRPtests (5 splits)" = 0.45, "GRPtests (1 split)" = 0.45, "PLStests" = 0.45, "BAGofT" = 0.45)

mle_lab <- d[d$method == "HL on MLE" & !is.na(d$mle_exists) & d$mle_exists < 0.995, ]
mle_lab$lab <- sprintf("%.0f%%", 100 * mle_lab$mle_exists)

pl <- ggplot(d, aes(kappa, y, colour = method, linetype = method, shape = method)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.03, ymax = 0.08, fill = "grey70", alpha = 0.25) +
  geom_hline(yintercept = 0.05, linetype = "13", colour = "grey45", linewidth = 0.35) +
  geom_line(aes(linewidth = method)) +
  geom_point(size = 1.5, stroke = 0.6) +
  geom_text(data = mle_lab, aes(label = lab), colour = "grey20", size = 2.3, family = "serif",
            vjust = -0.9, show.legend = FALSE) +
  scale_colour_manual(values = pal, breaks = LEV, name = NULL) +
  scale_linetype_manual(values = lty, breaks = LEV, name = NULL) +
  scale_shape_manual(values = shp, breaks = LEV, name = NULL) +
  scale_linewidth_manual(values = lwd, breaks = LEV, guide = "none") +
  scale_y_log10(limits = c(FLOOR, 1.25), breaks = c(0.001, 0.01, 0.05, 0.1, 0.3, 1),
                labels = c("0.001", "0.01", "0.05", "0.1", "0.3", "1")) +
  scale_x_continuous(breaks = c(0.05, 0.10, 0.25, 0.35, 0.45), labels = c(".05", ".10", ".25", ".35", ".45")) +
  facet_wrap(~ design, ncol = 3) +
  labs(x = expression(italic(p) / italic(n)), y = "rejection rate of correct models (log scale)") +
  theme_ek(base_size = 9) +
  theme(legend.position = c(0.84, 0.22), legend.key.width = unit(1.6, "lines"),
        legend.text = element_text(size = rel(0.9)), legend.spacing.y = unit(0, "pt"),
        strip.text = element_text(hjust = 0, size = rel(0.92)), panel.spacing = unit(0.9, "lines")) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(linewidth = 0.6)))

## GEOMETRY: the text width of the other figures (6.181in), so includegraphics[width=textwidth] applies
## no rescaling and the 9pt base type prints at 9pt.
OUP_W <- 6.181
for (dir in c("../Fig", "../bimj/Fig")) if (dir.exists(dir))
  ggsave(file.path(dir, "fig5_validity_map.pdf"), pl, width = OUP_W, height = OUP_W * 0.66, device = cairo_pdf)
ggsave("../Fig/fig5_validity_map.png", pl, width = OUP_W, height = OUP_W * 0.66, dpi = 220)

## the claims the figure carries in the text (Section 5.4)
r <- function(m, dsg, k) d$rate[d$method == m & grepl(dsg, d$design) & abs(d$kappa - k) < 1e-9]
for (dsg in levels(d$design)) for (k in c(0.05, 0.10)) for (m in c("SC.HL", "SC.EDGE")) {
  v <- r(m, substr(dsg, 4, 99), k); if (length(v)) stopifnot(v >= 0.03, v <= 0.08) }   # valid at p/n <= 0.10
for (dsg in c("correlation 0.4", "Dense signal, correlation 0.7", "correlation 0.8"))
  stopifnot(r("SC.EDGE", dsg, 0.25) <= 0.08, r("SC.HL", dsg, 0.25) <= 0.08)          # and at 0.25, dense balanced
stopifnot(r("SC.HL", "Sparse", 0.25) > 0.08,                                             # but not sparse
          all(d$rate[d$method == "PLStests" & grepl("25%", d$design)] == 1),             # PLStests, intercept
          !any(d$method == "HL on MLE" & abs(d$kappa - 0.45) < 1e-9 & !is.na(d$rate)))   # no MLE at 0.45
cat("figure 5 claims verified\n")
