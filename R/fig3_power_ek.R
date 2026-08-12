# Figure 3 for the Series C paper — the three power facts, EK house style.
# Replaces the base-R three-panel plot (plot_tau.R), which looked like default output.
#
# Data: manifest 8g / T2.2. df-CORRECTED tau (manifest E12 — the raw plug-in overstates
# by up to 90% because poly(eta.hat,3) spends 3 fitted degrees of freedom).
# All power cells: known-null Monte Carlo, B1 = 1000, MC SE <= 0.0158.
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
## Run from this directory:  cd paper_seriesC/R && Rscript <script>
## _ek_theme.R ships WITH the archive, so the figures rebuild without the author tree.
source("_ek_theme.R")

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
d$design <- ifelse(d$panel == "A", "design A  (p = 5)", "design B  (p = 100)")

A_COL <- "#0072B2"; B_COL <- "#D55E00"; E_COL <- "#009E73"
des_sc <- list(
  scale_colour_manual(values = c("design A  (p = 5)" = A_COL, "design B  (p = 100)" = B_COL),
                      name = NULL),
  scale_shape_manual(values  = c("design A  (p = 5)" = 19, "design B  (p = 100)" = 17),
                     name = NULL),
  scale_linetype_manual(values = c("design A  (p = 5)" = "solid", "design B  (p = 100)" = "22"),
                        name = NULL))

base <- function(p) p + ek_nominal() +
  geom_line(linewidth = 0.7) + geom_point(size = 1.9) +
  des_sc + ylab(NULL) +
  scale_y_continuous(limits = c(0, 0.63), breaks = seq(0, 0.6, 0.15)) +
  theme_ek(base_size = 9) +
  theme(legend.position = "none", legend.justification = c(0, 1),
        legend.key.width = unit(20, "pt"), legend.text = element_text(size = rel(0.8)),
        legend.background = element_blank(), legend.spacing.y = unit(0, "pt"),
        plot.title = element_text(hjust = 0, face = "bold", size = rel(0.95)),
        plot.subtitle = element_blank())

## B1 for every plotted cell (known-null Monte Carlo); used for the error bars.
B1 <- 1000
se_of <- function(p) sqrt(p * (1 - p) / B1)

## (A) index-aligned departure: the curves COLLAPSE on nu*sqrt(n)
## SORT within design: nu is NOT monotone in gamma (design B has nu = ... 0.417, 0.408
## ...), so drawing in data order makes the line double back on itself.
a <- d[d$alt == "cubic", ]
a <- a[order(a$design, a$nsqrtn), ]
pA <- base(ggplot(a, aes(nsqrtn, dec, colour = design, shape = design, linetype = design))) +
  labs(title = "A.  Index-aligned departure",
       x = expression(nu*sqrt(n))) +
  ylab("power (rejection rate)") +
  geom_errorbar(aes(ymin = pmax(dec - se_of(dec), 0), ymax = dec + se_of(dec)),
                width = 0.28, linewidth = 0.3, show.legend = FALSE) +
  annotate("text", x = 12.5, y = 0.11, label = "mean gap 0.017\n(MC s.e. 0.023)",
           size = 2.3, colour = "grey35", lineheight = 0.95)

## (B) coordinate-wise departure: design B CEILINGS
b <- d[d$alt == "quad", ]
b <- b[order(b$design, b$tau), ]
pB <- base(ggplot(b, aes(tau, dec, colour = design, shape = design, linetype = design))) +
  geom_hline(yintercept = 0.192, linetype = "13", colour = B_COL, linewidth = 0.4) +
  labs(title = "B.  Coordinate-wise departure",
       x = expression(paste("visible effect size  ", tau))) +
  annotate("text", x = 7.3, y = 0.225, hjust = 1, label = "ceiling 0.19",
           size = 2.4, colour = B_COL) +
  annotate("segment", x = 3.05, xend = 3.35, y = 0.470, yend = 0.470,
           colour = A_COL, linewidth = 0.35,
           arrow = arrow(length = unit(0.04, "in"), type = "closed", ends = "first")) +
  annotate("text", x = 3.45, y = 0.470, hjust = 0, label = "still rising",
           size = 2.4, colour = A_COL) +
  theme(legend.position = "none")   # same encoding as panel A; one legend is enough

## (C) the EDGE advantage REVERSES with dimension
lg <- rbind(
  data.frame(tau = a$tau, power = a$edge, panel = a$panel, basis = "EDGE"),
  data.frame(tau = a$tau, power = a$dec,  panel = a$panel, basis = "decile"))
lg$grp <- paste(lg$panel, lg$basis)
lg$design <- ifelse(lg$panel == "A", "design A  (p = 5)", "design B  (p = 100)")
lg <- lg[order(lg$grp, lg$tau), ]
pC <- ggplot(lg, aes(tau, power, colour = basis, linetype = design,
                     shape = basis, group = grp)) +
  ek_nominal() + geom_line(linewidth = 0.7) + geom_point(size = 1.9) +
  scale_colour_manual(values = c("EDGE" = E_COL, "decile" = "grey45"), name = NULL) +
  scale_shape_manual(values = c("EDGE" = 19, "decile" = 1), name = NULL) +
  scale_linetype_manual(values = c("design A  (p = 5)" = "solid",
                                   "design B  (p = 100)" = "22"), name = NULL) +
  scale_y_continuous(limits = c(0, 0.63), breaks = seq(0, 0.6, 0.15)) +
  labs(title = "C.  The EDGE advantage reverses",
       x = expression(paste("visible effect size  ", tau)), y = NULL) +
  theme_ek(base_size = 9) +
  theme(legend.position = "none", legend.justification = c(0, 1),
        legend.key.width = unit(20, "pt"), legend.text = element_text(size = rel(0.8)),
        legend.background = element_blank(), legend.box = "vertical",
        legend.spacing.y = unit(0, "pt"),
        plot.title = element_text(hjust = 0, face = "bold", size = rel(0.95)),
        plot.subtitle = element_blank())
## (no corner annotation in panel C: at the true OUP text width it collided with the
##  four-row legend, and the subtitle already carries the message)

pl <- pA | pB | pC
## GEOMETRY: build at the OUP Series C text width EXACTLY (446.70827pt = 6.181in),
## so includegraphics[width=textwidth] applies NO rescaling and the 9pt base type
## stays 9pt on the page. Building at Springer's EK_W2 and letting LaTeX shrink it was
## rendering the annotations at about 5pt.
OUP_W <- 446.70827/72.27
ggsave("../Fig/fig3_power.pdf", pl, width = OUP_W, height = OUP_W * 0.375,
       device = cairo_pdf)
ggsave("../Fig/fig3_power.png", pl, width = OUP_W, height = OUP_W * 0.375,
       dpi = EK_DPI)

## the collapse arithmetic, printed for the caption
ia  <- approx(a$nsqrtn[a$panel == "A"], a$dec[a$panel == "A"],
              xout = a$nsqrtn[a$panel == "B"], rule = 1)$y
gap <- abs(ia - a$dec[a$panel == "B"])
cat(sprintf("panel A mean |gap| = %.4f (MC SE of a single difference ~0.023)\n",
            mean(gap, na.rm = TRUE)))
## Guards. Each asserts a claim the figure or its caption actually makes, so a future
## change to the grid cannot silently falsify one of them.
stopifnot(mean(gap, na.rm = TRUE) < 0.03)                          # panel A: the collapse
stopifnot(max(d$dec[d$alt == "quad" & d$panel == "B"]) <= 0.20)    # panel B: the ceiling
decA <- d$dec[d$alt == "cubic" & d$panel == "A"]
edgeA <- d$edge[d$alt == "cubic" & d$panel == "A"]
decB <- d$dec[d$alt == "cubic" & d$panel == "B"]
edgeB <- d$edge[d$alt == "cubic" & d$panel == "B"]
## NOT all(edgeA > decA): it is FALSE at the first point (0.049 vs 0.057).
stopifnot(all(edgeB <= decB), sum(edgeA > decA) >= 4)              # panel C: the reversal
## B6 monotonicity threshold rel = gamma/sd(eta0) >= sqrt(6)/3; sd(eta0) = 0.794 (A), 1.50 (B).
relA <- d$gamma[d$alt == "cubic" & d$panel == "A"]/0.794
relB <- d$gamma[d$alt == "cubic" & d$panel == "B"]/1.50
stopifnot(sum(relA > sqrt(6)/3) == 1, sum(relB > sqrt(6)/3) == 3)  # as stated in section 4.1
cat("wrote ../Fig/fig3_power.pdf/.png
")
