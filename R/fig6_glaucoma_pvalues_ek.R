# Figure 6 for the paper -- THE GLAUCOMA MODEL, UNCORRECTED AGAINST CORRECTED. EK house style.
#
# One message: at every penalty and grouping the uncorrected test sits at or near its resolution floor,
# while the corrected tests scatter around the 5% line. The exact values are SI Table "tab:glaucoma".
# Data (archived): glaucoma_deep_results.rds (corrected, N_B = 499) and naive_glaucoma_both.csv
# (uncorrected, referred to the maximum likelihood covariance by 2e5 Monte Carlo draws; 0 = below 1e-5).
## Run from this directory:  Rscript fig6_glaucoma_pvalues_ek.R
suppressPackageStartupMessages({library(ggplot2)})
source("_ek_theme.R")
D <- c("../data", "../../gof-penalized-paper/data"); D <- D[dir.exists(D)][1]

corr <- readRDS(file.path(D, "glaucoma_deep_results.rds"))
naive <- read.csv(file.path(D, "naive_glaucoma_both.csv"))
key <- c("lambda.min" = "min", "0.5*1se" = "half", "half.1se" = "half", "lambda.1se" = "1se",
         "2*1se" = "twice", "twice.1se" = "twice")
corr$lam <- key[corr$lambda]; naive$lam <- key[naive$lambda]
d <- rbind(data.frame(lam = naive$lam, G = naive$G, test = "Uncorrected", p = naive$Om_p),
           data.frame(lam = corr$lam, G = corr$G, test = "SC.HL", p = corr$corr_dec),
           data.frame(lam = corr$lam, G = corr$G, test = "SC.EDGE", p = corr$corr_edge))
FLOOR <- 1e-5
d$y <- pmax(d$p, FLOOR)                          # below the Monte Carlo resolution: drawn on the floor
d$lam <- factor(d$lam, levels = c("min", "half", "1se", "twice"),
                labels = c("lambda[min]", "0.5*lambda['1se']", "lambda['1se']", "2*lambda['1se']"))
d$test <- factor(d$test, levels = c("Uncorrected", "SC.HL", "SC.EDGE"),
                 labels = c("A. Uncorrected", "B. SC.HL", "C. SC.EDGE"))

pl <- ggplot(d, aes(lam, y, shape = factor(G), colour = test)) +
  geom_hline(yintercept = 0.05, linetype = "13", colour = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = FLOOR, linetype = "solid", colour = "grey80", linewidth = 0.3) +
  geom_point(size = 2.1, stroke = 0.7, position = position_dodge(width = 0.5)) +
  scale_y_log10(limits = c(FLOOR * 0.8, 1), breaks = c(1e-5, 1e-4, 1e-3, 0.01, 0.05, 0.2, 1),
                labels = c("<1e-5", "1e-4", "0.001", "0.01", "0.05", "0.2", "1")) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_shape_manual(values = c(`5` = 1, `10` = 19, `20` = 2), name = "groups") +
  scale_colour_manual(values = c("A. Uncorrected" = "#56B4E9", "B. SC.HL" = "grey35", "C. SC.EDGE" = "#009E73"),
                      guide = "none") +
  facet_wrap(~ test, nrow = 1) +
  labs(x = "ridge penalty", y = expression(italic(p)*"-value (log scale)")) +
  theme_ek(base_size = 9) +
  theme(legend.position = "bottom", strip.text = element_text(hjust = 0), panel.grid.minor = element_blank())

W <- 6.181
for (dir in c("../Fig", "../bimj/Fig")) if (dir.exists(dir))
  ggsave(file.path(dir, "fig6_glaucoma_pvalues.pdf"), pl, width = W, height = 2.6, device = cairo_pdf)
ggsave("../Fig/fig6_glaucoma_pvalues.png", pl, width = W, height = 2.6, dpi = 220)

## the claims the figure carries in the text (Section 6.2)
u <- d$p[d$test == "A. Uncorrected"]; s <- d$p[d$test != "A. Uncorrected"]
stopifnot(all(u < 0.006), sum(u < FLOOR) >= 7,                     # uncorrected rejects everywhere
          abs(min(s) - 0.014) < 1e-9, abs(max(s) - 0.678) < 1e-9,  # corrected range across the grid
          abs(d$p[d$test == "B. SC.HL" & d$lam == "lambda['1se']" & d$G == 10] - 0.026) < 1e-9)
cat("figure 6 claims verified\n")
