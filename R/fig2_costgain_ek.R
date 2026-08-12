# Figure 2 for the Series C paper — THE COST AND THE GAIN. EK house style.
#
# The paper's most unusual virtue is that it reports a power cost it could have hidden.
# As a table that virtue is four numbers; as a shape it is two arrows of visibly
# different length pointing in opposite directions.
#
# Data: manifest 6.3 (T25_power_highB.R, B = 2000 per cell, MC SE <= 0.0112). GREEN.
## Run from this directory:  cd paper_seriesC/R && Rscript fig2_costgain_ek.R
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
source("_ek_theme.R")

A <- data.frame(
  alt   = rep(c("quadratic", "interaction"), each = 4),
  gamma = rep(c(0.25, 0.50, 0.75, 1.00), 2),
  dec   = c(0.0760, 0.1475, 0.2015, 0.2700,  0.0595, 0.1055, 0.1940, 0.2835),
  edge  = c(0.1150, 0.2615, 0.3755, 0.4150,  0.0750, 0.1885, 0.3555, 0.4615),
  mdec  = c(0.0890, 0.1715, 0.2350, 0.2960,  0.0630, 0.1290, 0.2275, 0.3265),
  medge = c(0.1205, 0.2675, 0.3835, 0.4080,  0.0870, 0.2065, 0.3780, 0.4885))
A$lab <- sprintf("%s  γ = %.2f", A$alt, A$gamma)
A$lab <- factor(A$lab, levels = rev(A$lab))

GREY <- "grey30"; GRN <- "#009E73"; ORG <- "#D55E00"

pA <- ggplot(A, aes(y = lab)) +
  ## COST and GAIN are drawn on separate half-rows, because the cost is so much the
  ## shorter of the two that overlaying them hides it -- which is the opposite of the
  ## point. Identical linewidth and arrow head, so the two lengths still compare honestly.
  geom_segment(aes(x = mdec, xend = dec, y = as.numeric(lab) + 0.17,
                   yend = as.numeric(lab) + 0.17), colour = ORG, linewidth = 0.6,
               arrow = arrow(length = unit(0.045, "in"), type = "closed")) +
  geom_segment(aes(x = dec, xend = edge, y = as.numeric(lab) - 0.17,
                   yend = as.numeric(lab) - 0.17), colour = GRN, linewidth = 0.6,
               arrow = arrow(length = unit(0.045, "in"), type = "closed")) +
  geom_point(aes(x = mdec, y = as.numeric(lab) + 0.17), shape = 1, size = 1.9, colour = "grey55") +
  geom_point(aes(x = dec,  y = as.numeric(lab) + 0.17), shape = 19, size = 1.9, colour = GREY) +
  geom_point(aes(x = dec,  y = as.numeric(lab) - 0.17), shape = 19, size = 1.9, colour = GREY) +
  geom_point(aes(x = edge, y = as.numeric(lab) - 0.17), shape = 19, size = 1.9, colour = GRN) +
  ## positions are numeric so the two arrows can sit on separate half-rows, so the y
  ## labels must be restored by hand -- otherwise the axis prints 2, 4, 6, 8.
  scale_y_continuous(breaks = seq_along(levels(A$lab)), labels = levels(A$lab),
                     expand = expansion(add = 0.55)) +
  scale_x_continuous(limits = c(0, 0.52), breaks = seq(0, 0.5, 0.1)) +
  labs(title = "A.  Design A: what the correction costs, and what the basis buys back",
       x = "power (rejection rate)", y = NULL) +
  theme_ek(base_size = 9) +
  theme(plot.title.position = "panel",
        plot.title = element_text(hjust = 0, face = "bold", size = rel(0.95)),
        panel.grid.major.y = element_blank())

B <- data.frame(gamma = c(0, 0.5, 1.0, 1.5),
                dec   = c(0.0640, 0.0595, 0.0635, 0.0740),
                edge  = c(0.0495, 0.0485, 0.0600, 0.0630))
se <- function(p) sqrt(p*(1-p)/2000)
nullband <- c(B$dec[1] - 2*se(B$dec[1]), B$dec[1] + 2*se(B$dec[1]))

Bl <- rbind(data.frame(gamma = B$gamma, r = B$dec,  basis = "decile"),
            data.frame(gamma = B$gamma, r = B$edge, basis = "EDGE"))
pB <- ggplot(Bl, aes(gamma, r, colour = basis, shape = basis, linetype = basis)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = nullband[1], ymax = nullband[2],
           fill = "grey60", alpha = 0.16) +
  ek_nominal() +
  geom_errorbar(aes(ymin = pmax(r - 2*se(r), 0), ymax = r + 2*se(r)),
                width = 0.06, linewidth = 0.3, show.legend = FALSE) +
  geom_line(linewidth = 0.7) + geom_point(size = 2.0) +
  scale_colour_manual(values = c(decile = GREY, EDGE = GRN), name = NULL) +
  scale_shape_manual(values = c(decile = 19, EDGE = 17), name = NULL) +
  scale_linetype_manual(values = c(decile = "solid", EDGE = "22"), name = NULL) +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05)) +
  labs(title = "B.  Design B: no detectable power over this range",
       x = expression(gamma), y = "power (rejection rate)") +
  theme_ek(base_size = 9) +
  theme(legend.position = "none", plot.title.position = "panel",
        plot.title = element_text(hjust = 0, face = "bold", size = rel(0.95)))

OUP_W <- 446.70827/72.27
pl <- pA / pB + plot_layout(heights = c(1.35, 1))
ggsave("../Fig/fig2_costgain.pdf", pl, width = OUP_W, height = OUP_W*0.78, device = cairo_pdf)
ggsave("../Fig/fig2_costgain.png", pl, width = OUP_W, height = OUP_W*0.78, dpi = EK_DPI)

## ---- guards -----------------------------------------------------------------------
## The cost arrow must point BACKWARDS at every alternative and the gain arrow FORWARDS
## at all but the first (E22: EDGE does NOT beat decile everywhere).
stopifnot(all(A$dec < A$mdec))                       # the correction always costs
stopifnot(sum(A$edge > A$dec) == nrow(A))            # on design A the basis always gains
stopifnot(mean(A$edge - A$dec) > mean(A$mdec - A$dec))  # the gain exceeds the cost
stopifnot(max(B$dec) - min(B$dec) < 0.02)            # design B really is flat
cat(sprintf("mean cost %.4f ; mean gain %.4f ; gain/cost = %.1f x\n",
            mean(A$mdec - A$dec), mean(A$edge - A$dec),
            mean(A$edge - A$dec)/mean(A$mdec - A$dec)))
cat("wrote ../Fig/fig2_costgain.pdf/.png\n")
