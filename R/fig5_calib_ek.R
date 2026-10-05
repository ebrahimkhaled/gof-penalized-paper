# Figure 5 for the paper — GlaucomaM calibration, EK house style.
# Replaces the default base-R plot. One declarative message: the ridge fit is too
# FLAT, and that flatness is what the uncorrected test is reading as misfit.
suppressPackageStartupMessages({library(glmnet); library(TH.data); library(ggplot2)})
## Run from this directory:  cd paper_seriesC/R && Rscript <script>
## _ek_theme.R ships WITH the archive, so the figures rebuild without the author tree.
source("_ek_theme.R")

data("GlaucomaM", package = "TH.data")
X <- scale(as.matrix(GlaucomaM[, setdiff(names(GlaucomaM), "Class")]))
X <- X[, apply(X, 2, function(c) all(is.finite(c)) & sd(c) > 0), drop = FALSE]
y <- as.numeric(GlaucomaM$Class == "glaucoma"); n <- nrow(X); p <- ncol(X)
set.seed(1); cv <- cv.glmnet(X, y, family = "binomial", alpha = 0)
g  <- glmnet(X, y, family = "binomial", alpha = 0, lambda = cv$lambda.1se)
ph <- as.numeric(predict(g, newx = X, type = "response"))

G <- 10
grp  <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
idx  <- split(seq_len(n), grp)
pbar <- vapply(idx, function(I) mean(ph[I]), 0.0)
obs  <- vapply(idx, function(I) mean(y[I]),  0.0)
k    <- vapply(idx, function(I) sum(y[I]),   0.0)
m    <- vapply(idx, function(I) length(I),   0.0)
# Wilson interval
z <- 1.96; cen <- (k + z^2/2)/(m + z^2); hw <- z*sqrt(k*(m-k)/m + z^2/4)/(m + z^2)
## Direction of the shrinkage pattern, stated so the colours cannot be mislabelled:
##   LOW-risk deciles  -> penalty pulls the prediction UP toward the mean
##                     -> predicted TOO HIGH -> observed sits BELOW the diagonal
##   HIGH-risk deciles -> penalty pulls the prediction DOWN toward the mean
##                     -> predicted TOO LOW  -> observed sits ABOVE the diagonal
d <- data.frame(pbar, obs, lo = pmax(cen-hw,0), hi = pmin(cen+hw,1),
                side = ifelse(obs > pbar, "predicted too low", "predicted too high"))
stopifnot(all(d$side[d$pbar < 0.5] == "predicted too high"),
          all(d$side[d$pbar > 0.6] == "predicted too low"))   # the shrinkage signature
## Exact split, asserted so the alt text cannot drift again: SIX deciles below the
## diagonal and FOUR above, and the sign changes exactly once (no zig-zag). Decile 6 sits
## almost on the line (obs 0.526 vs pred 0.532) and falls in the [0.5,0.6] window the
## guard above deliberately leaves unchecked -- which is how the 5/4 miscount survived.
stopifnot(sum(d$obs < d$pbar) == 6, sum(d$obs > d$pbar) == 4,
          sum(diff(sign(d$obs - d$pbar)) != 0) == 1)
lp <- log(pmin(pmax(ph,1e-8),1-1e-8)/(1-pmin(pmax(ph,1e-8),1-1e-8)))
slope <- coef(glm(y ~ lp, family = binomial))[2]

ek_blue <- "#0072B2"; ek_red <- "#D55E00"; ek_grey <- "grey45"

pl <- ggplot(d, aes(pbar, obs)) +
  # above the diagonal = predicted too low (blue); below = predicted too high (orange)
  annotate("polygon", x = c(0,1,0),   y = c(0,1,1),   fill = ek_blue, alpha = 0.045) +
  annotate("polygon", x = c(0,1,1),   y = c(0,1,0),   fill = ek_red,  alpha = 0.045) +
  ek_diag() +
  geom_linerange(aes(ymin = lo, ymax = hi, colour = side),
                 linewidth = 0.55, alpha = 0.85, show.legend = FALSE) +
  geom_line(colour = "grey60", linewidth = 0.35, linetype = "22") +
  # SHAPE is redundant with COLOUR on purpose: the figure's whole message is which side
  # of the diagonal a decile falls on, and colour alone loses that in greyscale printing
  # and for colour-blind readers. Filled circle = predicted too high, open triangle = too low.
  geom_point(aes(colour = side, shape = side), size = 2.5, stroke = 0.9, fill = "white") +
  scale_colour_manual(values = c("predicted too low"  = ek_blue,
                                 "predicted too high" = ek_red), name = NULL) +
  scale_shape_manual(values = c("predicted too low"  = 24,   # filled-outline triangle
                                "predicted too high" = 21),  # filled-outline circle
                     name = NULL) +
  # annotation sits in the empty corner but is coloured to match the points it describes
  # (the two side keys moved to the caption: RSS puts keys there, not in the image)
  annotate("text", x = 0.50, y = 1.035, hjust = 0.5, size = 2.75, colour = "grey20",
           label = sprintf("apparent calibration slope %.2f: too flat", slope)) +
  coord_equal(xlim = c(0,1), ylim = c(0,1.06), expand = FALSE) +
  labs(x = "mean predicted probability within decile",
       y = "observed event rate within decile") +
  theme_ek(base_size = 9) +
  theme(legend.position = "none",
        plot.margin = margin(6, 8, 4, 4))

## GEOMETRY: the paper sets this at 0.72 x textwidth of the Biometrical Journal page
## (0.72 x 423.72342pt = 4.222in). BUILD AT THAT SIZE -- do not build large and let
## \includegraphics shrink it, or the type shrinks with the plot.
BJ_W <- 423.72342/72.27
ggsave("../Fig/fig5_calibration.pdf", pl, width = 0.72*BJ_W, height = 0.72*BJ_W*1.015, device = cairo_pdf)
ggsave("../Fig/fig5_calibration.png", pl, width = 0.72*BJ_W, height = 0.72*BJ_W*1.015, dpi = EK_DPI)
cat(sprintf("n=%d p=%d kappa=%.3f slope=%.3f lambda.1se(theory)=%.1f\n",
            n, p, p/n, slope, cv$lambda.1se*n))
cat("wrote ../Fig/fig5_calibration.pdf/.png\n")
