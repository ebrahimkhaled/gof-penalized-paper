# Replication archive
[![DOI](https://zenodo.org/badge/1331749712.svg)](https://doi.org/10.5281/zenodo.21900114)

**Shrinkage invalidates the Hosmer–Lemeshow test: goodness of fit for penalized
logistic regression, with an application to glaucoma diagnosis**

Ebrahim Khaled Ebrahim (ORCID 0009-0006-7839-8778)
Department of Applied Statistics, Faculty of Business, Alexandria University, Egypt

Submitted to *Journal of the Royal Statistical Society, Series C*.

---

## What is here

Everything needed to reproduce every number, table and figure in the paper and its
supplement. All computation is in R.

```
paper.tex, supplement.tex      the manuscript sources
oup-authoring-template.cls     the OUP class (build: pdflatex twice, no bibtex)
Fig/                           figures as they appear in the paper
R/                             all scripts + _ek_theme.R (self-contained)
*.R                            the simulation and application scripts
*_pvalues*.csv                 per-replication p-values — power study + rivals only
                               (T25_power_pvalues.csv, T24_rivals_pvalues.csv)
*_results.rds                  saved objects for every run
```

## Data

The `GlaucomaM` data are **not redistributed here**. They are publicly available in the R
package `TH.data` on CRAN and are used without modification:

```r
install.packages("TH.data"); data("GlaucomaM", package = "TH.data")
```

The two secondary datasets are `GlaucomaMVF` (package `ipred`) and `Sonar`
(package `mlbench`), also from CRAN.

## Reproducing the paper

Current numbering: **Tables** 1 oracle grouping · 2 corrected size · 3 power · 4 power cost ·
5 glaucoma grid. **Figures** 1 penalty path · 2 cost and gain · 3 power · 4 glaucoma calibration.
Figure files and their scripts carry the same numbers.

| Paper item | Script | Cost |
|---|---|---|
| Proposition 1 numerical check | `corrected_test.R` | seconds |
| Covariance ordering Omega_K >= Omega_MLE | `check_omega_sign.R` | seconds |
| Figure 1, the penalty path | `lambda_sweep.R`, `hd_sweep.R` -> `R/fig1_penaltypath_ek.R` | minutes |
| Table 1, oracle vs fitted grouping | `oracle_knownnull.R` | minutes |
| Table 2, corrected size | `T23_rerun.R` | 11 min, 22 workers |
| Table 2, the lambda=50 row (pooled) | `verify_E13.R` | 15 min, 22 workers |
| Tables 3-4 and Figure 2, power and its cost | `T25_power_highB.R` -> `R/fig2_costgain_ek.R` | ~37 min, 22 workers |
| Figure 3, the visible effect size | `T22_tau.R`, `plot_tau.R` -> `R/fig3_power_ek.R` | minutes |
| Table 5 and Figure 4, glaucoma | `glaucoma_deep.R` -> `R/fig4_calib_ek.R` | minutes |
| Uncorrected variants on GlaucomaM | `naive_glaucoma_both.R` | seconds |
| Second dataset (GlaucomaMVF) | `extract_T26b.R` | minutes |
| The one-fit screen (S2.5) | `shrinkage_screen.R` | minutes |
| Conditioning and the generator | `glaucoma_conditioning.R`, `check_generator.R` | seconds |
| Table S3, residual geometry | `omega_spectrum.R` | 26 s, 22 workers |
| S4.3, overshoot intervention | `overshoot_test.R` | 4 min, 22 workers |
| Table S1, the failed repairs | `round2.R`, `prepivot_cost.R` | minutes |

All four figure scripts run from `R/` and source `_ek_theme.R`, which ships with this archive;
each ends in `stopifnot()` guards asserting the claims its figure makes.

Parallel scripts use PSOCK clusters sized `detectCores() - 2` and set their stream with
`parallel::clusterSetRNGStream`; the stream seed is written at the top of each script and
listed in Table S3 of the supplement. **Results do not depend on the number of workers**,
because each replication is seeded individually by `set.seed(seed)` inside the worker.

## Two things a replicator should know

**1. The penalty scale.** Two conventions are in play. The theory scale is
`lambda` in `-loglik + 0.5*lambda*||beta||^2`; `glmnet` minimises
`-loglik/n + 0.5*lambda_g*||beta||^2`. So `lambda = n * lambda_g`. Every table states
which scale it uses. The identity was verified numerically to 1.6e-9.

**2. Seed blocks are not interchangeable.** `T23_rerun.R` dispatches
`seeds <- seq_len(B) + 1000L`, i.e. seeds 1001–2000, *not* 1–1000. `verify_E13.R` uses
1–1500 and 100001–101500, which is why those blocks are independent of it and can be
pooled. If you re-run a cell on seeds 1–1000 expecting to match Table 2, you will not.

## Per-replication p-values

Two experiments archive their individual p-values rather than only the rejection rate:
the power study (`T25_power_pvalues.csv`) and the head-to-head comparison
(`T24_rivals_pvalues.csv`). For those, any cell can be re-derived and any Monte Carlo
standard error recomputed independently, and any threshold re-applied. The remaining
experiments archive rejection indicators and their fitted objects (`*_results.rds`),
which reproduce each reported rate but do not permit re-thresholding.

For example, to recover the size-adjusted power of the EDGE basis in design B:

```r
d   <- read.csv("T25_power_pvalues.csv")
b   <- subset(d, panel == "B")
thr <- quantile(b$p_pre_edge[b$gamma == 0], 0.05, type = 1)
tapply(b$p_pre_edge <= thr, b$gamma, mean)
```

## Software

The EDGE basis is implemented as `edge.gof()` in the R package `ebrahim.gof` (2.4.0) on
CRAN. **The shrinkage correction itself is not yet in a released version of that package** —
the reference implementation is in this archive, self-contained so that reproduction does
not depend on any package version. The EDGE statistic computed here was verified identical
to `ebrahim.gof::edge.gof` to a difference of 0.00e+00 on three datasets (`verify_edge.R`).

Session: R 4.4.x on Windows 11, packages `glmnet`, `TH.data`, `ggplot2`, `patchwork`,
`parallel`.

## Licence

Code: MIT. Manuscript text and figures: CC BY 4.0.
