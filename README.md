# Replication archive

[![DOI](https://zenodo.org/badge/1331749712.svg)](https://doi.org/10.5281/zenodo.21900114)

**The Hosmer–Lemeshow test after shrinkage: a valid goodness-of-fit test for penalized
logistic regression, applied to glaucoma diagnosis**

Ebrahim Khaled Ebrahim (ORCID 0009-0006-7839-8778) and Ahmed El-Kotory
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
R/                             the two figure scripts + _ek_theme.R (self-contained)
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

Tables are named by content as well as number, because the numbering shifts if a table is
added or moved. Current numbering is Table 1 uncorrected size, 2 oracle grouping,
3 corrected size, 4 power, 5 power cost, 6 glaucoma grid.

| Paper item | Script | Cost |
|---|---|---|
| Proposition 1 numerical check | `corrected_test.R` | seconds |
| Covariance ordering Ω_K ⪰ Ω_MLE | `check_omega_sign.R` | seconds |
| Table 1, uncorrected size | `lambda_sweep.R`, `hd_sweep.R` | minutes |
| Table 2, oracle vs fitted grouping | `oracle_knownnull.R` | minutes |
| Table 3, corrected size | `T23_rerun.R` | 11 min, 22 workers |
| Table 3, the λ=50 row (pooled) | `verify_E13.R` | 15 min, 22 workers |
| Tables 4 and 5, power and its cost | `T25_power_highB.R` | ~37 min, 22 workers |
| Figure 1, data then rendering | `T22_tau.R` → `R/fig1_power_ek.R` | minutes |
| Table 6 and Figure 2, glaucoma | `glaucoma_deep.R`, `R/fig2_calib_ek.R` | minutes |
| Uncorrected variants on GlaucomaM | `naive_glaucoma_both.R` | seconds |
| Table S2, residual geometry | `omega_spectrum.R` | 26 s, 22 workers |
| §S4.3, overshoot intervention | `overshoot_test.R` | 4 min, 22 workers |
| Table S1, the failed repairs | `round2.R`, `cand_C.R` | minutes |

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

The corrected test is implemented as `gof.pen()` in the R package `ebrahim.gof` on CRAN;
it accepts a fitted `glmnet` object directly. The scripts here use a self-contained
implementation so that the archive does not depend on a package version; it was verified
identical to `ebrahim.gof::edge.gof` to a difference of 0.00e+00 on three datasets
(`verify_edge.R`).

Session: R 4.4.x on Windows 11, packages `glmnet`, `TH.data`, `ggplot2`, `patchwork`,
`parallel`.

## Citing this archive

Concept DOI (always resolves to the latest version): **10.5281/zenodo.21900114**

## Licence

Code: MIT. Manuscript text and figures: CC BY 4.0.
