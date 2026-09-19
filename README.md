# Replication archive
[![DOI](https://zenodo.org/badge/1331749712.svg)](https://doi.org/10.5281/zenodo.21900114)

**Shrinkage invalidates the Hosmer–Lemeshow test: goodness of fit for ridge
logistic regression, with an application to glaucoma diagnosis**

Ebrahim Khaled Ebrahim (ORCID 0009-0006-7839-8778)
Department of Applied Statistics, Faculty of Business, Alexandria University, Egypt

Mohammad Arashi (ORCID 0000-0002-5881-9241)
Department of Statistics, Faculty of Mathematical Sciences, Ferdowsi University of Mashhad, Iran

An earlier version of the paper is available as a preprint: https://arxiv.org/abs/2609.06413.

---


## Which DOI to use

| Purpose | DOI |
|---|---|
| Cite the archive generally (always resolves to the newest release) | `10.5281/zenodo.21900114` — the concept DOI, and what the badge above points at |
| Cite the exact code behind a manuscript | the **version** DOI of the release that manuscript names |

The manuscript cites release `v1.0.4`. Earlier releases do not contain every script listed
below; do not reproduce from them.


## What is here

Everything needed to reproduce every number, table and figure in the paper and its
Supporting Information. All computation is in R.

```
R/                 every script, including the four figure scripts and _ek_theme.R
data/              saved results (*_results.rds), per-replication p-values
                   (T25_power_pvalues.csv, T24_rivals_pvalues.csv) and summary tables
Fig/               figures as they appear in the paper
sessionInfo.txt    R and package versions used
```

Scripts write their results to the working directory; the archived copies of those results
are in `data/`.

## Data

The `GlaucomaM` data are **not redistributed here**. They are publicly available in the R
package `TH.data` on CRAN and are used without modification:

```r
install.packages("TH.data"); data("GlaucomaM", package = "TH.data")
```

The secondary datasets are `GlaucomaMVF` (package `ipred`) and `Sonar`
(package `mlbench`), also from CRAN.

## Reproducing the paper

Numbering in the current manuscript: **Propositions** 2.1, 3.1 and 4.1. **Tables** 1 oracle
grouping · 2 corrected size · 3 power · 4 power cost · 5 glaucoma grid. **Figures** 1 penalty
path · 2 cost and gain · 3 power · 4 glaucoma calibration. Figure files and their scripts carry
the same numbers. The Supporting Information numbers its sections and tables S1, S2, ….

| Paper item | Script | Cost |
|---|---|---|
| **The test itself, as one callable function** | `shrink_gof.R` -> `shrink.gof()` | seconds |
| Proposition 2.1 numerical check | `corrected_test.R` | seconds |
| Covariance ordering Omega_K >= Omega_MLE | `check_omega_sign.R` | seconds |
| Figure 1, the penalty path | `lambda_sweep.R`, `hd_sweep.R`, `T38_hd_sweep_0124.R` -> `R/fig1_penaltypath_ek.R` | minutes |
| Table 1, oracle vs fitted grouping | `oracle_knownnull.R` | minutes |
| Table 2, corrected size | `T23_rerun.R` | 11 min, 22 workers |
| Table 2, the lambda=50 row, pooled (Table S6) | `verify_E13.R` | 15 min, 22 workers |
| Tables 3-4 and Figure 2, power and its cost | `T25_power_highB.R` -> `R/fig2_costgain_ek.R` | ~37 min, 22 workers |
| Figure 3, the visible effect size | `T22_tau.R`, `plot_tau.R` -> `R/fig3_power_ek.R` | minutes |
| Section 4.4 and Section S5.2, residual-prediction tests | `T24_rivals.R` | not recorded |
| Table 5 and Figure 4, glaucoma | `glaucoma_deep.R` -> `R/fig4_calib_ek.R` | minutes |
| Uncorrected variants on GlaucomaM | `naive_glaucoma_both.R` | seconds |
| Second dataset, GlaucomaMVF (Table S11) | `T26_app.R` -> `extract_T26b.R` | minutes |
| The one-fit screen (Section S3) | `shrinkage_screen.R` | minutes |
| Conditioning and the generator (Section S6) | `glaucoma_conditioning.R`, `check_generator.R` | seconds |
| Table S1, the failed repairs | `round2.R`, `cand_C.R` | minutes |
| Table S3, residual geometry | `omega_spectrum.R` | 26 s, 22 workers |
| Section S4.3, overshoot intervention | `overshoot_test.R` | 4 min, 22 workers |
| Section S4.5, cost and the price of prepivoting | `prepivot_cost.R` | minutes |
| EDGE statistic against `ebrahim.gof::edge.gof` | `verify_edge.R` | seconds |
| Section 5.1 and Table S7, the bootstrap size | `T30_keycells_NB499.R` | 218 min, 22 workers |
| Section 5.3 and Table S8, a cross-validated penalty | `T31_cv_lambda_size.R`, `T32_glaucoma_design_size.R`, `T37_glaucoma_T1_replication.R` -> `summarise_cv_lambda.R` | 106 + 32 + 17 min, 22 workers; summary seconds |
| Section 3.1 and Table S9, a smooth non-ridge penalty | `T34_logcosh_size.R` -> `T34a_logcosh_check.R` | 131 min, 22 workers |
| Section 3.1 and Table S10, sensitivity to the inverse of F | `T33_finv_stability.R` | 9 min |
| Section S6, spread over bootstrap streams | `T35_seed_spread.R`, `T36_draw_construction.R` | 5 min |

Development checks kept for completeness, whose results the current manuscript does not report
or reports only through a later run: `confirm_P1.R` (first size confirmation of the prepivoted
procedure, superseded by Table 2), `oracle_index.R` (first oracle-index experiment, superseded by
`oracle_knownnull.R`), `T26_app_a.R` (re-measurement of development constants), `power_study.R`,
`edge_basis_study.R`.

All four figure scripts run from `R/` and source `_ek_theme.R`, which ships with this archive;
each ends in `stopifnot()` guards asserting the claims its figure makes.

Parallel scripts use PSOCK clusters sized `detectCores() - 2` and set their stream with
`parallel::clusterSetRNGStream`; the stream seed is written at the top of each script and
listed in Table S12 of the Supporting Information. **Results do not depend on the number of
workers**, because each replication is seeded individually by `set.seed(seed)` inside the
worker.

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

Ten experiments archive their individual p-values rather than only the rejection rate: the
power study (`data/T25_power_pvalues.csv`), the head-to-head comparison
(`data/T24_rivals_pvalues.csv`), and the later studies — `T30_keycells`
(both bootstrap sizes, one row per replicate), `T31_cv_lambda`, `T32_glaucoma_design`,
`T33_finv_stability`, `T34_logcosh`, `T35_seed_spread`, `T36_draw_construction` and
`T37_glaucoma_T1`. For those, any cell can be re-derived, any Monte Carlo standard error
recomputed independently, and any threshold re-applied. The remaining experiments archive
rejection indicators and their fitted objects (`*_results.rds`), which reproduce each
reported rate but do not permit re-thresholding.

For example, to recover the size-adjusted power of the EDGE basis in design B:

```r
d   <- read.csv("data/T25_power_pvalues.csv")
b   <- subset(d, panel == "B")
thr <- quantile(b$p_pre_edge[b$gamma == 0], 0.05, type = 1)
tapply(b$p_pre_edge <= thr, b$gamma, mean)
```

## Software

The corrected tests are on CRAN as `shrink.gof()` in the R package `ebrahim.gof`, from
version 2.6.0 onwards. This archive also keeps `R/shrink_gof.R`, a self-contained reference
implementation (base R and `stats` only), so that reproduction does not depend on any
package version.

```r
source("R/shrink_gof.R")
shrink.gof(X, y, lambda = 295.3, G = 10, basis = c("edge", "decile"), B = 499)
#   SC.HL    the correction on the decile (Hosmer-Lemeshow) grouping
#   SC.EDGE  the correction on the EDGE basis
```

`lambda` is on the theory scale, `lambda = n * lambda_glmnet`. Checked against Table 5 of the
paper on `GlaucomaM` at `lambda.1se` with `G = 10`: `SC.HL` returns 0.026, matching the
published value exactly, and `SC.EDGE` returns 0.038 against a published 0.034 — a difference
well inside bootstrap error at `B = 499`. The uncorrected p-value is 0.00000, as published.

The EDGE basis is implemented as `edge.gof()` in `ebrahim.gof`. The EDGE statistic computed
here was verified identical to `ebrahim.gof::edge.gof` (version 2.4.0) to a difference of
0.00e+00 on three datasets (`verify_edge.R`).

Session: see `sessionInfo.txt`. `GRPtests`, used only by `T24_rivals.R`, was archived on CRAN
on 8 May 2022 and must be installed from the archived tarball.

## Licence

Code: MIT. Figures: CC BY 4.0.
