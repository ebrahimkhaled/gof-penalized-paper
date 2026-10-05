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

The manuscript cites release `v1.0.7`; its version DOI is the one given in the manuscript's Data
Availability Statement and on the Zenodo page of this release. Earlier releases do not contain every script listed
below; do not reproduce from them.


## What is here

Everything needed to reproduce every number, table and figure in the paper and its
Supporting Information. All simulation and estimation is in R; the tables of Section S8 and the
data behind Figure 4 are assembled from the archived p-values by two short Python scripts
(`T40_validity_map_data.py`, `T40_si_tables.py`; pandas, NumPy, SciPy).

```
R/                 every script, including the five figure scripts and _ek_theme.R
data/              saved results (*_results.rds), per-replication p-values
                   (T25_power_pvalues.csv, T24_rivals_pvalues.csv) and summary tables
Fig/               figures as they appear in the paper
sessionInfo.txt    R, package and Python versions used
dev/               helpers used only on the authors' machine (job chaining, pausing, writing the
                   SI section); they contain local paths and are NOT needed to reproduce anything
```

Most scripts write their results to the working directory; the archived copies of those results
are in `data/`. The grid scripts (`T40_highdim_grid.R`, `T41_mle_baseline.R`, `T42_intercept.R`)
write straight to `data/` and skip every (cell, replicate) already there, so run from this archive
they report "0 jobs to run". To spot-check a setting, send a run to a scratch file instead:

```
cd R
T40_OUT=check.csv T40_ONLY=15 T40_REPS=50 T40_NW=4 Rscript T40_highdim_grid.R recheck
```

and compare its `fp_y`/`fp_x` fingerprints and p-values with the same rows of
`data/T40_recheck_pvalues.csv` (cell 15, replicates 1001-1050). `T40_OUT` works the same way in
`main` mode; `smoke` (the default mode) runs one replicate of every kind of cell.

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
grouping · 2 corrected size · 3 power · 4 what each test checks and misses · 5 glaucoma grid.
**Figures** 1 penalty path · 2 cost and gain · 3 power · 4 validity map · 5 glaucoma calibration.
Figure files and their scripts carry the printed figure numbers. The Supporting Information numbers its sections and tables
S1, S2, …; conditions and proofs are its Section S1.

| Paper item | Script | Cost |
|---|---|---|
| **The test itself, as one callable function** | `shrink_gof.R` -> `shrink.gof()` | seconds |
| Proposition 2.1 numerical check | `corrected_test.R` | seconds |
| Covariance ordering Omega_K >= Omega_MLE | `check_omega_sign.R` | seconds |
| Figure 1, the penalty path | `lambda_sweep.R`, `hd_sweep.R`, `T38_hd_sweep_0124.R` -> `R/fig1_penaltypath_ek.R` | minutes |
| Table 1, oracle vs fitted grouping | `oracle_knownnull.R` | minutes |
| Table 2, corrected size | `T23_rerun.R` | 11 min, 22 workers |
| Table 2, the lambda=50 row, pooled (Table S6) | `verify_E13.R` | 15 min, 22 workers |
| Table 3 and Figure 2, power and its cost (Section 5.2) | `T25_power_highB.R` -> `R/fig2_costgain_ek.R` | ~37 min, 22 workers |
| Figure 3, the visible effect size | `T22_tau.R`, `plot_tau.R` -> `R/fig3_power_ek.R` | minutes |
| Section 4.4 and Section S6.2, residual-prediction tests in design B | `T24_rivals.R` | not recorded |
| **Section 5.4, Figure 4 and Table 4: the grid over p/n against the packaged tests** (protocol `PREDECLARATION_T40_highdim.md`, Amendments 1–7, hashes in `.sha256`) | `T40_highdim_grid.R main`, `T41_mle_baseline.R`, `T42_intercept.R main`, second blocks `T40_highdim_grid.R recheck` and `T42_intercept.R recheck`, BAGofT `T40_bag_runner.R` (one process per test, `T40_bag_one.R`) -> `T40_summarise.R`, `T40_size_adjusted.R`, `T40_validity_map_data.py` -> `R/fig4_validity_map_ek.R` | about 27 h + 7 h + 3.3 h + 0.5 h, 22 workers; BAGofT 42–86 min per test |
| Section S8 and Tables S12–S14, the full grid | `T40_si_tables.py` (writes the tables from the data) | seconds |
| Section S8, the generator diagnostics | `T43_generator_inflation_diag.R`, `T44_signal_check.R`, `T45_oracle_generator.R` (`T44_rescaled_generator_pilot.R` is the failed first attempt, kept) | 1 h on 2 workers |
| Table 5 and Figure 5, glaucoma | `glaucoma_deep.R` -> `R/fig5_calib_ek.R` | minutes |
| Uncorrected variants on GlaucomaM | `naive_glaucoma_both.R` | seconds |
| Second dataset, GlaucomaMVF (Table S11), without `lora` | `T39_mvf_nolora.R` | minutes |
| The same with `lora` as a predictor (reported as a check) | `T26_app.R` -> `extract_T26b.R` | minutes |
| The one-fit screen (Section S4) | `shrinkage_screen.R` | minutes |
| Conditioning and the generator (Section S7) | `glaucoma_conditioning.R`, `check_generator.R` | seconds |
| Table S1, the failed repairs | `round2.R`, `cand_C.R` | minutes |
| Table S3, residual geometry | `omega_spectrum.R` | 26 s, 22 workers |
| Section S5.3, overshoot intervention | `overshoot_test.R` | 4 min, 22 workers |
| Section S5.5, cost and the price of prepivoting | `prepivot_cost.R` | minutes |
| EDGE statistic against `ebrahim.gof::edge.gof` | `verify_edge.R` | seconds |
| Section 5.1 and Table S7, the bootstrap size | `T30_keycells_NB499.R` | 218 min, 22 workers |
| Section 5.3 and Table S8, a cross-validated penalty | `T31_cv_lambda_size.R`, `T32_glaucoma_design_size.R`, `T37_glaucoma_T1_replication.R` -> `summarise_cv_lambda.R` | 106 + 32 + 17 min, 22 workers; summary seconds |
| Section 3.1 and Table S9, a smooth non-ridge penalty | `T34_logcosh_size.R` -> `T34a_logcosh_check.R` | 131 min, 22 workers |
| Section 3.1 and Table S10, sensitivity to the inverse of F | `T33_finv_stability.R` | 9 min |
| Section S7, spread over bootstrap streams | `T35_seed_spread.R`, `T36_draw_construction.R` | 5 min |

The grid's per-replicate p-values are `data/T40_main_pvalues.csv`, `data/T40_recheck_pvalues.csv`,
`data/T41_mle_pvalues.csv`, `data/T42_intercept_pvalues.csv`, `data/T42_recheck_pvalues.csv` and
`data/T40_bagoft_pvalues.csv`; every row carries a fingerprint of its dataset (`fp_y`, `fp_x`), and
replicate r of cell c is generated after `set.seed(c * 1e6 + r)`. The BAGofT p/n = 0.45 cells
(29, 31, 33) were still running when this release was made (Amendment 7); the rows present for them
are partial and are not used in the manuscript.

Development checks kept for completeness, whose results the current manuscript does not report
or reports only through a later run: `confirm_P1.R` (first size confirmation of the prepivoted
procedure, superseded by Table 2), `oracle_index.R` (first oracle-index experiment, superseded by
`oracle_knownnull.R`), `T26_app_a.R` (re-measurement of development constants), `power_study.R`,
`edge_basis_study.R`.

All five figure scripts run from `R/` and source `_ek_theme.R`, which ships with this archive;
each ends in `stopifnot()` guards asserting the claims its figure makes.

Parallel scripts use PSOCK clusters sized `detectCores() - 2` and set their stream with
`parallel::clusterSetRNGStream`; the stream seed is written at the top of each script and
listed in Table S15 of the Supporting Information. **Results do not depend on the number of
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

The grid of Section 5.4 (`T40_main`, `T40_recheck`, `T41_mle`, `T42_intercept`, `T42_recheck`,
`T40_bagoft`, each `_pvalues.csv`, every row with a data fingerprint) and ten earlier experiments
archive their individual p-values rather than only the rejection rate: the
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

The self-contained copy `R/shrink_gof.R` also returns the condition number of F and the generator
diagnostics and accepts an eigenvalue floor (`finv.floor`); the CRAN function does not yet.

Session: see `sessionInfo.txt`. `GRPtests`, used by `T24_rivals.R` and the grid, was archived on
CRAN on 8 May 2022 and must be installed from the archived tarball; `PLStests` 0.1.1 was archived
on 2 September 2025 and is installed the same way (mirror github.com/cran/PLStests).

## Licence

Code: MIT. Figures: CC BY 4.0.
