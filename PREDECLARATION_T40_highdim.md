# T40 — The corrected test across p/n, against packaged high-dimensional rivals

Pre-declaration, written 2026-10-02 before any T40 data exist. Its SHA-256 is recorded in the project tracker
when it is frozen; any later change is a dated amendment at the end of this file, never an edit above it.
Background and sources: `paper_seriesC/bimj/HIGHDIM_LITERATURE_2026-10-02.md`.

## 1. Questions

1. How do the uncorrected and corrected grouped tests behave as p/n grows through the proportional regime,
   including past the point where the maximum likelihood estimate stops existing (about p/n = 0.39 at signal 1.5)?
2. How do they compare, on the same datasets, with the packaged high-dimensional goodness-of-fit tests?
3. Do the answers depend on the covariate correlation or on a sparse signal?

Each test is judged against its OWN null. The rivals refit their own models (GRP: its lasso; PLStests: a post-lasso
refit without intercept; BAGofT: a cross-validated ridge at every split), so a rival's "size" is the rate at which
it rejects data generated from a correct logistic model, not a statement about the deployed ridge fit.

## 2. Design

- n = 400 throughout. p ∈ {20, 40, 100, 140, 180}, so p/n ∈ {0.05, 0.10, 0.25, 0.35, 0.45}.
- Covariates N_p(0, Σ) with Σ_jk = ρ^|j−k|.
- Signal strength fixed at sd(η₀) = γ = 1.5 in every cell, so only p/n, ρ and the signal pattern move.
  - **Dense**: β ∝ the design-B cycle (0.35, −0.30, 0.25, −0.20, 0.15) repeated over p, rescaled to γ = 1.5.
  - **Sparse**: β ∝ (1, 1, 1, 1, 1, 0, …, 0), rescaled to γ = 1.5 (the literature's support, at our signal strength).
- No intercept (event fraction one half, as in the paper).
- Departures, each of unit variance and added to η₀ with amplitude a:
  - **index**: a·(z² − 1)/√2 with z = η₀/γ (aligned with the index);
  - **coord**: a·(x₁² − 1)/√2 (one coordinate);
  - **inter**: a·x₁x₂/√(1 + ρ²) (an interaction).
- Panels:
  - **Main (M)**: dense, ρ = 0.7, all five p; null plus each departure at a ∈ {0.5, 1.0}.
  - **Correlation (R)**: dense, ρ ∈ {0.4, 0.8}, all five p; null plus each departure at a = 1.0.
  - **Sparse (S)**: sparse, ρ = 0.7, all five p; null plus each departure at a = 1.0.
- The penalty is chosen by 10-fold `cv.glmnet(alpha = 0)` on each dataset (one-standard-error rule), converted
  to the theory scale λ = n·λ_g, and held fixed inside the bootstrap, as a practitioner would.

## 3. Methods (all on the same dataset in each replicate)

| id | test | reference |
|---|---|---|
| HL | decile statistic on the ridge fit, G = 10 | textbook χ²_{G−2} |
| EDGE_u | EDGE statistic on the ridge fit | the MLE covariance I − UF⁻¹Uᵀ, by Monte Carlo (20,000 draws) — the strongest uncorrected form, as in Table 5 |
| SC.HL, SC.EDGE | the paper's corrected tests | prepivoted parametric bootstrap from π(β̃), N_B = 499 |
| GRP5, GRP1 | `GRPtests::GRPtest(fam = "binomial")` at its default nsplits = 5 and at 1 | package |
| PLS | `PLStests::PLStests(family = "binomial")`, defaults; primary p-value `T_cauchy` (also stored: `T_alpha`, `T_hmp`, `T_beta`) | package |
| BAG | `BAGofT::BAGofT(testGlmnet(y ~ ., alpha = 0), data)` at its published defaults (nsplits = 100, nsim = 100, ne = ⌊5√n⌋); p-value `p.value` (also `p.value3`; never `pmin`) | package |

BAG runs in **panel M only** (cost: about 51 min per test at p = 100). Every other method runs in every cell.
No rival is run at a non-default setting except GRP1, which is reported beside the default, not instead of it.

## 4. Replications and seeds

- Null cells: 1000 replications; departure cells: 500. BAG: 150 replications in every panel-M cell
  (its own column; the other methods keep their full counts on the same seeds).
- Replicate r of cell c uses `set.seed(c·10⁶ + r)` inside the worker (Mersenne-Twister); results do not depend
  on the number of workers. Each row stores sum(y) and a checksum of X[1:3, 1:3] as a data fingerprint.
- Every replicate is written to disk as it finishes; a restarted run skips completed (cell, replicate) pairs.

## 5. Outcomes and rules

- Rejection at p < 0.05. Rates with exact binomial 95% intervals; Monte Carlo SE stated beside every verdict.
- A test "holds its level" in a cell when the rate lies in [0.03, 0.08] (the paper's range). A rate outside the
  range or on its edge is re-estimated on an independent block (seeds offset by 5·10⁵) before it is discussed.
- No verdict of "no difference" between two tests is drawn without the SE of the difference; ceiling or floor
  cells (all rates > 0.95 or < 0.06 under a departure) are reported but not used to rank tests.
- Also recorded per replicate: λ̂, whether `glm` converges and its largest |coefficient| (MLE existence proxy),
  the generator's range, and the wall time of each method.

## 6. Predictions (declared now; each will be marked confirmed, refuted or inconclusive)

- P1. HL rejects correct models at a rate near 1 in every cell with p/n ≥ 0.10.
- P2. SC.HL and SC.EDGE hold their level for p/n ≤ 0.25 in panel M (where the paper already has evidence).
  For p/n = 0.35 and 0.45 no prediction: this is the open question.
- P3. Under `coord`, the power of the grouped tests falls as p/n grows (dilution); under `index` it falls less.
- P4. GRP is more powerful than the grouped tests under `coord`; the grouped tests are more powerful under `index`.
- P5. GRP's size depends on ρ and nsplits (the literature disagrees: ≤ 0.06 vs 0.31–0.42).

## 7. Budget

A timing pilot (5 replicates of the most expensive cell of each kind) runs first; if the projected wall time
on 22 workers exceeds 7 days, the BAG replication count is cut (never its splits) and the cut is recorded here
as an amendment before the main run starts.

## Amendments

**Amendment 1 — 2026-10-02, after the timing pilot and before any main-run data.** The pilot measured the
cheap methods at 6 s (p = 20) to 44 s (p = 180) per replicate, about 370 CPU-hours for the whole grid, and
BAGofT at about 51 minutes per test at p = 100 (roughly 200 times the corrected test's 13.7 s). BAGofT at the
Section 4 counts would need about 4,400 CPU-hours (8 days on 22 workers). On the author's instruction BAGofT is
therefore run **after** the main grid, on a reduced set of cells: panel M, p ∈ {20, 100, 180}, the null and the
`index` and `coord` departures at a = 1.0 — 9 cells × **100 replicates**, always at its published defaults
(nsplits = 100, nsim = 100, ne = ⌊5√n⌋); never fewer splits. Its measured cost per test is reported beside the
corrected test's as a practical limitation. Every other part of the protocol is unchanged.

**Amendment 2 — 2026-10-02, while T40 runs and before any T40 result has been read.** On the author's
instruction a maximum-likelihood baseline is added as a separate script, `R/T41_mle_baseline.R`, on exactly the
T40 datasets (it reads T40's own generator and seeds). For every replicate of every cell, the unpenalized logistic
model with intercept is fitted by `glm.fit` (maxit = 100) and two classical tests are applied to it: the textbook
decile HL statistic against χ²_{G−2} (**HL_mle**), and the EDGE statistic against the MLE covariance
I − UF⁻¹Uᵀ by Monte Carlo (**EDGE_mle**, 20,000 draws). The MLE is declared to **exist** in a replicate when
`glm.fit` converges, every fitted probability lies in (10⁻⁸, 1 − 10⁻⁸) and the largest |coefficient| is below 50;
rates are reported among the replicates where it exists, beside the fraction of replicates where it does not
(where a practitioner would have no valid MLE test at all).
Prediction **P6**: HL_mle and EDGE_mle reject correct models above 0.08 at p/n ≥ 0.25 (as in the paper's
design B, 0.178 and 0.284) and are near nominal at p/n = 0.05; the MLE fails to exist in most replicates at
p/n = 0.45.

**Amendment 3 — 2026-10-02, execution only, before any BAGofT result exists.** The timing pilot's BAGofT
workers lost their cluster socket connection after about two hours. A standalone diagnostic ran one test at
each of p = 20, 100, 180 without error (37, 73 and 111 minutes under load; peak memory about 220 MB), so the
fault is the long-lived socket, not the method. BAGofT is therefore run with **one R process per test**
(`R/T40_bag_one.R`, driven by `R/T40_bag_runner.R`, 22 at a time). A process that ends without writing its row
is retried once; a second failure is recorded as a failure and counted in the report. Nothing about the design,
the replicate counts, the seeds or BAGofT's settings changes.

**Amendment 4 — 2026-10-02, after an interim look at the complete panel-M cells, at the author's request.**
The interim look showed PLStests leading on two departures. The T40 truth has no intercept, which matches
PLStests' intercept-free refit, so the design may favour it. **T42** (`R/T42_intercept.R`, cells 201–212) repeats
panel M with intercept α = −1.530973, which makes the event fraction exactly 0.25 under the null: dense β,
ρ = 0.7, p ∈ {20, 100, 180}, the null (1000 replicates) and `index`, `coord`, `inter` at a = 1.0 (500 each);
every T40 method and the T41 MLE baseline on the same datasets; no BAGofT. Because this amendment follows an
interim look, T42 is reported as **confirmatory of a specific concern raised by the data**, not as part of the
original design. Size-adjusted power (Section 5 convention: threshold = the 5th percentile of the method's
p-values in the matched null cell) is added to the analysis for all panels.
Prediction **P7**: with an intercept in the truth, PLStests rejects correct models well above 0.08 at every p;
SC.HL and SC.EDGE stay within [0.03, 0.08] at p/n ≤ 0.25; HL stays far above the range.

**Amendment 5 — 2026-10-02, after an interim look at 18 of the 20 null cells.** The interim look put four
corrected-test size estimates above 0.08 (ρ = 0.8 at p = 140 and 180; sparse signal at p = 100). The paper's rule
is that a rate outside the range, or on its edge, is re-estimated on an independent block before it is discussed.
To avoid re-testing only the cells that looked bad (which would bias the pooled estimate downward by regression
to the mean), the second block covers **every null cell with p/n ≥ 0.25 in every panel**: cells 15, 22, 29, 44,
48, 52, 64, 68, 72, 84, 88, 92, 1000 further replicates each, replicates 1001–2000 (seeds cell·10⁶ + 1001…, never
used by the main run), all methods, written to `data/T40_recheck_pvalues.csv` (`T40_highdim_grid.R recheck`).
Reporting: the two blocks are shown separately and pooled (2000 replicates per cell); a cell is called liberal
if the pooled rate's exact 95% interval lies above 0.05 and outside-range if the pooled rate exceeds 0.08.
Rivals' sizes are pooled the same way. The run order becomes T40 main → T42 → recheck → BAGofT.

**Amendment 6 — 2026-10-03, after T42 completed, at the author's request.** The T42 null cell at p/n = 0.25
(cell 205) gave SC.HL 0.087 and SC.EDGE 0.096, outside the range. Under the same rule it gets an independent
second block: replicates 1001–2000 (`T42_intercept.R recheck`, `data/T42_recheck_pvalues.csv`), reported
separately and pooled exactly as in Amendment 5. The p/n = 0.05 cell (inside the range) and the p/n = 0.45 cell
(about 0.20, more than 15 standard errors out, so a second block cannot change its reading) are not re-run.
Run order: recheck → T42 recheck → BAGofT.

**Amendment 7 — 2026-10-05, at the author's request, after the p/n = 0.05 and 0.25 BAGofT cells were complete
and before any p/n = 0.45 BAGofT test had finished.** BAGofT took a median of 42 and 86 minutes per test at
p/n = 0.05 and 0.25 (a projected 150 minutes or more at 0.45, about 1.5 days for the remaining 300 tests). The
submitted paper reports BAGofT at p/n = 0.05 and 0.25 only, the range over which the corrected test is claimed
valid; the p/n = 0.45 BAGofT cells keep running and will be reported in full, whatever they show, at revision or
in the archive. No other cell, method or setting changes.
