"""Tidy size table behind the validity-map figure (T40 + T41 + T42 + the Amendment 5/6 second blocks).

One row per design x p/n x method: rejection rate of correct models at p < 0.05, its replicate count and
exact 95% interval. Where a second block exists (null cells with p/n >= 0.25) the two blocks are pooled, as
Amendments 5 and 6 pre-declare. The classical HL on the maximum likelihood fit is counted only over replicates
where that estimate exists, and the fraction where it exists is kept beside it.
"""
import os
import numpy as np, pandas as pd
from scipy.stats import beta

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
main = pd.concat([pd.read_csv(os.path.join(D, "T40_main_pvalues.csv")),
                  pd.read_csv(os.path.join(D, "T42_intercept_pvalues.csv"))])
second = pd.concat([pd.read_csv(os.path.join(D, "T40_recheck_pvalues.csv")),
                    pd.read_csv(os.path.join(D, "T42_recheck_pvalues.csv"))])
mle = pd.read_csv(os.path.join(D, "T41_mle_pvalues.csv"))

# null cells: T40 design -> (label, p list); T42 cells 201/205/209
designs = {"Dense, correlation 0.7": [1, 8, 15, 22, 29], "Dense, correlation 0.4": [36, 40, 44, 48, 52],
           "Dense, correlation 0.8": [56, 60, 64, 68, 72], "Sparse, correlation 0.7": [76, 80, 84, 88, 92]}
P = [20, 40, 100, 140, 180]
cellmap = {c: (d, P[i]) for d, cs in designs.items() for i, c in enumerate(cs)}
cellmap.update({201: ("Dense, 25% events", 20), 205: ("Dense, 25% events", 100), 209: ("Dense, 25% events", 180)})

def ci(k, n):
    return (beta.ppf(.025, k, n - k + 1) if k > 0 else 0.0, beta.ppf(.975, k + 1, n - k) if k < n else 1.0)

rows = []
for c, (dsg, p) in cellmap.items():
    x = pd.concat([main[main.cell == c], second[second.cell == c]])
    for meth, col in [("SC.HL", "p_schl"), ("SC.EDGE", "p_scedge"), ("HL", "p_hl"),
                      ("GRPtests (5 splits)", "p_grp5"), ("GRPtests (1 split)", "p_grp1"), ("PLStests", "p_pls")]:
        v = x[col].dropna(); k = int((v < .05).sum()); n = len(v); lo, hi = ci(k, n)
        rows.append(dict(design=dsg, p=p, kappa=p / 400, method=meth, rate=k / n, n=n, lo=lo, hi=hi, mle_exists=np.nan))
    # classical HL on the MLE: T41 for T40 cells, the T42 file's own columns for T42 cells
    if c > 200:
        m = x[["mle_exists", "p_hl_mle"]]
    else:
        m = mle[mle.cell == c][["mle_exists", "p_hl_mle"]]
    ex = m.mle_exists.astype(bool); v = m.p_hl_mle[ex].dropna()
    k = int((v < .05).sum()); n = len(v)
    lo, hi = ci(k, n) if n else (np.nan, np.nan)
    rows.append(dict(design=dsg, p=p, kappa=p / 400, method="HL on MLE", rate=k / n if n else np.nan, n=n,
                     lo=lo, hi=hi, mle_exists=ex.mean()))

# BAGofT: dense design under correlation 0.7 only, 100 replicates per cell, p/n 0.05 and 0.25 (Amendments 1, 7)
bag = pd.read_csv(os.path.join(D, "T40_bagoft_pvalues.csv"))
for c, p in ((1, 20), (15, 100)):
    v = bag[(bag.cell == c) & (bag.rep <= 100)].p_bag.dropna(); k = int((v < .05).sum()); n = len(v)
    if n < 100:
        continue
    lo, hi = ci(k, n)
    rows.append(dict(design="Dense, correlation 0.7", p=p, kappa=p / 400, method="BAGofT", rate=k / n, n=n,
                     lo=lo, hi=hi, mle_exists=np.nan))

out = pd.DataFrame(rows).sort_values(["design", "method", "kappa"])
out.to_csv(os.path.join(D, "T40_validity_map.csv"), index=False)
print(out.pivot_table(index=["design", "kappa"], columns="method", values="rate").round(3).to_string())
