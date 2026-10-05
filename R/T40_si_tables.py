"""LaTeX tables for the Supporting Information section on the p/n grid (T40-T42), generated from the data so
no value is typed by hand. Writes data/T40_si_tables.tex: (1) size of every method in every correct-model cell,
(2) the two independent blocks of every re-run cell, (3) size-adjusted power at the stronger departure."""
import os
import numpy as np, pandas as pd

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
v = pd.read_csv(os.path.join(D, "T40_validity_map.csv"))
DES = ["Dense, correlation 0.4", "Dense, correlation 0.7", "Dense, correlation 0.8",
       "Sparse, correlation 0.7", "Dense, 25% events"]
SHORT = dict(zip(DES, ["dense 0.4", "dense 0.7", "dense 0.8", "sparse 0.7", "25\\% events"]))
METH = ["SC.HL", "SC.EDGE", "HL", "HL on MLE", "GRPtests (5 splits)", "GRPtests (1 split)", "PLStests"]
HEAD = ["\\texttt{SC.HL}", "\\texttt{SC.EDGE}", "HL", "HL (MLE)", "GRP5", "GRP1", "PLS"]

def f(x):
    return "---" if pd.isna(x) else "%.3f" % x

out = []
out.append("\\begin{table}[htbp]\n\\centering\n\\caption{Rejection rate of correctly specified models at nominal "
           "level $0.05$, $n = 400$, over the grid of Figure~\\mref{fig:validity}. HL (MLE) is the classical test on "
           "the maximum likelihood fit, over the datasets in which that estimate exists; the last column gives that "
           "share. Cells at $p/n \\ge 0.25$ pool two independent blocks of $1000$ replicates "
           "(Table~\\ref{tab:t40blocks}); the others have $1000$.}\n\\label{tab:t40size}\n\\small\n"
           "\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{llrrrrrrrr}\n\\toprule\n"
           "Design & $p/n$ & " + " & ".join(HEAD) + " & MLE exists \\\\\n\\midrule")
for i, dsg in enumerate(DES):
    x = v[v.design == dsg]
    for k in sorted(x.kappa.unique()):
        y = x[x.kappa == k].set_index("method")
        me = y.loc["HL on MLE", "mle_exists"]
        out.append("%s & %.2f & %s & %.2f \\\\" % (SHORT[dsg], k, " & ".join(f(y.loc[m, "rate"]) for m in METH), me))
    if i < len(DES) - 1:
        out.append("\\addlinespace")
out.append("\\bottomrule\n\\end{tabular}\n\\end{table}\n")

# second blocks
main = pd.concat([pd.read_csv(os.path.join(D, "T40_main_pvalues.csv")), pd.read_csv(os.path.join(D, "T42_intercept_pvalues.csv"))])
sec = pd.concat([pd.read_csv(os.path.join(D, "T40_recheck_pvalues.csv")), pd.read_csv(os.path.join(D, "T42_recheck_pvalues.csv"))])
lab = {15: ("dense 0.7", .25), 22: ("dense 0.7", .35), 29: ("dense 0.7", .45), 44: ("dense 0.4", .25),
       48: ("dense 0.4", .35), 52: ("dense 0.4", .45), 64: ("dense 0.8", .25), 68: ("dense 0.8", .35),
       72: ("dense 0.8", .45), 84: ("sparse 0.7", .25), 88: ("sparse 0.7", .35), 92: ("sparse 0.7", .45),
       205: ("25\\% events", .25)}
out.append("\\begin{table}[htbp]\n\\centering\n\\caption{The two independent blocks of every re-run correct-model "
           "cell (replicates $1$--$1000$ and $1001$--$2000$, seeds disjoint), declared before the second block was "
           "run. The blocks differ by less than two standard errors of their difference in every cell except the "
           "dense design under correlation $0.8$ at $p/n = 0.35$ on the decile basis, where the first block was the "
           "higher by $0.027$.}\n\\label{tab:t40blocks}\n\\small\n"
           "\\begin{tabular}{llrrrr}\n\\toprule\n& & \\multicolumn{2}{c}{\\texttt{SC.HL}} & "
           "\\multicolumn{2}{c}{\\texttt{SC.EDGE}} \\\\\n\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\n"
           "Design & $p/n$ & block 1 & block 2 & block 1 & block 2 \\\\\n\\midrule")
for c, (dsg, k) in lab.items():
    a, b = main[main.cell == c], sec[sec.cell == c]
    out.append("%s & %.2f & %.3f & %.3f & %.3f & %.3f \\\\" % (dsg, k, (a.p_schl < .05).mean(), (b.p_schl < .05).mean(),
                                                               (a.p_scedge < .05).mean(), (b.p_scedge < .05).mean()))
out.append("\\bottomrule\n\\end{tabular}\n\\end{table}\n")

# size-adjusted power (T40 from T40_size_adjusted.csv; T42 computed here with the same convention)
s = pd.read_csv(os.path.join(D, "T40_size_adjusted.csv")); s = s[s.a == 1].copy()
s["dsg"] = [("sparse 0.7" if g == "sparse" else "dense %.1f" % r) for g, r in zip(s.signal, s.rho)]
t = pd.read_csv(os.path.join(D, "T42_intercept_pvalues.csv"))
t["p"] = np.array([20, 100, 180])[(t.cell - 201) // 4]; t["dep"] = np.array(["null", "index", "coord", "inter"])[(t.cell - 201) % 4]
M = {"SC.HL": "p_schl", "SC.EDGE": "p_scedge", "GRP5": "p_grp5", "GRP1": "p_grp1", "PLS": "p_pls"}
rows = []
for p in [20, 100, 180]:
    nul = t[(t.p == p) & (t.dep == "null")]
    for dep in ["index", "coord", "inter"]:
        g = t[(t.p == p) & (t.dep == dep)]
        r = {"dsg": "25\\% events", "kappa": p / 400, "dep": dep}
        for m, c in M.items():
            thr = np.quantile(nul[c], .05, method="inverted_cdf"); r[m] = (g[c] <= thr).mean()
        rows.append(r)
s = pd.concat([s[["dsg", "kappa", "dep"] + list(M)], pd.DataFrame(rows)])
DEPN = {"index": "along risk", "coord": "one covariate", "inter": "interaction"}
out.append("\\begin{table}[htbp]\n\\centering\n\\caption{Size-adjusted power at departure strength $a = 1$ "
           "($500$ replicates per cell): each test is thresholded at the fifth percentile of its own $p$-values in "
           "the matching correct-model cell. In the design with $25\\%$ events PLStests rejects every correct model, "
           "so its $p$-values are all near zero and its size-adjusted power carries no information there "
           "(shown as ---).}\n\\label{tab:t40power}\n\\small\n\\setlength{\\tabcolsep}{4pt}\n"
           "\\begin{tabular}{lllrrrrr}\n\\toprule\nDesign & Departure & $p/n$ & \\texttt{SC.HL} & \\texttt{SC.EDGE} & "
           "GRP5 & GRP1 & PLS \\\\\n\\midrule")
order = ["dense 0.4", "dense 0.7", "dense 0.8", "sparse 0.7", "25\\% events"]
first = True
for dsg in order:
    for dep in ["index", "coord", "inter"]:
        x = s[(s.dsg == dsg) & (s.dep == dep)].sort_values("kappa")
        if not len(x):
            continue
        if not first:
            out.append("\\addlinespace")
        first = False
        for _, r in x.iterrows():
            pls = "---" if dsg.startswith("25") else "%.3f" % r["PLS"]
            out.append("%s & %s & %.2f & %.3f & %.3f & %.3f & %.3f & %s \\\\" % (
                dsg, DEPN[dep], r.kappa, r["SC.HL"], r["SC.EDGE"], r["GRP5"], r["GRP1"], pls))
out.append("\\bottomrule\n\\end{tabular}\n\\end{table}\n")

with open(os.path.join(D, "T40_si_tables.tex"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(out))
print("written", sum(1 for l in out if l.endswith("\\\\")), "rows")
