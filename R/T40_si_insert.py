"""Insert (or replace) the SI section on the p/n grid in bimj/src/supplement_bimj_v3.tex, with the tables that
T40_si_tables.py generated from the data. Idempotent: an existing section with the same label is replaced."""
import io, os, re, subprocess, sys
sys.path.insert(0, r"C:\Users\ebrah\.gemini\Projects")
from safe_write import safe_write

HERE = os.path.dirname(os.path.abspath(__file__))
subprocess.run([sys.executable, os.path.join(HERE, "T40_si_tables.py")], check=True)
tables = io.open(os.path.join(HERE, "..", "data", "T40_si_tables.tex"), encoding="utf-8").read()
SP = r"C:\Users\ebrah\.gemini\Projects\PDFs\arashi_proposal\paper_seriesC\bimj\src\supplement_bimj_v3.tex"
s = io.open(SP, encoding="utf-8").read()

section = r"""\section{Behaviour across $p/n$ and comparison with packaged tests}\label{sec:t40}

% T40, T41, T42 and the second blocks; protocol PREDECLARATION_T40_highdim.md with Amendments 1-6.
\paragraph{Protocol.} The grid of Section~\mref{sec:grid} of the main text was declared, with its predictions,
before any of it was run, and its file hashes were recorded; six dated amendments followed, each recorded
before the data it concerns. In order: (1) BAGofT \citep{zhang2023} runs last, on $p/n = 0.05$, $0.25$ and
$0.45$ only, at $100$ replicates, because one test takes $37$ to $111$ minutes; (2) the classical test on the
maximum likelihood fit is added on the same datasets; (3) BAGofT runs one R process per test after a socket
failure in the timing pilot; (4) the design with an intercept ($25\%$ events), added after an interim look
raised the concern that a truth without an intercept favours PLStests, and reported as confirmatory of that
concern; (5) an independent second block for \emph{every} correct-model cell at $p/n \ge 0.25$, so that the
cells that looked worst are not the only ones re-run; and (6) the same for the $p/n = 0.25$ cell of the
intercept design; and (7) the submitted paper reports BAGofT at $p/n = 0.05$ and $0.25$, the range over
which the corrected test is claimed valid, because one test took a median of $42$ and $86$ minutes there;
its $p/n = 0.45$ cells were still running at submission and will be reported in full.

\paragraph{BAGofT.} BAGofT was run with \texttt{testGlmnet(alpha = 0)}, a ridge fit tuned by
cross-validation, at its package defaults ($100$ splits, $100$ bootstrap draws), one R process per test, on
the first $100$ datasets of the null, along-risk and one-covariate cells of the dense design under
correlation $0.7$; its $p$-value is the package's \texttt{p.value}. Every row's data fingerprint matches the
main run. On correct models it rejected $0.04$ at $p/n = 0.05$ and $0.06$ at $0.25$ (exact interval
$[0.022, 0.126]$). Against misfit along the predicted risk it rejected $0.07$ and $0.00$, where
\texttt{SC.EDGE} on the same datasets rejected $0.79$ and $0.28$; against misfit in one covariate it
rejected $0.59$ and $0.43$, close to GRPtests ($0.63$ and $0.24$ with five splits). The median time per test
was $42$ minutes at $p = 20$ and $86$ minutes at $p = 100$, against $3$ and $13$ seconds for the corrected
test on the same datasets.

\paragraph{Designs.} $n = 400$ and $p \in \{20, 40, 100, 140, 180\}$. Covariates are Gaussian with
autoregressive correlation $\rho$; the dense signal has coefficients cycling through $0.35, -0.30, 0.25,
-0.20, 0.15$, the sparse one five equal non-zero coefficients, each scaled so that the true linear predictor
has standard deviation $1.5$. Departures at strength $a$ add to the linear predictor a standardized quadratic
in the true index (along the predicted risk), $(x_1^2 - 1)/\sqrt{2}$ (one covariate), or a standardized
$x_1x_2$ (interaction). The penalty is chosen by the one-standard-error rule of ten-fold cross-validation in
every dataset and held fixed in the bootstrap ($499$ draws), with $G = 10$. GRPtests and PLStests run at their
package defaults (GRPtests also with one split). Replicate $r$ of cell $c$ uses seed $10^{6}c + r$, and every
row stores a fingerprint of its data, which was checked across the files that share datasets.

\paragraph{Findings beyond the main text.} GRPtests with five splits returns $p = 1$ on $32\%$ of
correct-model datasets, and its one- and five-split versions reach opposite verdicts on the same data in
$11.8\%$ of them. PLStests returns $p < 10^{-4}$ on $1.9\%$ of correct models without an intercept and on all of
them with one. The bootstrap generator of the corrected test pins a probability at zero or one in fewer than
$0.3\%$ of datasets, so that is not why the test becomes liberal at large $p/n$; the spread of its linear
predictor exceeds the truth's by a factor that grows from $0.96$ at $p/n = 0.05$ to $1.10$--$1.21$ at
$p/n = 0.45$, largest for the sparse signal, the order in which the size fails.

\paragraph{An oracle check of the generator.} To separate the strength of the generator's signal from its
direction, we re-ran three correct-model cells on the same datasets ($60$ replicates each, $199$ bootstrap
draws) with generators that a user could not build, because they use the true signal strength: the
debiased direction rescaled to that strength, and the ridge direction rescaled to it, each with the
intercept re-solved to match the observed events. In the sparse design at $p/n = 0.45$ the decile-basis
rejection rate of correct models was $0.233$ with the usual generator on these replicates, $0.167$ with the
debiased direction at the true strength and $0.050$ with the ridge direction; in the $25\%$-event design at
$p/n = 0.45$ it was $0.267$, $0.117$ and $0.000$; and in the dense design at $p/n = 0.25$, where the test is
valid, $0.033$, $0.033$ and $0.000$. With $60$ replicates the standard error near $0.05$ is about $0.03$, so
this is a diagnostic, not an estimate: correcting the strength alone roughly halves the excess, and the
smoother direction removes it at the price of a conservative test.

""" + tables + "\n\n"

s = re.sub(r"\\section\{Behaviour across \$p/n\$ and comparison with packaged tests\}.*?(?=\\section\{)", "", s, flags=re.S)
anchor = "\\section{Reproducibility}"
assert s.count(anchor) == 1
s = s.replace(anchor, section + anchor)
if "{zhang2023}" not in s.split("\\begin{thebibliography}")[1]:
    item = ("\\bibitem[Zhang et~al.(2023)]{zhang2023}\nZhang, J., Ding, J. and Yang, Y. (2023).\n"
            "Is a classification procedure good enough? A goodness-of-fit assessment tool for classification\nlearning.\n"
            "\\emph{Journal of the American Statistical Association}, \\textbf{118}(542), 1115--1125.\n"
            "\\url{https://doi.org/10.1080/01621459.2021.1979010}.\n\n")
    s = s.replace("\\end{thebibliography}", item + "\\end{thebibliography}")
safe_write(SP, s, keep_backup=True)
print("SI section inserted")
