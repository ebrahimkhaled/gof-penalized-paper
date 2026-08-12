# Cover letter — JRSS Series C

*Paste into the submission system's cover-letter box, or send as a PDF.*
*Rewritten 2026-08-12 to match the current manuscript — the previous version carried the
old title and three claims the paper has since corrected.*

---

To the Editors,
*Journal of the Royal Statistical Society, Series C (Applied Statistics)*

Dear Editors,

Please consider the enclosed manuscript, **"Shrinkage invalidates the Hosmer–Lemeshow
test: goodness of fit for penalized logistic regression, with an application to glaucoma
diagnosis"**, for publication in Series C.

**The applied problem.** A diagnostic model for glaucoma is built from 62 morphometric
measurements of the optic nerve head on 196 eyes. The measurements are repeated geometric
summaries of the same disc, so they are collinear by construction — the largest pairwise
correlation is 0.996 — and the maximum likelihood estimate does not exist. The analyst
must penalize; there is no alternative. Having fitted a ridge model that discriminates
well (AUC 0.905), the analyst applies the check that clinical prediction papers use almost
universally, and obtains *p* < 10⁻⁵ at the cross-validated penalty, with rejection at the
5% level in all twelve penalty-and-grouping configurations. Read at face value, that says
the model is badly misspecified and should be discarded.

**My contribution is that the reading is not available.** The test's reference
distribution is derived under maximum likelihood estimation — precisely what these data do
not admit. Under a penalty the grouped standardized residuals acquire a non-centrality
induced by shrinkage, and in simulation designs of comparable dimension the uncorrected
test rejects *correctly specified* models between 92% and 100% of the time. I derive the
corrected law, show that subtracting an estimate of the shrinkage non-centrality restores
the ordinary maximum likelihood reference exactly to first order, and obtain a valid
procedure by prepivoting. On the glaucoma data the corrected analysis returns *p* ≈ 0.03
rather than *p* < 10⁻⁵: the evidence for misfit is real but modest, and weaker than the
standard analysis suggested by more than three orders of magnitude.

**Why Series C.** The paper is motivated by, and returns to, real clinical data, and its
conclusion is about what a clinician should do differently. The residual misfit the
corrected test finds has a recognizable shape — the predicted probabilities are too flat,
not mis-ordered — so the model remains usable for the risk ranking a screening triage
requires, and the remedy where probabilities are read as numbers is recalibration rather
than reconstruction. A second dataset from the same clinical question, `GlaucomaMVF`,
gives a *different* verdict under the same mechanism: there the correction reverses the
finding completely. Two outcomes from one mechanism is, I think, the right advertisement
for a correction that does not simply manufacture acceptance. Following Series C's
guidance I have placed the proofs and the full simulation grids in web-based supporting
materials.

**Something a reader can use immediately.** Section 2.5 gives a screen that costs one fit
and no refits — roughly 150 times cheaper than the test itself — which tells an analyst
whether their own fit is in the affected regime. It compares two quantities already
computed inside the correction, and across twelve configurations its rule of thumb is
right in eleven.

**On scope and honesty of reporting.** I am explicit about the limits. The theory is
first order, treats the penalty as non-random with λ = o(n), and holds the dimension
fixed; validity at *p/n* = 0.25, and on the glaucoma data at 0.32, rests on prepivoting
and on simulation rather than on the asymptotics, and I say so. I give the attenuation
law governing what any grouped test can detect once the linear predictor must be
estimated, and present it as a contribution rather than an apology. I report that the
correction **costs power** — 2.5 percentage points on the decile basis, *p* = 0.0007 —
which I could have omitted; I also show that choosing a better basis recovers about four
and a half times that. I report the four neighbouring constructions that fail, with the
size of each. I report one simulation cell that is detectably conservative, together with
the fact that my first estimate of it did not replicate and that the pooled estimate over
4000 replications is the one I quote. I report a head-to-head comparison in which a
competing test **beats mine decisively** on one class of alternative, and show that the
two tests reject on disjoint datasets. And I report a dataset I rejected as an example
because cross-validation had selected an essentially intercept-only model there, which
would have flattered the method.

The work is original and is not under consideration elsewhere. I declare no conflict of interest and received no funding. The
`GlaucomaM` and `GlaucomaMVF` data are publicly available in the R packages `TH.data` and
`ipred`; all code, seeds and intermediate results are available at
https://github.com/ebrahimkhaled/gof-penalized-paper and deposited at
https://doi.org/10.5281/zenodo.21900115. The EDGE basis is implemented in the R package
`ebrahim.gof` on CRAN.

I would be glad to suggest referees if that would help.

Yours sincerely,

**Ebrahim Khaled Ebrahim** (corresponding author)
Department of Applied Statistics, Faculty of Business, Alexandria University, Egypt
ebrahimkhaled@alexu.edu.eg · ORCID 0009-0006-7839-8778


---

## Checklist before sending

- [x] Zenodo DOI 10.5281/zenodo.21900115 in the letter, the paper and the supplement
- [x] Title matches the manuscript
- [ ] Upload: `paper.pdf`, `supplement.pdf`, source `.tex` + `Fig/` + `oup-authoring-template.cls` + `oup-abbrvnat.bst`
- [ ] Alt text is already in the source via `\figalttext` — check the system does not ask for it separately
- [ ] Suggested referees, if invited to name any
- [x] `-ize` / `-ise` DECIDED: **keep `-ize`.** "Penalized" is valid Oxford/UK usage per *Chambers*,
      it is used consistently (checked: no stray `-ise` forms), and it is in the **title** — a
      copy-editor conversion to "penalised" after acceptance would change the searchable string
      irreversibly. If you want production to leave it alone, add one line to the letter:
      *"Spellings follow Chambers; I have used the -ize forms deliberately and would ask that they
      be retained, since the title string is already indexed."*
- [x] Software claim checked against the installed package: `ebrahim.gof` 2.4.0 exports
      `edge.gof()` but **no `gof.pen()`** — §7 corrected accordingly
