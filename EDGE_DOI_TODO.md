# When the EDGE arXiv ID lands — do exactly this

The Series C paper cites the EDGE companion, which is **still on arXiv hold** (arXiv Support
confirmed 2026-08-12: no decision, moderators reminded, no updates given on held submissions).
Everything is already written so that the ID drops into **one place**.

---

## The single edit

`paper_seriesC/paper.tex`, in the bibliography:

```latex
\bibitem[Ebrahim and El-Kotory(2026)]{ebrahim2026edge}
Ebrahim, E.K. and El-Kotory, A. (2026).
EDGE: a closed-form directed test for the calibration of probabilistic binary classifiers.
Submitted to \emph{Advances in Data Analysis and Classification}.
```

Replace the last line with:

```latex
\emph{arXiv preprint} arXiv:<ID>. \doi{10.48550/arXiv.<ID>}.
```

or, if `\doi` is not defined in the class, plainly:

```latex
\emph{arXiv preprint} arXiv:<ID>, doi 10.48550/arXiv.<ID>.
```

**Then recompile twice and check `citeundef=0`. That is the whole job.**

⛔ **Do not write an arXiv number that does not resolve.** Open `arxiv.org/abs/<ID>` and confirm the
title and author list match before pasting. The 2026-07-24 audit on this project found seven wrong
DOIs; this is the same failure mode.

---

## Why nothing else needs to change

**The paper does not depend on the companion for its definitions.** This was deliberate:

- **Eq. (5) of §2.1 defines the EDGE statistic in full** — the orthogonal cubic polynomial basis in
  the group-mean fitted probability, and the quadratic form
  $(Z^{\top}r)^{\top}(Z^{\top}Z)^{-1}(Z^{\top}r)$. A reader can implement it from the paper alone.
- **The CRAN package `ebrahim.gof` is the checkable pointer** — §7 states that the EDGE statistic
  computed by this paper's code agrees exactly with the package on three datasets.
- So `ebrahim2026edge` is cited as **corroboration and provenance**, not as definition. If the hold
  never clears, the paper still stands; the citation just stays as "submitted".

**Where it appears (3 places, all safe):** §2.1 eq. (5) attribution, §7 Software, and the
bibliography entry itself.

---

## If a referee objects before the ID arrives

RSS discourages citing unpublished work. The pre-emptive answer, if you want it in the cover letter:

> The EDGE basis is one of two bases we report, and the paper is self-contained with respect to it:
> equation (5) gives its definition, and the implementation is on CRAN in `ebrahim.gof`. The
> companion manuscript is cited for provenance only, and none of the results depend on it.

That sentence is **not** currently in the cover letter — add it only if you want to raise the issue
first rather than answer it if asked.

---

## Related, and separate

The **Zenodo** side is already done: concept DOI `10.5281/zenodo.21900114`, version DOI
`10.5281/zenodo.21900115`. The paper, supplement and cover letter all cite the **version** DOI,
matching the rest of the family; the README badge points at the concept DOI. Nothing there is
waiting on anything.
