# Tochnog GNU — Convergence & Solver Work

**Documentation package for the original Tochnog GNU development team**
(fernando lorenzo, osman buyukusik, and all who kept TOCHNOG alive).

Date: 2026-08-28 · Branch: `documentation-improvement`

---

## Who we are

We are a continuation of the Tochnog GNU development. We took the
sourceforge 2014 fork (the last open-source release, maintained by
F. Lorenzo) and have been working on it intensively since August 2026:
implementing features, fixing bugs, and converging the GNU line towards
the features documented in the Tochnog Professional manual.

This package summarises what we have done. Everything is committed in the
repository, with per-feature documentation (user + developer manuals) and
analytical verification.

## The two stories

1. **[The convergence project](01-convergence.md)** — ~200 features from the
   Professional manual implemented into the GNU line, with analytical
   validation, tests, and documentation.

2. **[The solver finding](02-solver-finding.md)** — a deep, reproducible bug
   in the core solver of the open-source line: the staggered u-σ scheme
   converges to the wrong fixed point (section moments 0.23× the statics)
   and the Bi-CG stopping criteria declare success without converging.
   We diagnosed it, fixed it, and validated the fixes against the Tochnog
   Professional binary.

## The offer

- **The code**: the full repository is available; every change is committed
  with documentation. We are happy to share patches, the full history, or
  a tarball.
- **The documentation**: per-feature user/developer manuals, the full
  verification log, the diagnosis document.
- **Collaboration**: the solver finding is publishable (see
  [03-papers.md](03-papers.md) for the three paper lines we are developing).
  The original developers would be natural co-authors — it is your code.
  We would love to discuss it.

## Contact

This package was prepared by the current maintainers of the GNU fork.
Reach us through the repository (or via the shared drive where this
package lives).

---

## Index

| Document | Contents |
|---|---|
| `01-convergence.md` | The convergence project: methodology, scope, verification, current state |
| `02-solver-finding.md` | The solver saga: the bug, the diagnosis, the fixes, the validation against Professional |
| `03-papers.md` | The three paper lines and the evidence behind them |
| (repo) `SEGUIMIENTO-CONVERGENCIA.md` | Full verification log of every implemented feature (Spanish) |
| (repo) `DIAG-SOLVE-MIXTO.md` | The paper-grade solver diagnosis (English, 600+ lines) |
| (repo) `VALIDACION-PROFESIONAL.md` | The GNU ↔ Professional comparison harness and tables |
| (repo) `manual-user/`, `manual-developer/` | Per-feature manuals (English) |
