# Three paper lines (research agenda)

The solver finding and the convergence work support three publishable
studies. Full evidence in `ProjectDocs/papers/PAPER-LINES.md` and
`ProjectDocs/DIAG-SOLVE-MIXTO.md`.

## Line 1 — "The staggered fixed point is not the solution"

Solver honesty and the fixed point of operator-split u-σ schemes in FE
codes. Evidence: the measured fixed point (0.2315×), the dishonest Bi-CG
stopping criteria (three false-success exits), the flat-residual identity
for symmetric structures, the Bi-CG ≡ SuperLU byte-identical parity, and
the fixes. Target: computational mechanics journal (CMAME / IJNME, or a
shorter empirical paper).

## Line 2 — Sensitivity of solver formulations on the same system

A systematic comparison on one system: Bi-CG vs SuperLU (byte-identical —
the star result), honest vs dishonest stopping, CG on the SPD matrix,
MINRES, and the effect of the outer staggered scheme. The comparison table
is complete and reproducible.

## Line 3 — Pedagogical: hand calculations catch rubbish FE results

A case study in verification by hand calculation: the one-line statics
check |M_end| + |M_center| = pL²/8 for a fixed-fixed beam, applied to a
seemingly sophisticated 2D plane-strain road base model whose integrated
moments violate it by a factor of ~5. The same model gives exact statics in
the Professional — the detector works, the defect is software-specific.
Lesson: complex models can produce completely rubbish solutions; a
one-line hand calculation detects it. Target: engineering education journal
(IJEE / EJEE).

---

The original Tochnog GNU developers are natural co-authors of these
studies — it is their code, their scheme, and their tests that made the
investigation possible. We would welcome the discussion.
