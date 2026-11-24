# Bi-CG direct-LU retry band width (default solver correctness)

## What it does

When the default iterative solver (Bi-CG) breaks down, Tochnog retries the
same linear system with a direct LU factorization. That retry builds the
system in LAPACK band storage, whose size is set by the maximum distance
between two coupled degrees of freedom. This document describes the fix to
that band width.

## Why it is useful

The band width was underestimated by `nprinc-1` unknowns per coupled node
pair, so for the last-to-first degree-of-freedom coupling of adjacent nodes
(for example the pressure <-> velocity consolidation term of a
`-materi -groundflow` model) the band fill wrote **outside** the allocated
band matrix. The out-of-bounds writes silently corrupted the matrix, and the
direct-LU retry then returned a wrong solution instead of the exact one.

Symptom: a coupled consolidation model (`examp2`, a 1D column with
`inertia_apply -yes` and `groundflow_consolidation_apply -yes`) returned a
pore pressure of `-0.7206` instead of the correct `-0.6345` (the value of
the Professional and of the jan-2014 GNU), with the total stress identical
in both. Forcing the direct solver (`options_solver -matrix_lapack` or
`-matrix_superlu`) already gave the correct value, confirming the retry was
the faulty path.

## Input syntax

No new keyword. The behaviour is automatic:

```
end_initia
inertia_apply -yes
groundflow_consolidation_apply -yes
```

Any model whose default Bi-CG solve breaks down and falls back to the
direct-LU retry benefits; models that converge with Bi-CG are unaffected.

## Parameters

None.

## Minimal example

`external-downloads/sfnet/extracted/test/examp2.dat` with the two records
above after `end_initia`. The post-point pore pressure at t=2.2 is now
`-0.634519741323` (target `-0.63 +/- 0.05`), matching the Professional
`-0.6345197413352` to ten significant digits.

## Notes

- The band width of the Cuthill-McKee ordering path (used by
  `options_solver -matrix_lapack`) already had the correct `(1+diff)*nprinc`
  bound; only the Bi-CG path was short.
- The fix makes the band matrix exact, so the direct-LU retry is now a
  reliable fallback; it is not a change to the physics or to the
  consolidation formulation.
