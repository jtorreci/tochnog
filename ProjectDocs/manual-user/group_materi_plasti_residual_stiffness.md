# group_materi_plasti_residual_stiffness

## Description

Material record (manual Professional 6.735): `factor` is the part of the
original shear stiffness used as stiffness in plastic points. Default 1
(record absent): the full elastic stiffness.

Only the material matrix of the equilibrium iterations changes; the
stresses are still computed by the plasticity law. A converged result is
therefore the same with or without the record. A factor below 1 brings the
matrix of a yielding point closer to its real (softer) stiffness, which can
reduce the equilibrium iterations near failure, or give a better state when
a step ends without convergence. Very small factors make the matrix
ill-conditioned.

A point is plastic when, in the current stress update, its yield test
started plastic iterations or a direct cut (`_mohr_coul_direct`,
`_tension_direct`) changed its stress. The volumetric (bulk) part of the
stiffness is kept; the rest is multiplied by `factor`.

## Syntax

```
group_materi_plasti_residual_stiffness index factor
```

- `index`: element group.
- `factor`: larger than 0 (0 or less is an error).

## Example

```
group_materi_plasti_mohr_coul 0 0.5 1. 0.
group_materi_plasti_residual_stiffness 0 0.1
```

Tests: `validation-suite/test-2026/group_materi_plasti_residual_stiffness1`
(50 iterations: converged state unchanged), `_2` (one iteration: the state
changes), `_3` (no record), `_error1` (factor 0).

## Differences with the Professional version

- The Professional 02-08-2026 accepts the record without any measurable
  effect (identical results and iterations with factor 1, 0.1 and 0.01, also
  in unconverged steps). The GNU follows the manual
  (DIVERGENCES PLASTI-RESIDUAL-STIFFNESS, owner decision 2026-11-06).
- The Professional accepts factor 0; the GNU stops with an error.
- Hypoplastic and other laws without a yield test never count as plastic.
