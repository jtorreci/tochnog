# matrix_pardiso (solver type `-matrix_pardiso`)

## What it does

`-matrix_pardiso` asks for the PARDISO direct sparse solver, the default
solver of Tochnog Professional (manual 6.1047, `solver` /
`control_solver`). PARDISO is a commercial library and is not compiled in
the GNU. Since F4 T9 (2026-10-10) the GNU solves `-matrix_pardiso` with
its own direct sparse solver, SuperLU, and falls back to the iterative
Bi-CG when the SuperLU factorization fails. Before, `-matrix_pardiso`
always ran the Bi-CG.

## Why it is useful

Inputs written for the Professional keep their direct-solver behaviour.
Badly conditioned 3D models that the Bi-CG cannot solve in reasonable
time finish: the Professional `force11` (one 3D solve) needed more than
an hour with the Bi-CG and takes about 10 s with SuperLU.

The GNU default solver stays the Bi-CG (`-matrix_iterative_bicg`): a
direct factorization does not scale to large 3D models (the Professional
`large1`, 30 x 30 x 30 hex8, needs more than 4 GB with SuperLU, while the
Bi-CG solves it in 0.65 GB).

## Syntax

```
solver -matrix_pardiso
```

or per control index:

```
control_solver 10 -matrix_pardiso
```

## Behaviour

- The first time it prints `Warning: -matrix_pardiso is not compiled in
  the GNU - solving with the SuperLU direct solver (Bi-CG if SuperLU
  fails)`.
- SuperLU runs under a memory budget: 70 % of the room left in the
  process (the address-space limit, e.g. `ulimit -v`, or the available
  memory, whichever is smaller). If the factorization runs out of that
  budget, or the matrix is exactly singular, the solve is repeated with
  the Bi-CG (and its direct band LU retry) and the message `Warning:
  SuperLU failed for -matrix_pardiso - falling back to the iterative
  Bi-CG solver (for the rest of the calculation)` is printed; every later
  solve of the run uses the Bi-CG.
- The `solver_pardiso_*` records are still accepted and ignored.
- The environment variable `TN_PARDISO_SUPERLU_MB` replaces the budget by
  the given number of megabytes (diagnostics and tests).

## Example

```
echo -no
number_of_space_dimensions 3
materi_velocity
materi_displacement
materi_stress
end_initia
solver -matrix_pardiso
...
```

Tests: `validation-suite/test-2026/pardiso_superlu1.dat` (SuperLU solves
it) and `pardiso_bicg_fallback1.dat` (run with `TN_PARDISO_SUPERLU_MB=4`:
SuperLU runs out of memory and the Bi-CG gives the same analytic result).

## Differences with the Professional version

- **Default solver** (class B, deliberate, decided 2026-10-10). GNU:
  Bi-CG/CG iterative solver. Professional 02-08-2026: PARDISO direct
  solver. A direct solver does not scale to large 3D models in the GNU
  (SuperLU on `large1` exceeds 4 GB). Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **`-matrix_pardiso`** is solved with SuperLU (class B): direct solve,
  same solution up to round-off; the PARDISO options
  (`solver_pardiso_*`) have no effect.
