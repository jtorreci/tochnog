# solver_matrix_symmetric (coupled velocity-pressure systems)

## Description

`solver_matrix_symmetric <switch>` (manual Professional 6.1052):

```
solver_matrix_symmetric -yes
```

With `-yes`, "if needed, matrices are symmetrized so that less memory
will be needed and a symmetrical equation solver can be used" (manual
6.1052). It is the companion of the `bicg` solver family: on systems
that are naturally symmetric (e.g. the bounda_alternate staggered
schemes of manual 6.19) it selects the cheaper plain CG path.

The GNU does NOT symmetrize the assembled matrix: the per-solve
symmetry measurement decides between plain CG and honest Bi-CG. When
the user asserts `-yes` but the assembled system is NOT symmetric (e.g.
a fully coupled consolidation solve where velocities and pore pressures
are solved simultaneously), the GNU prints a warning and runs honest
Bi-CG on the as-assembled system instead of CG — plain CG on a
non-symmetric system diverges:

```
Warning: solver_matrix_symmetric -yes but the measured system is NOT
symmetric - running Bi-CG on the as-assembled matrix instead of CG (the
Professional symmetrizes the matrix in this case)
```

As a consequence the GNU result does not depend on the switch: `-yes`
and `-no` give the same answer (only the solver path and the warning
differ).

## Usage

```
solver_matrix_symmetric -yes
```

Place it in the data part; it is a flat (index-less) record of the
`solver_*` family. See also `solver`, `solver_bicg_error`,
`solver_bicg_stop`.

## Parameters

| Parameter | Meaning |
|---|---|
| `switch` | `-yes`: use the symmetric solver (CG) when the system measures symmetric; Bi-CG on the as-assembled system otherwise (warning). `-no` (default): CG when measured symmetric, Bi-CG otherwise. |

## Example

The guard tests `validation-suite/test-2026/k24_consol_3d1.dat` (with
`-yes`) and `k24_consol_3d2.dat` (with `-no`) solve a 3D Biot
consolidation column (brick 10x2x2 nodes, E = 10, k = 0.1, capacity 1,
unit edge load on x = 0, drained face x = 1) to t = 2.2; both give
pres -0.7083 and disx 0.05421 at node 1 (0, -0.5, -0.5), the
Professional result with `-no`.

## Differences with the Professional version

- **Coupled consolidation (u-p) systems with `-yes`** (DIVERGENCES K24,
  class A, project decision 2026-10-11). The Biot velocity-pressure
  matrix is not symmetric. The Professional 02-08-2026 symmetrizes it
  and solves a different system; the GNU solves the true one. On the
  column above the Professional with `-yes` gives pres -0.678332 at
  t = 2.2 (GNU and Professional `-no`: -0.7083 / -0.7085); on the corpus
  test `large2` it gives -0.672711 at x = 0 (GNU -0.704550, Professional
  `-no` -0.704603, a 4.5 % gap). The Professional `-yes` result is
  non-physical:
  - undrained (k = 1e-8) its pres drifts from -0.80318 (t = 0.1) to
    -0.83404 (t = 5) with no flow; the analytic 1D undrained value is
    -(c + 1/E)/(c + 3/E) = -0.8462 and the GNU and the Professional
    `-no` give -0.8477 (test `k24_consol_3d_undrained1/2`);
  - with capacity 0.1 it gives a positive pres (+0.025, +0.968
    undrained) under compression, and with capacity 0 pres ~0;
  - oedometric column: GNU -0.77306, Professional `-yes` -0.76984, 1D
    finite-element reference -0.7756.
  The storage equation of manual 2.4.1, c dh/dt = k lap(h) + div(v),
  supports the GNU. Inputs with `-no`, or whose system measures
  symmetric, give the same results in both programs.
