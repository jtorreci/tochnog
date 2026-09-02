# solver_matrix_symmetric (coupled velocity-pressure systems)

## Description

`solver_matrix_symmetric <switch>` (manual Professional 6.1053):

```
solver_matrix_symmetric -yes
```

With `-yes`, "if needed, matrices are symmetrized so that less memory
will be needed and a symmetrical equation solver can be used" (manual
6.1053). It is the companion of the `bicg` solver family: on systems
that are naturally symmetric (e.g. the bounda_alternate staggered
schemes of manual 6.19) it selects the cheaper plain CG path.

The GNU does NOT re-symmetrize the assembled matrix: the per-solve
symmetry measurement decides between plain CG and honest Bi-CG. When
the user asserts `-yes` but the assembled system is NOT symmetric (e.g.
a fully coupled consolidation solve where velocities and pore pressures
are solved simultaneously), the GNU prints a warning and runs honest
Bi-CG on the as-assembled system instead of CG — plain CG on a
non-symmetric system diverges. This converges where the previous
behaviour (forced CG) diverged; large2 (coupled 3D consolidation under
a surface load) runs to completion.

## Usage

```
solver_matrix_symmetric -yes
```

Place it in the data part; it is a flat (index-less) record of the
`solver_*` family. See also `solver`, `solver_bicg_error`,
`solver_bicg_stop`.

## Note

The GNU result on a measured-asymmetric system is the solution of the
TRUE (as-assembled) system, not of its symmetrization; the two differ
only in the anti-symmetric part of the matrix.
