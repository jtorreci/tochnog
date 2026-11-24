# group_truss_memory

## Description

`group_truss_memory index memory_type` (manual Professional 6.777) selects
the kinematics of the truss elements of an element group.

| memory_type | Behaviour |
|---|---|
| `-updated` | geometrically nonlinear: the axis follows the current nodes, the force grows with `EA/L * dL` (current length) |
| `-total` | geometrically nonlinear in the reference configuration: the axis follows the current nodes, `N = N0 + EA (L - L0)/L0` (engineering strain of the initial length) |
| `-updated_linear`, `-total_linear` | geometrically linear: the axis and the stiffness of the initial geometry |

## Default (no record)

The manual (6.777) gives no default. Since 2026-10-11 (F4 T6d, BACKLOG
C11) the GNU follows the Professional 02-08-2026:

- with `materi_displacement` initialised: `-total_linear` (manual 4.12:
  a total Lagrange formulation is used when `materi_displacement` is
  initialised). Pro `dynamic3`: node 2 `disx = 60.44025`; the old GNU
  default `-updated` gave 32.66;
- without `materi_displacement` (velocity-integrated truss):
  `-updated_linear`. A truss whose free end moves by (-0.5, 1) has force
  -50 along its initial axis, as the Professional; `-updated` gives
  10.459.

The linear types use the length projected on the initial axis WITH its
sign, so a truss that is pushed through its other node (inversion) keeps
the correct force (before 2026-10-11 the GNU took the absolute value:
`dynamic3` with `-total_linear` gave node 2 `disx = -5962.6`).

## Accepted types

The manual lists `-updated_linear`, `-updated` and `-total_linear`; the
Professional 02-08-2026 binary also accepts `-total` and treats it as
described above (test `other/truss14`). Since 2026-10-06 the GNU does the
same; before, `-total` behaved as `-total_linear`.

For the nonlinear types (`-updated`, `-total`) the iteration matrix also
holds the geometric stiffness `N/L (I - n n)` of a truss in TENSION. A
straight truss loaded across its axis therefore rotates to the load
instead of stopping with "iterative solver broke down". The converged
state is unchanged by this matrix term.

## Example

```
group_type 0  -truss
group_truss_memory 0  -total
group_truss_area 0  1.0
group_truss_elasti_young 0  1.0
bounda_dof   10  1 -velx -vely
bounda_force 20  2 -velx
bounda_time  20  1.
bounda_force 30  2 -vely
bounda_time  30  1.
```

A unit truss along x with node 1 fixed: force 1 along x (step 20), then
force 1 along y added (step 30). The truss aligns with the resultant:
`N = sqrt(2)`, node 2 `disx = 0.70711`, `disy = 1.70711` (= Professional).
Own tests: `validation-suite/test-2026/truss_total_rotation1.dat` and
`truss_updated_transverse1.dat` (Pro `truss2`); defaults and inversion:
`truss_memory_default_displacement1`, `truss_memory_default_velocity1`,
`truss_total_linear_inversion1`, `truss_dynamic_no_record1`.

## Input checks

As the Professional 02-08-2026 (owner decision 2026-10-11), an explicit
`group_truss_memory` record is checked against `materi_displacement`:

- `-updated` or `-updated_linear` with `materi_displacement` initialised:
  input error "group_truss_memory -updated or -updated_linear is not
  allowed when materi_displacement is initialised";
- `-total` or `-total_linear` without `materi_displacement`: input error
  "group_truss_memory -total or -total_linear needs materi_displacement
  to be initialised".

`check_error -no` skips both checks. The sfnet input `truss5` used
`-updated` with `materi_displacement`; the corpus copy was adapted (record
dropped, default `-total_linear`, target -93.8 from the Professional). The
copy inside `validation-suite/test.zip` keeps the original input and now
stops with the first error. Tests: `truss_memory_updated_displacement_error1`,
`truss_memory_total_velocity_error1`.

## Differences with the Professional version

- None known for the default and the input checks. The GNU messages name
  both records; the Professional prints "materi_displacement should (not)
  be initialised".

## Notes

- A straight truss with ZERO force under a transverse load still has a
  singular first iteration: the GNU prints "iterative solver broke down"
  and continues (BACKLOG C12); after the first update the force builds up
  and the iteration converges (Pro `truss2`/`truss3`/`truss14`).
- In compression the geometric term is left out of the matrix
  (iteration stability); the force itself is exact.
