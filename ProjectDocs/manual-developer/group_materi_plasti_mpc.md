# group_materi_plasti_mpc

## Files and functions

- `database.cc`: `GROUP_MATERI_PLASTI_MPC` (INTEGER switch) and
  `GROUP_MATERI_PLASTI_MPC_FACTOR` (DOUBLE), both length 1, class MATERI,
  requiring `GROUP_TYPE` (F5 T3 C2). Enums appended in `tochnog.h` /
  `tochnog-mod.h`.
- `group.cc`, `group_materi_plasti_boundary_evaluate()`: after the
  bounda / element-group detection, when no wall was found and the group
  has `-yes`, scans every active `MPC_NODE_NUMBER` record; a positive entry
  that is a node of the element sets `plasti_on_boundary = 2`.
- `plasti.cc`, `plasti_rule()`: `plasti_on_boundary==2` reads
  `GROUP_MATERI_PLASTI_MPC_FACTOR` instead of
  `GROUP_MATERI_PLASTI_BOUNDARY_FACTOR` (default
  `BOUNDARY_REDUCTION_FACTOR` 2/3).
- `stress.cc`, `materi_direct_cutoff()`: same selection for the plane
  direct cut-off (default 1 there, pre-existing; BACKLOG C114).

## Implementation details

- The value 1 (bounda / element group) wins over 2 because the mpc scan
  runs only when no wall was found (Professional 02-08-2026: the bounda
  factor wins).
- The scan is O(number of mpc records) per element evaluation and runs only
  for groups with `-yes`.
- The other wall consumers (`_wall` records, camclay, hypo, tension) test
  `plasti_on_boundary` as a boolean, so 2 behaves as 1 there.

## Verification

Professional 02-08-2026 probes (t3/c2): two quad4, mpc on the nodes of
element 2 only: element 1 unchanged in both codes; `-no` and a factor
without `-yes` have no effect; bounda + mpc = bounda alone. The reduced
values differ by the missing cohesion reduction (BACKLOG C114).

## Pending

- Cohesion reduction and the direct-path defaults: BACKLOG C114.
