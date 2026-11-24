# mesh_correct

## Files and functions

- `mesh.cc` — `mesh_correct_interfaces()`: called from `top.cc`
  `step_start()` on every step, right after `interface_convert()`.
- `database.cc` — `MESH_CORRECT` (INTEGER, 1, `no_index`; F5 T4 U5d).
  Enum after `ELEMENT_BOUNDARY` in `tochnog.h` / `tochnog-mod.h`.

## Implementation details

With `MESH_CORRECT -yes`, for every `-quad4` (length 5) and `-hex8`
(length 9) element whose group has a `GROUP_INTERFACE` record (the test
`elem.cc` uses to call `interface_element()`), the distances from node 1
to nodes 2 (xi), 3 (eta) and, for the hex8, 5 (zeta) are taken from the
`NODE` coordinates. The interface direction (eta / zeta) is kept unless
another direction is strictly shorter; then the node permutation measured
on the Professional 02-08-2026 is written back (see the user page) and
`mesh_has_changed()` is called. Idempotent: after the permutation the
interface direction is the shortest.

Verified on the Professional 02-08-2026: quad4 orders 3 4 5 6, 3 5 4 6,
3 4 6 5, 4 3 5 6, 3 6 4 5 and hex8 openings along xi, eta, zeta (all nine
equal).

## Pending

- `-quad9` / `-hex27` interface elements.
- The reversed second side (pairs not facing) is not corrected, as in the
  Professional.
