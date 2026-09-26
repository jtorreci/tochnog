# geometry_list

Manual Professional 6.531 and section 9.5. User page:
[../manual-user/geometry_list.md](../manual-user/geometry_list.md).

## Implementation

- `database.cc` (`db_initialize`, next to `GEOMETRY_LINE`): record
  `GEOMETRY_LIST`, type `INTEGER`, `data_length = DATA_ITEM_SIZE`,
  `fixed_length = 0` (variable length), `data_class = GEOMETRY`.
  Enum `GEOMETRY_LIST` inserted after `GEOMETRY_LINE` in both
  `tochnog.h` and `tochnog-mod.h` (same position, headers kept in sync).
- Because the data class is `GEOMETRY`, every generic consumer that accepts
  `-<record> <index>` when `db_data_class(...)==GEOMETRY` accepts
  `-geometry_list` without special-casing: `bounda.cc` (bounda_dof /
  bounda_force, lines ~514 and ~633), `area.cc`, `force.cc`, `filter.cc`,
  `data.cc` (data_put / data_delete), `post.cc`, `print_node.cc`
  (control_print_node_geometry), and `geometry_set` members.
- `geometry.cc` `geometry()`: new branch `entity==-GEOMETRY_LIST` next to
  `-NODE_BOUNDARY`. `in_geometry = 1` iff `inod>=0` and `inod` is one of the
  list values (linear scan of `db_int`/`db_len`). `penetration = 0`,
  `projection = coord` (the node's own coordinates, so
  `control_mesh_adjust_geometry` does not move list nodes), `normal` stays
  zero, `factor` stays 1 (the `geometry_bounda_sine_*` modifiers still
  apply after the branch like for every entity). Coordinate-only calls
  (`inod<0`) are never in the geometry.
- `geometry.cc` `node_geometry_present_calculate()`: `GEOMETRY_LIST` added
  to the scanned entity list between `GEOMETRY_LINE` and `GEOMETRY_POINT`
  (alphabetical, matches the Professional output order).
- `change.cc` `change_geometry()`: explicit "Sorry" exit for a
  `-geometry_list` target, placed before the `db(...GET)` into the
  one-element `idum` buffer (an integer list would overflow it).
- `adjust.cc`: no guard needed (see projection above); only
  `GEOMETRY_SET` keeps its existing refusal.

No external libraries. No hardcoded parameters.

## Verification (Professional binary 02-08-2026)

The only Professional usage example (`test/other/tochnog_abaqus.dat`) is a
non-runnable include fragment, so dedicated tests were built. Scratch runs,
`LC_ALL=C`, clean `build_safe.sh --clean` build.

| test | GNU | Professional |
|------|-----|--------------|
| `geometry_list1` (bounda_dof on lists) | node 4/8 `disx=0.03`, node 8 `temp=10`, nodes 1/5 `temp=0` and `velx=0`, nodes 1-4 `vely=0` | identical prescribed values |
| `geometry_list_delete1` (`control_mesh_delete_geometry` list {3}) | no element deleted | no element deleted |
| `geometry_list_delete2` (list {3 4 7 8}) | element 3 deleted, nodes 4/8 removed | element 3 deleted, nodes 4/8 removed |
| `node_geometry_present` with a list + a point | node 1 `-geometry_list 1 -geometry_point 5`, node 8 `-geometry_list 8` | identical |

Equivalence check: replacing every list with an equivalent
`geometry_line`/`geometry_point` gives a byte-identical `node_dof` field in
BOTH binaries (GNU list = GNU geometric, Pro list = Pro geometric). The
remaining interior differences between GNU and Pro in `geometry_list1`
(for example `sigxx` 9.998 vs 10.0, node 2 `temp` 2.758 vs 2.751) are
therefore a pre-existing numerical difference of the model, not of
`geometry_list`.

Self-checking copies (`target_item`/`target_value`, rc=0 on both binaries;
a perturbed target gives rc=1 on the GNU) live in `validation-suite/test-2026/`
(gitignored). `geometry_list1.dat` is the model shown in the user page with
these targets:

```
target_item 0 -node_dof 4 -disx
target_value 0 0.03 1.e-8
target_item 1 -node_dof 8 -temp
target_value 1 10. 1.e-8
target_item 2 -node_dof 5 -temp
target_value 2 0. 1.e-8
target_item 3 -node_dof 5 -velx
target_value 3 0. 1.e-8
target_item 4 -node_dof 3 -vely
target_value 4 0. 1.e-8
```

The mesh of all three tests: 8 nodes `(0..3, 0..1)`, 3 `-quad4`
(`1 2 5 6`, `2 3 6 7`, `3 4 7 8`), `-materi -condif`, young 1000,
poisson 0, lists `1={1 5}`, `2={1 2 3 4}`, `3={4 8}`, `6={8}`,
`bounda_dof 0 -geometry_list 1 -velx`, `1 -geometry_list 2 -vely`,
`2 -geometry_list 3 -velx` (time 0→0, 1→0.03), `3 -geometry_list 6 -temp`
(time 0→0, 1→10), `4 -geometry_list 1 -temp`, one step of 1.
`geometry_list_delete1.dat` adds `geometry_list 7 3` +
`control_mesh_delete_geometry 0 -geometry_list 7` (targets node 4 `disx=0.03`,
node 8 `temp=10`); `geometry_list_delete2.dat` uses
`geometry_list 7 3 4 7 8` (targets node 7 `disx=0`, `temp=0`).

## Related finding: `-node <n>` as geometry entity (manual 9.4)

Not implemented, not supported: `node` has data class `NODE`, so
`bounda_dof 10 -node 6 -temp` is not taken as a geometry. In the GNU the
record parses and is silently ignored (pure condif strip, all `temp=0`,
rc=0). The Professional 02-08-2026 also parses it but applies nothing
(Pardiso zero pivot on the same model, while the same model with
`-geometry_list` solves). Left as is.

## Pending / notes

- `geometry_set` combined with `bounda_time`: in a probe the Professional
  gave `temp=1` everywhere (ignoring the time table) while the GNU gave 10;
  the same divergence appears with purely geometric set members, so it is a
  `geometry_set` difference, not a `geometry_list` one. Not investigated.
- `change_geometry` on a list: GNU stops with "Sorry"; the Professional
  behaviour was not measured.
- An undefined list index stops with the database error of `db_int`.
- The membership test is a linear scan per node (lists are short records).
