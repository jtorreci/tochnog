# mpc_geometry options

`mpc_geometry`, `mpc_geometry_method`, `mpc_geometry_switch`,
`mpc_geometry_tolerance`, `mpc_geometry_dof` (manual Professional
6.866-6.870; F5 T4 U7). `mpc_geometry`, `_dof` and `_switch` were
registered without consumption before U7 (the Professional `mpc2` passed
by coincidence, a rigid translation).

## Files and functions

- `mpc.cc`
  - `mpc_geometry_generate()`: called from `mpc_node_apply()` after
    `mpc_element_group_generate()`, every equilibrium iteration; returns
    while `MPC_GEOMETRY_MESH_FINGERPRINT` equals `mpc_mesh_fingerprint()`.
    On a change it deletes the previous range (`mpc_node_number`,
    `mpc_node_factor`, `mpc_from_geometry`) and writes the new records after
    the highest `mpc_node_number`, per active `mpc_geometry` index:
    reads the method (`-method0` default; another word is a `db_error`),
    the tolerance (default 1.e-4, negative is a `db_error`), the dofs (`GET`:
    missing record = the Professional error "Error detected for data item :
    mpc_geometry_dof"; an unknown dof is a `db_error`) and, for `-method0`,
    the switches (`GET`, at least `ndim` values `-yes` / `-no`).
  - `mpc_geometry_nodes()`: the nodes of a geometry in node order
    (`geometry()`, `NODE_START_REFINED`, `PROJECT_EXACT`).
  - The ties are written with `mpc_tie_put()` with `weight==NULL` (one
    master, no factor record = factor 1, as the Professional `.dbs`) and the
    `MPC_FROM_GEOMETRY` marker.
- `database.cc` — `mpc_geometry_method` (INTEGER, 1), `mpc_geometry_tolerance`
  (DOUBLE, 1), both `data_required = MPC_GEOMETRY`; marker
  `mpc_from_geometry`; internal `mpc_geometry_mesh_fingerprint` (fixed 3,
  `external = 0`, no index). The word `-method0` (`METHOD0`, new enum
  before `METHOD1`).
- `tochnog.h` / `tochnog-mod.h` — enums after `MPC_GEOMETRY_SWITCH`.

## Implementation details

- Generated records are user-like for `mpc_node_apply()` (not in the
  `mpc_linear_quadratic` or `mpc_element_group` ranges, no
  `mpc_from_bounda`): a free slave with free masters is condensed in the
  solve, a bounded one keeps its bounda value.
- `-method0`: slave = each node of the second geometry, master = the first
  node of the first geometry (node order) with `|x_i - y_i| <= tolerance`
  on every switched coordinate, not the slave itself.
- `-method1`: master = the first node of the geometry.
- `-method2`: for node k the master is the highest-numbered later node with
  all coordinates within the tolerance (so a coincident group ties to its
  highest node; a pair gives the Professional records).
- Record order: per index, per slave node, the dofs in the order of
  `mpc_geometry_dof` (Professional order).
- Verification: Professional 02-08-2026 example `other/mpc2` and probes
  g0-g10 (method0 switch / tolerance / dof order, method1, method2 pairs,
  missing dof / switch errors). See SEGUIMIENTO-CONVERGENCIA.md, row F5 T4
  U7.

## Pending refactorings

- O(nodes x nodes) matching for `-method0` / `-method2`.
- The hard-coded default tolerance 1.e-4.

## Differences with the Professional version

`-method2` with three or more coincident nodes ties all to the highest node
(the Professional stops with a double-slave error): DIVERGENCES
MPC-GEOMETRY-METHOD2-GROUP (class C). The Professional earthquake examples
are blocked by BACKLOG C129 (`bounda_time` with 2500 values), not by
`mpc_geometry`.
