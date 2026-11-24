# mpc_geometry options

`mpc_geometry` with `mpc_geometry_method`, `mpc_geometry_switch`,
`mpc_geometry_tolerance` and `mpc_geometry_dof` (manual Professional
6.866-6.870).

## Description

`mpc_geometry index geometry_0 index_0 geometry_1 index_1` generates
`mpc_node_number` records that make dofs of nodes in geometries equal, for
example the right edge of a soil column doing the same as the left edge
(periodic boundary of an earthquake column). Each slave node gets one
record per dof of `mpc_geometry_dof`, with one master and factor 1. The
records are written to the `.dbs` with an `mpc_from_geometry -yes` marker
and are created again when the mesh changes.

| Record | Meaning | Default |
|--------|---------|---------|
| `mpc_geometry_method index method` (6.867) | `-method0`: every node of the second geometry gets the dofs of the node of the first geometry with the same position. `-method1`: give only the first geometry; its first node is the master, all its other nodes are slaves. `-method2`: give only the first geometry; nodes of it with equal coordinates become equal. | `-method0` |
| `mpc_geometry_switch index switch_x switch_y [switch_z]` (6.868) | `-method0`: only the coordinates with `-yes` are compared (in 2D two switches). Required with `-method0`. | - |
| `mpc_geometry_tolerance index tolerance` (6.869) | Coordinates that differ less than `tolerance` count as equal. | `1.e-4` |
| `mpc_geometry_dof index dof_0 dof_1 ...` (6.870) | The dofs that are made equal, e.g. `-velx -vely`. Required. | - |

Details: with `-method0` the master is the first node (in node order) of
the first geometry that matches; a node in both geometries is not tied to
itself. With `-method2` all coordinates are compared and the slaves of a
group of coincident nodes are tied to the highest-numbered node of the
group (for two nodes: the lower node is the slave, as the Professional).

## Usage

```
mpc_geometry        10 -geometry_line 10 -geometry_line 20
mpc_geometry_dof    10 -velx -vely
mpc_geometry_switch 10 -no -yes
```

## Example

The Professional example `mpc2` (copy in `validation-suite/test-2026`):
one quad4 with the left edge in `geometry_line 10` and the right edge in
`geometry_line 20`; the left nodes 1 and 3 are prescribed `disx` 1 and the
right nodes 2 and 4 get the records `2 -velx 1 -velx`, `2 -vely 1 -vely`,
`4 -velx 3 -velx`, `4 -vely 3 -vely`, so they also move 1. Tests:
`mpc2.dat`, `mpc_geometry_method1/2.dat`, `mpc_geometry_tolerance1/2.dat`,
`mpc_geometry_error1/2.dat`.

## Differences with the Professional version

- `mpc_geometry_dof` and (for `-method0`) `mpc_geometry_switch` are
  required, as in the Professional (the manual gives no default).
- `-method2` with three or more coincident nodes: the Professional ties a
  node to each higher coincident node and stops ("nodal dofs cannot be
  slave in more than one multi point constraint"); the GNU ties all of them
  to the highest one.
- The Professional earthquake examples `earthquake_1` ... `_3` use
  `mpc_geometry` but do not run in the GNU yet: their `bounda_time` table
  has 2500 values (BACKLOG C129).
