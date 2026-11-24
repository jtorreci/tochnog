# control_mesh_map

## Description

Maps the quadratic volume elements to linear elements and back, as the
Professional manual 6.210 describes, so that a large linear calculation is
first done with the cheap linear elements and the quadratic calculation
starts from that solution:

- the first `control_mesh_map index -yes` turns every `-hex27` into the
  `-hex8`, every `-tet10` into the `-tet4` and every `-prism15` into the
  `-prism6` of its corner nodes (a `-hex20` too: the GNU turns it into a
  `-hex27` at the start of the run, so it becomes a `-hex8` as well);
- calculate with the linear elements (a direct solver);
- the next `control_mesh_map index -yes` turns them back into their
  quadratic elements. The mid nodes get the linear solution, interpolated
  from the corner nodes (all node dofs: velocities, displacements,
  stresses, ...), and the elements start from these nodal values;
- calculate with the quadratic elements (an iterative solver, for which
  the linear solution is the first guess).

Each further map toggles again. Other elements (`-hex8`, `-tet4`,
`-prism6`, `-prism18`, `-prism12`, 2D `-quad9` / `-tria6`, trusses, ...)
are left as they are. The map runs once, when its own control index is
processed, also when that index has no `control_timestep`, and before the
timestep when the index has one.

The manual advises the option for large linear calculations only, with
`global_element_dof_apply -no`; as in the Professional, the run does not
check it.

## Usage

```
control_mesh_map index switch
```

## Parameters

| Record | Meaning |
|--------|---------|
| `control_mesh_map` (6.210) | `-yes`: map at this control index (to linear when the mesh has quadratic elements, else back). `-no` (default): nothing. |

The `.dbs` holds `element_map index name nodes`, the quadratic element of
a mapped element (= Professional), and the GNU record `element_map_node`
(the coordinates of its nodes).

## Example

The Professional `test/other/map1` (a ring of `-hex27`):

```
global_element_dof_apply -no
...
control_mesh_map            10  -yes
control_timestep            20 1.0 1.0
control_solver              20 -matrix_pardiso
control_mesh_map            40  -yes
control_timestep            50 1.0 1.0
control_solver              50 -matrix_iterative_bicg
```

Self-checking inputs in `validation-suite/test-2026` (a column of two
`-hex27` under self-weight, values = Professional 02-08-2026):
`control_mesh_map1` (final state), `control_mesh_map2` (after the map
back), `control_mesh_map3` (after the linear calculation); meshes
`control_mesh_map_mixed1`, `_toggle1`, `_no1`, `_keep1` (`-prism18`),
`_keep2` (2D `-quad9`), `_tet10_1`, `_prism15_1`, `_hex20_1`; mid node
interpolation `control_mesh_map_tet10_2`, `_prism15_2`.

## Notes

- The quadratic results after the map back are not those of a run with
  only the quadratic elements: the stresses of the linear solution are
  interpolated to the integration points and the quadratic calculation
  adds increments to them. On the coarse column of `control_mesh_map1`
  the displacement of a mid node differs by up to 15 % (-0.01714 against
  -0.01491), on the fine ring of `other/map1` by 2e-5; in the GNU as in
  the Professional (BACKLOG C134; owner decision 2026-11-11: kept, restart the
  quadratic calculation by hand for the plain quadratic answer).

## Differences with the Professional version

- While the elements are linear the GNU deletes the mid nodes (every GNU
  mesh change deletes the nodes no element uses) and creates them again on
  the map back, with the stored coordinates; the Professional keeps them,
  unused, in the `.dbs` (DIVERGENCES MESH-MAP-MID-NODES).
- `-tet10` and `-prism15`: the GNU takes the corner nodes (`-tet4 1 3 6
  10`, `-prism6 1 2 3 4 5 6`); the Professional takes the nodes `1 4 6 10`
  and `1 3 3 4 5 6` (a mid node, a repeated node) (DIVERGENCES
  MESH-MAP-TET10, MESH-MAP-PRISM15).
- A `-hex20` comes back as a `-hex27` (the GNU converts every `-hex20` to
  `-hex27`); the Professional restores the `-hex20` (DIVERGENCES
  MESH-MAP-HEX20).
- The Professional example `other/map1` misses its target in the GNU
  (`-disz` of node 1245: GNU -1.7795e-05, Professional -1.9017e-05, target
  -1.9e-5 +- 1e-6) also without the map: the plain `-hex27` run of that
  ring already differs (BACKLOG C135).
- With `global_element_dof_apply -no` (the manual's advice) the `-hex27`
  column differs from the Professional also without the map (nodal
  stresses up to 1.8 of 20; BACKLOG C136); with the default `-yes` the GNU
  equals the Professional at every stage.
