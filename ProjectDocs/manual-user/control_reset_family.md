# control_reset family (Sprint 8 additions)

## Description

Region/node restrictions for `control_reset_dof` (manual Professional
6.351-6.356), applied when at least one of them is present (without
them, all nodes are treated — previous GNU behaviour):

- `control_reset_geometry index -<geometry_entity> <index>` — reset only
  on nodes of elements COMPLETELY inside the geometry.
- `control_reset_node index node_0 node_1 ...` — reset only on nodes of
  elements with ALL their nodes in the list.
- `control_reset_element_group index group_0 group_1 ...` — restrict the
  elements considered to those of the listed element groups (combines
  with the two above).

Interface resets (independent records):

- `control_reset_interface index -<geometry_entity> <index>` — reset ALL
  accumulated interface data (normal strain + tangential forces) of the
  interface elements in the geometry.
- `control_reset_interface_strain index -<geometry_entity> <index>` —
  reset ONLY the normal strain; the tangential force history (the
  "remembered" stresses) is kept: new strains start at 0 and new
  stresses grow from the remembered ones through the stiffness.
- `control_reset_element_dof index <-yes|-no>` — registered; -yes (only
  element_dof/element_intpnt_dof, not node_dof) is not yet wired
  (documented partial).

Typical use (manual): zero displacements/strains after gravity in an
updated formulation, or restart interface memory after excavation.

## Example

```
control_reset_dof            1  -hisv0
control_reset_value_constant 1  0.
control_reset_geometry       1  -geometry_brick 1
```

Test `creset_geom`: hisv0 0.5 -> 0 only in the left column (inside the
brick); the right column keeps 0.5 — without the filter both reset
(the A/B fails). Test `creset_iface`: the interface normal-strain
history reads exactly 0.0 after `control_reset_interface_strain`.
