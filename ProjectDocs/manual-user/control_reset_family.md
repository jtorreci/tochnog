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
  accumulated interface data (strains + stresses + tangential forces) of
  the interface elements in the geometry.
- `control_reset_interface_strain index -<geometry_entity> <index>` —
  reset the accumulated strains to 0 but REMEMBER the stresses (manual
  Professional 6.355: "the interface stresses at this moment of
  resetting will be remembered... the new interface stresses are
  calculated from the interface stresses at this moment of resetting
  plus stress due to additional deformation"). The accumulated normal
  stress lives in its own history (`element_interface_force_norm`,
  internal), so the reset zeroes the strains only: a constant load does
  NOT re-compress the interface after the reset (interface10 of the
  corpus: displacement −6e-4 and sigma_n −6 stay, strain record ≈ 0 —
  identical to the Professional).
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
