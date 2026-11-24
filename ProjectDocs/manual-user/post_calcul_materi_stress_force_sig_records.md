# post_calcul -materi_stress -force: the `*_sig` node_dof output records

The `post_calcul -materi_stress -force` family (manual Professional 6.913)
computes per-node section forces and moments of isoparametric elements with a
single element over the structure thickness (sheet piles, tunnel shells,
rings). For every node of the section it produces a set of OUTPUT items stored
in the flat `node_dof_calcul` record. Each item has a NAME, and the names can
be used directly as the value number of a `target_item` record.

## The item names

The plot vectors of the family are drawn in the structure THICKNESS direction
(so that postprocessors show clear vectors). The item names encode the vector
group and the component:

| group | 2D items | 3D items |
|---|---|---|
| normal force | `-norx_sig` `-nory_sig` `-nors_sig` | `-norx_sig` `-nory_sig` `-norz_sig` `-nors_sig` |
| shear force | `-shex_sig` `-shey_sig` `-shes_sig` | `-shex_sig` `-shey_sig` `-shez_sig` `-shes_sig` |
| moment (2D) | `-momx_sig` `-momy_sig` `-moms_sig` | — |
| moment 1 (3D) | — | `-mom1x_sig` `-mom1y_sig` `-mom1z_sig` `-mom1s_sig` |
| moment 2 (3D) | — | `-mom2x_sig` `-mom2y_sig` `-mom2z_sig` `-mom2s_sig` |

The `x`/`y` (and `z` in 3D) components are the GLOBAL plot components: the
vector is drawn in the structure thickness direction, so the components
themselves are NOT the physical components of the force/moment (manual 6.913:
"the components by themselves are not the real physical components of the
force or moment; they are only convenient values for getting clear plots in
postprocessors"). The `s` component is the SIZE of the vector — the real
physical size of the force or moment, usable for design. The `s` items carry
the SIGN of the dominant plot component (compression negative for the normal
force); the shear `s` items are always positive (the manual: "for the shear
force, however, always a positive value is calculated by Tochnog").

## Usage in targets

The names are the value number of a `target_item` against the flat
`node_dof_calcul` record of a node:

```
target_item  10  -node_dof_calcul 7 -nory_sig
target_value 10  -12.34 1.e-2
target_item  20  -node_dof_calcul 7 -shey_sig
target_value 20  +100.  1.e-2
target_item  30  -node_dof_calcul 7 -momy_sig
target_value 30  -5000. 1.e-2
```

The item names are recognized by the parser (they are registered names) and
`exit_tn` resolves them to the slot of the item in the flat `node_dof_calcul`
layout through the generated post_calcul item names. The numeric form
(`target_item N -node_dof_calcul <node> <slot>`) keeps working.

## Output

The same items are written by the prints of the family
(`control_print -node_dof_calcul`, `control_print_materi_stress_force`, VTK
`-node_dof_calcul` blocks) under the names without the leading dash
(`nory_sig`, `moms_sig`, ...).

## Related records

- `post_calcul_materi_stress_force_reference_point` — the reference point that
  orients the plot vectors (the sign of the directional components).
- `post_calcul_materi_stress_force_plot_switch` — inverts the drawing direction
  of an item vector.
- `post_element_force_result` — the section resultants of the
  `post_element_force` family (5 values, no plot components).
