# groundflow_phreatic_level_multiple (familia)

## Description

The same as `groundflow_phreaticlevel`, but now several groundwater levels can
be specified. For each `groundflow_phreatic_level_multiple` a separate value
for `index` is used.

This option typically can be used if you have in vertical direction
non-permeable layers separating the total domain in independent parts with each
its own groundwater level.

You can specify with one of `groundflow_phreatic_level_multiple_element`,
`groundflow_phreatic_level_multiple_element_group`,
`groundflow_phreatic_level_multiple_element_geometry` or
`groundflow_phreatic_level_multiple_node` the parts of the domain that belong
to the groundwater level of `groundflow_phreatic_level_multiple` with the same
index. Only one of these records can be used, you cannot combine them.

With `groundflow_phreatic_level_multiple_n` you specify `nx ny` in 3D again.
In the `group_type` for elements which should get the static groundflow
pressure you need to add `-groundflow`.

## Usage

```
groundflow_phreatic_level_multiple <index> <water_level ...>
groundflow_phreatic_level_multiple_element_group <index> <element_group_0> ...
groundflow_phreatic_level_multiple_static <index> <switch>
```

## Records

| Record | Meaning                                                     |
|--------|-------------------------------------------------------------|
| `groundflow_phreatic_level_multiple` | Groundwater level (index + value or table).       |
| `groundflow_phreatic_level_multiple_element` | Element numbers for the level with the same index. |
| `groundflow_phreatic_level_multiple_element_group` | Element group numbers for the level with the same index. |
| `groundflow_phreatic_level_multiple_element_geometry` | Element geometry numbers for the level with the same index. |
| `groundflow_phreatic_level_multiple_n` | `nx ny` for the 3D level table.               |
| `groundflow_phreatic_level_multiple_node` | Node numbers for the level with the same index. |
| `groundflow_phreatic_level_multiple_static` | `-yes` sets the total pressure equal to the static pressure in the nodes of the level. |

Only one of `_element`/`_element_group`/`_element_geometry`/`_node` can be used
per index; they cannot be combined.

## Example

```
groundflow_phreatic_level_multiple 0  3.
groundflow_phreatic_level_multiple_element_group 0  0
groundflow_phreatic_level_multiple_static 0  -yes
groundflow_phreatic_level_multiple 1  1.
groundflow_phreatic_level_multiple_element_group 1  1
groundflow_phreatic_level_multiple_static 1  -yes
```

Two independent groundwater levels: group 0 under a level at height 3 and group
1 under a level at height 1. With `_static -yes` the pore pressure in the nodes
of each level is set to the static pressure `dens*g*(water_level - coord)`
(e.g. -3 and -1 at the bottom with density 1 and gravity (0,-1)), without
solving the hydraulic heads.

## Behaviour notes (2026-09-07, ground8 verification)

- The level records may use ANY index (e.g. 10/20); the family is not
  restricted to index 0.
- A level WITHOUT `groundflow_phreatic_level_multiple_static -yes` applies the
  free-surface condition of a single `groundflow_phreaticlevel` to the nodes of
  its own domain: nodes of the domain at or above the phreatic line are dry and
  get their hydraulic-pressure dof bounded to 0 (total pressure 0 there), which
  confines the saturated flow domain below the level. Below the level the total
  pressure evaluates to the hydrostatic profile `dens*g*(water_level - coord)`
  (capped to the atmospheric pressure above the level), matching the
  Professional .dbs of ground8 of the corpus (hydrostatic under each level,
  0 elsewhere and in the dry zones).
- Post points inside a level domain resolve the level of the element that
  contains the point, so `post_point_dof_calcul ... -to_pres` etc. reproduce
  the same profile as the nodes.

## Verified with

- ground8 of the corpus (two levels at y=-50 and y=-10 over two element
  groups, mechanics + flow): post-point total pressures -10 at both targets,
  rc=0; profile digit-consistent with the Professional binary .dbs
  (25-10-2023). See the SEGUIMIENTO-CONVERGENCIA.md register.
