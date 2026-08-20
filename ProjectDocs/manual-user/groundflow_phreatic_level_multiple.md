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
