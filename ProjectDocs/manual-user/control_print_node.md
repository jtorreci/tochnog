# control_print_node

## Description

`control_print_node` prints NODAL data records to plain ASCII files
(manual Professional 6.330). The data item can be `node_dof`,
`node_dof_calcul`, or ANY nodal data record whose name starts with
`node` (e.g. `node_rhside`). The files contain lines like
`x y z <value>` (in 1D only `x`), one line per node. **One file is
written per selected part**.

With `number_0 number_1 ...` you select which parts of the data record
are printed:

- numbers `0, 1, ...`: the value at that position of the record;
- for `node_dof`: dof labels, e.g. `-velx -vely`;
- for `node_dof_calcul`: post_calcul labels (e.g. `-materi_velocity`,
  which selects the average operator `avel`).

The manual example is `control_print_node index -node_dof -velx -vely`
(the `-velx -velx` of the manual is a typo): it produces the files
`velx.index` and `vely.index` with `x y velx` and `x y vely` columns.

The companion records are:

- [`control_print_node_angular`](control_print_node_angular.md) — print an angle instead of the coordinates.
- [`control_print_node_angular_middle`](control_print_node_angular_middle.md) — middle point of the angle axes.
- [`control_print_node_geometry`](control_print_node_geometry.md) — restrict the printing to nodes on a geometry.
- [`control_print_node_sort`](control_print_node_sort.md) — sort the printed lines.
- [`control_print_node_zero`](control_print_node_zero.md) — suppress zero valued results.

## Uso

```
control_print_node 10  -node_dof -velx -vely
control_timestep 10  1.0 1.0
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `10`      | Index of the control record. Must match an active `control_timestep`. |
| `data_item_name` | Any nodal record whose name starts with `node` (e.g. `-node_dof`, `-node_dof_calcul`). |
| `number_0 ...` | Parts to print: numbers or labels. If none is given, ALL parts of the record are printed. |

## Output

- One file per selected part: `label.<index>` for dof labels (e.g.
  `velx.10`, consistent with `control_print_dof_line`) and post_calcul
  items (e.g. `avel.10`); `<record>_<n>.<index>` for numeric parts
  (e.g. `node_dof_0.10`). Each call appends the current state.
- Lines `x y z <value>` (in 1D only `x`). With
  [`control_print_node_angular`](control_print_node_angular.md) the
  first column is the angle in degrees instead of the coordinates.
- Nodes whose record is not active, and parts beyond the record length,
  are omitted.

## Example

A 2D model with `vely` prescribed to `0` at the bottom and `-0.01` at
the top:

```
control_print_node 10  -node_dof -velx -vely
control_timestep 10  1.0 1.0
```

Produces `velx.10` (all zeros) and `vely.10`:

```
0 0 0
1 0 0
0 1 -0.01
1 1 -0.01
```

## Differences with the Professional version

- File naming of the numeric parts is a GNU decision
  (`<record>_<n>.<index>`); the manual only shows the label case.
- `post_calcul_label` does not exist in the GNU: a post_calcul label is
  matched exactly on the underlying unknown (e.g. `-materi_velocity`
  selects the `avel` item) or as a substring of the item label
  (e.g. `-sigyy` matches `asigyy`).
- No parts given means "all parts" (the manual does not state a
  default).
