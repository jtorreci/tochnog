# control_print_gid_other

Manual Professional 6.308.

## Description

With `-no` the "other" results are not written: the mesh deformation (`mesh deform`) and the
`element_group` result of [control_print_gid_group](control_print_gid_group.md). The dof
fields, the post_calcul fields, `coord` and the beam/truss vectors stay. Default `-yes`.

The records apply to the GiD files of `control_print_gid` (`<base>.flavia.msh`, `<base>.flavia.res`). The mesh is written once, at the first GiD print; the results at every print. Overview of the family: [print_gid_options](../manual-developer/print_gid_options.md).

## Syntax

```
control_print_gid_other  index  switch
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Index of a `control_print_gid` record. |
| `switch` | `-yes` (default) or `-no`. |

## Example

```
control_print_gid        10  -yes
control_print_gid_other  10  -no
```

## Differences with the Professional version

- The Professional 02-08-2026 also drops results that the GNU GiD writer never writes
  (truss_force, beam_force_moment_*, bounda_force, node_boundary, *_bounda_dof,
  *_node_rhside*). What both write is dropped alike.

## Tests

`validation-suite/test-2026/`: `control_print_gid_other1` (checked by `scripts/build_safe.sh` in a scratch copy).
