# control_print_gid_group

Manual Professional 6.304 and print_gid_group 6.984.

## Description

With `-yes` the element group of every element is written as the element result
`element_group`, on the Gauss points of each element type (`GaussPoints "element_<type>"`: one
internal point for the isoparametric elements, two points with the nodes included for truss,
beam, truss_beam and spring2). In GiD a contour fill shows the groups.
[control_print_gid_other](control_print_gid_other.md) `-no` drops it.

The global `print_gid_group` is the default of every `control_print_gid_group`.

The records apply to the GiD files of `control_print_gid` (`<base>.flavia.msh`, `<base>.flavia.res`). The mesh is written once, at the first GiD print; the results at every print. Overview of the family: [print_gid_options](../manual-developer/print_gid_options.md).

## Syntax

```
control_print_gid_group  index  switch
print_gid_group  switch
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Index of a `control_print_gid` record. |
| `switch` | `-yes`: write `element_group`; `-no`: do not. |

## Example

```
control_print_gid        10  -yes
control_print_gid_group  10  -yes
```

## Differences with the Professional version

- Without either record the GNU writes no `element_group` result (GiD files unchanged); the
  manual and the Professional 02-08-2026 write it by default.
- The GNU writes one internal Gauss point for the isoparametric elements (the Professional writes
  the element integration points with given natural coordinates); the value is the same.

## Tests

`validation-suite/test-2026/`: `control_print_gid_group1`, `print_gid_group1` (checked by `scripts/build_safe.sh` in a scratch copy).
