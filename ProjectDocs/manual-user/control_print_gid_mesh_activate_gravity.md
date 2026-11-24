# control_print_gid_mesh_activate_gravity

Manual Professional 6.305 and print_gid_mesh_activate_gravity 6.985.

## Description

Accepted (`-yes` or `-no`). The manual only says "See also mesh_activate_gravity_time". The GNU
GiD mesh always holds every element, also the ones that `mesh_activate_gravity_time` has not
activated yet at the first print.

The records apply to the GiD files of `control_print_gid` (`<base>.flavia.msh`, `<base>.flavia.res`). The mesh is written once, at the first GiD print; the results at every print. Overview of the family: [print_gid_options](../manual-developer/print_gid_options.md).

## Syntax

```
control_print_gid_mesh_activate_gravity  index  switch
print_gid_mesh_activate_gravity  switch
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Index of a `control_print_gid` record. |
| `switch` | `-yes` or `-no`. |

## Example

```
control_print_gid                        10  -yes
control_print_gid_mesh_activate_gravity  10  -yes
```

## Differences with the Professional version

- Measured on the Professional 02-08-2026 with `mesh_activate_gravity_time`: by default (and
  with `-no`) its mesh holds only the elements already activated at the first GiD print; with
  `-yes` it holds every element. The GNU always writes every element (the Professional `-yes`)
  and the switch has no effect. Open question to the owner.

## Tests

`validation-suite/test-2026/`: `control_print_gid_mesh_activate_gravity1` (checked by `scripts/build_safe.sh` in a scratch copy).
