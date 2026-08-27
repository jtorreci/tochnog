# post_calcul_materi_stress_force_element_group

## Description

`post_calcul_materi_stress_force_element_group` specifies the element
groups for which `post_calcul -materi_stress -force` (manual
Professional 6.913) determines the normal force, shear force and
moment(s). It is the MANDATORY configuration record of the -force
family: without it the calculation aborts with a clear error. The
option is meant for isoparametric elements (-quad4, -quad9, -hex8,
-hex27) with a single element over the structure thickness (sheet
piles, tunnel shells, ...).

## Input syntax

```
post_calcul_materi_stress_force_element_group element_group_0 element_group_1 ...
post_calcul -materi_stress -force
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `element_group_i` | Element group whose elements produce forces and moments. At least one group is required. |

## Example

```
post_calcul_materi_stress_force_element_group 0
post_calcul -materi_stress -force
```

## Notes

- The record belongs to the `post_calcul` configuration: it is stored
  in data class POST (no index, shared by all `-force` records).
- The numerical integration is NOT implemented yet (lot 1 =
  infrastructure: registration, dispatch, validation, print structure).
  Until lot 2/3 land, the calculated values are 0 and a
  "not yet implemented" notice is printed once per run.
