# post_calcul_materi_stress_force_thickness_switch

## Description

`post_calcul_materi_stress_force_thickness_switch` (manual
Professional 6.917): in 3D Tochnog normally assumes that the SHORTEST
element direction in the evaluated side is the structure thickness
direction. If that is not the case (e.g. very short elements in the
tunnel length direction), set the switch to `-yes` for the groups that
must use the LONGEST element direction as thickness direction. One
switch per element group of
[`post_calcul_materi_stress_force_element_group`](post_calcul_materi_stress_force_element_group.md).

## Input syntax

```
post_calcul_materi_stress_force_thickness_switch switch_0 switch_1 ...
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `switch_i` | `-yes` (longest direction is thickness) or `-no` (default, shortest direction). One switch per element group. |

## Notes

- The number of switches must equal the number of element groups;
  otherwise an error is raised.
- Consumed by the numerical integration (lots 2/3), not yet
  implemented: the record is validated and stored.
