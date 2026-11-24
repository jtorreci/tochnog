# control_change_dataitem_apply

## Description

`control_change_dataitem_apply` allows you to enable or disable the
`change_dataitem*` records for the timestep record with the same index.

If the switch is set to `-no`, any `change_dataitem*` data in the input
file is ignored for that control index. This is useful, for example, to
stop applying a prescribed displacement after a certain time step.

## Uso

Place it in the data part, with the same `icontrol` index as the
`change_dataitem` records you want to control:

```
control_change_dataitem_apply 20  -no
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `20`      | Index of the control record. Must match the `change_dataitem` and `control_timestep` records. |
| `-no`     | Ignore all `change_dataitem*` records for this control index. |

Without this record (or with `-yes`), the `change_dataitem*` records are
applied (default).

## Example

Stop prescribing the displacement of a geometry line after control step 20:

```
change_dataitem                 20  -geometry_line 1 1 -use
change_dataitem_time            20  0. 10. 100. 0.
control_change_dataitem_apply   20  -no
```

With `-no`, the `change_dataitem` record above is ignored.
