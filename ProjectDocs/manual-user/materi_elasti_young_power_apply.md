# materi_elasti_young_power_apply

## Description

`materi_elasti_young_power_apply` (manual Professional 6.801) is a
switch that ignores any nonlinearity in the young modulus dependent on
a power law: with `-no`, the **constant young of the
`group_materi_elasti_young_power` record (`E0`) is applied at all
times**. This is convenient for testing an input file linearly, without
commenting out the record.

The record is registered under the GNU control convention
`control_materi_elasti_young_power_apply`; the manual name
`materi_elasti_young_power_apply` is accepted as an alias.

## Uso

```
group_materi_elasti_young_power 0
                        1000.   ( E0 )
                        500.    ( E1 )
                        800.    ( E2 )
                        3000.   ( E3 )
                        1.      ( p1 )
                        1.      ( alpha )

materi_elasti_young_power_apply 0  -no
```

With the default (`-yes` or record absent) the power law is applied;
with `-no` the modulus is the constant `E0` of the record.

## Parámetros

| Record | Parameters | Meaning |
|--------|------------|---------|
| `control_materi_elasti_young_power_apply` (alias `materi_elasti_young_power_apply`) | `-yes` / `-no` | `-no`: the power-law nonlinearity is ignored and the constant young `E0` of the young_power record is applied at all times. |

## Física

Same physical model as `group_materi_elasti_young_power`, with the
pressure-dependent term switched off: `E = E0` (no `p/p1` term, no
`E2`/`E3` caps needed). The stiffness is the linear elastic one built
with `E0`.

## Validation

- `myoung6_apply.dat`: same oedometer as `myoung6.dat` with the gate
  set to `-no`: `E = E0 = 1000` constant -> `sigma_zz = -1.3462` and
  `sigma_xx = -0.5769` EXACT (identical to the plain linear base,
  A/B against `myoung6.dat` whose fixed point gives
  `sigma_zz = -1.6657`).

## Notas

- Applies only to `group_materi_elasti_young_power` (the
  `group_materi_elasti_young_polynomial` strain-dependent law is not
  affected).
