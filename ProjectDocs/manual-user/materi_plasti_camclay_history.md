# materi_plasti_camclay_history

## Description

`materi_plasti_camclay_history` (manual Professional 4.16) adds the two
history variables of the Modified Cam Clay model to the `node_dof`
records:

- `cchis0` — the void ratio `e0`;
- `cchis1` — the preconsolidation pressure `p0`.

No number follows the keyword: the two variables are fixed. The initial
values must be set by the user (typically with `control_reset_dof`
before the loading steps, together with the initial stresses), exactly
as the Professional manual prescribes: "The history variables e, p0
need to be initialized by materi_plasti_camclay_history record (and
given initial values in node_dof records)".

With this record declared, `group_materi_plasti_camclay` takes the
Professional layout `index M kappa lambda` (manual 6.689, three
values) and tracks `p0` explicitly as a history variable. The legacy
GNU layout (a fourth parameter `N`, `p0` derived from the current
pressure through the normal compression line) is still accepted when
four values are given AND `cchis1` is not initialized.

## Uso

In the initialization part:

```
materi_velocity
materi_displacement
materi_strain_total
materi_strain_plasti
materi_stress
materi_plasti_camclay_history
end_initia
```

Data part (initial state of an isotropically consolidated sample):

```
group_materi_plasti_camclay    0  0.882  0.031  0.088   (M kappa lambda)
group_materi_elasti_camclay_g  0  3000.
group_materi_elasti_camclay_pressure_min  0  1.

control_reset_dof               1 -sigxx -sigyy -sigzz
control_reset_value_constant    1 -207.
control_reset_dof               2 -cchis0
control_reset_value_constant    2 0.627
control_reset_dof               3 -cchis1
control_reset_value_constant    3 207.0
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `cchis0` | Void ratio `e` history (node_dof name). |
| `cchis1` | Preconsolidation pressure `p0` history (node_dof name). |

The values are updated by the Cam Clay law during the calculation
(volume change for `e`, plastic volumetric hardening for `p0`).

## Verification

`validation_13` (undrained triaxial compression of Weald clay, strain
controlled): rc=0 — final `p0` = 262.3 matches the Professional binary
(25-10-2023, .dbs) and the target 264.5±10; `e` stays 0.627 (constant
volume).
