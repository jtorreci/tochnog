# strain_settlement_parameters

## Description

`strain_settlement_parameters` adds an extra **vertical settlement creep
strain** for soil dumping: after a material is dumped (activated), it
continues to settle over time. The vertical creep strain of a soil particle
is assumed to be (saturating power law; the manual OCR is ambiguous, this
law was decided 2026-08-18):

```
eps_zz(t) = Ar * (t/t_ref)^n / (t_plus + (t/t_ref)^n)
```

with:
- `Ar` = reference creep strain rate (reference_creep_strain_rate),
- `t_plus` = time_plus,
- `t_ref` = reference time,
- `n` = power constant,
- `t` = the time elapsed after the material has become active (dumping).

The horizontal creep strains are `eps_xx = eps_yy = lateral_factor * eps_zz`.

`strain_settlement_element_group` selects the element groups (or `-all`).

This record should be combined with `mesh_activate_gravity_time` +
`mesh_activate_gravity_time_strain_settlement -yes` (which determines when
each element becomes active).

## Uso

Place it in the data part:

```
strain_settlement_parameters 0  0.0 1.0 0.01 1.0 1.0 0.0
```

Parameters: `time_global_start time_plus reference_creep_strain_rate
reference_time power_n lateral_factor`.

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `time_global_start` | Global time when creep starts. |
| `time_plus` | Time constant of the creep law. |
| `reference_creep_strain_rate` | Ar (creep rate). |
| `reference_time` | Reference time t_ref. |
| `power_n` | Power constant n. |
| `lateral_factor` | Horizontal/vertical strain ratio (0 for uniaxial). |

## Dependencia de parámetros

`strain_settlement_diagram dof_0 value_0 dof_1 value_1 ...` +
`strain_settlement_diagram_dof dof` + `strain_settlement_diagram_number n`
make a parameter of `strain_settlement_parameters` depend on a dof value
(1 = time_plus, 2 = reference_creep_strain_rate, 3 = reference_time,
4 = power_n, 5 = lateral_factor) via the table.

## Related

- `mesh_activate_gravity_time` (+ `_time_strain_settlement -yes`) — the
  element activation that sets the "time after dumping".

## Estado de implementación

- **Implementado**: `strain_settlement_parameters` (creep applied in
  `set_deften_etc`, `strain_settlement_creep` in materi.cc),
  `strain_settlement_element_group`, and the `strain_settlement_diagram*`
  parameter dependence. Validated with `strain_settle` (disy grows with the
  creep) and `strain_settle_diag` (Ar doubled via the diagram).
