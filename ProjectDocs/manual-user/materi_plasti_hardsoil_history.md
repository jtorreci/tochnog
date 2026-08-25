# materi_plasti_hardsoil_history

## Description

`materi_plasti_hardsoil_history` (manual Professional 4.22) adds the
history variable **abs(p)** (the maximum pressure history) of the
Hardening-Soil model to the `node_dof` records. It is the switch of the
elastic law between first loading (`E50`) and unloading/reloading
(`Eur`): while the current |p| is the new maximum the material is on
first loading, when |p| drops below the maximum it is unloading.

This is the SAME concept as `materi_stress_pressure_history` (manual
4.50): both initialize the shared `sph` dof (same basename, same
running-max update in the solver). Initializing either record is
sufficient for the hardsoil loading/unloading switch.

## Uso

In the initialization part:

```
materi_plasti_hardsoil_history
```

## Parámetros

None (switch-only initialization).

## Validation

- `mhardsoil_elast` / `mhardsoil_elast2`: during first loading the sph
  dof tracks the growing |p| (`0.0732` / `0.0538` EXACT, the running
  maximum over the stress dofs).
- `mhardsoil_unload`: the sph dof stays at the peak `0.3333` while
  unloading (the maximum is kept) and the unloading branch (`Eur`)
  engages exactly at the first unloading step.
