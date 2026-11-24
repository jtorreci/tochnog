# control_materi_plasti_hardsoil_gammap_initial

## Description

`control_materi_plasti_hardsoil_gammap_initial` (manual Professional
6.146) adds an **extra initial contribution to gamma_p** exactly such
that the hardsoil yield function is zero-valued at the initial stress
state. This is convenient to start a calculation with deviatoric
stresses which would be OUTSIDE the yield surface without this extra
contribution (the HS surface passes through `q = 0` at `gamma_p = 0`:
any initial deviatoric stress is outside).

With the switch set to `-yes`, the first timestep of the
`control_timestep` record with the same index computes

```
gamma_p_extra = q/(E50*(1 - q/qa)) - 2*q/Eur
```

at the INITIAL stress state (i.e. the value of the yield function at
`gamma_p = 0`), stores it in the record
`element_intpnt_materi_plasti_hardsoil_gammap_initial` and adds it to
`gamma_p` inside the hardsoil yield function from then on. The initial
state is then ON the surface (`f = 0`) and no initial plastic return
occurs.

Default switch is `-no`.

## Uso

```
control_materi_plasti_hardsoil_gammap_initial 0  -yes
```

## Parámetros

| # | Parameter | Meaning |
|---|-----------|---------|
| 1 | switch | `-yes` creates the extra initial gamma_p, `-no` does nothing (default) |

## Validation

- `mhardsoil_gp0` (initial deviatoric stress `sigxx = -2`, no loading,
  control `-yes`): `gamma_p_extra = 2/(1000*(1-2/38.49)) - 4/3000 =
  0.0007763` EXACT in `element_intpnt_materi_plasti_hardsoil_gammap_initial`,
  the yield function is zero at the start, the initial stress stays
  `sigxx = -2.0` EXACT and `kappa` stays 0.
- `mhardsoil_gp0_off` (A/B without the control): the same state is
  outside the surface, the return relaxes the deviatoric stress
  (`sigxx -> -1.225`) and `kappa` grows to `0.000293`.
