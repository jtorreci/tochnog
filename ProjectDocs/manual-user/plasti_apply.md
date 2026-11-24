# plasti_apply

## Description

`plasti_apply switch` (manual Professional 6.899, global, no index): with
`-no`, any plasticity data in the input file is ignored for all timesteps,
so the calculation runs "just linear" (elastic) without commenting out the
plasticity records. The per-timestep form is `control_plasti_apply index
switch` (6.257, see `control_plasti_apply.md`).

## Usage

```
plasti_apply switch
control_plasti_apply index switch
```

| Parameter | Meaning |
|-----------|---------|
| `switch` | `-yes` (default) applies plasticity, `-no` ignores the plasticity data |
| `index` | (control form) the `control_timestep` index the switch applies to |

## What is ignored (measured on the Professional 02-08-2026)

With `-no` the result equals the model with the plasticity records removed
(or, for laws that are also the stiffness, the frozen state):

- yield-surface models: `group_materi_plasti_mohr_coul`, `_vonmises`,
  `_druck_prag`, `_tension`, `_camclay` (the preconsolidation history does
  not evolve), `_hardsoil` (plastic part; the `group_materi_elasti_hardsoil`
  stiffness stays), and every other law evaluated by the plastic yield test;
- direct cut-offs: `group_materi_plasti_mohr_coul_direct`,
  `_tension_direct`, `_compression_direct` (including `_wall`/`_visco`
  variants);
- hypoplasticity (`group_materi_plasti_hypo_*`): the stress and the
  hypoplastic history (void ratio) stay frozen;
- incremental models (`_tskh`, `_aitskh`, `_camclay_incremental`): elastic
  response instead;
- slide springs (`slide_plasti_*`) when no `control_slide_plasti_apply`
  record is given (manual 6.371).

Plastic strains reached in earlier timesteps are kept (a later `-no`
timestep continues elastically from them).

Not affected: the truss/spring/interface/beam plasticity records
(`group_truss_plasti`, `group_spring_plasti`,
`group_interface_materi_plasti_*`) the per-group switch
`group_plasti_apply` is described in `group_plasti_apply.md`; see the
differences section.

## Precedence

`control_plasti_apply` with the index of the current timestep overrules
`plasti_apply` (measured on the Professional 02-08-2026); default `-yes`.

## Example

```
group_materi_plasti_mohr_coul 0 0.5 1. 0.
plasti_apply -no
```

Self-checking inputs in `validation-suite/test-2026/`: `plasti_apply1.dat`
(Mohr-Coulomb -> elastic sigyy -21.978), `plasti_apply2.dat`
(Professional `hypo1.dat` -> sigyy -200, hyhis0 0.6), `plasti_apply3.dat`
(Professional `validation_13.dat` Cam clay -> cchis1 207),
`control_plasti_apply1.dat` / `control_plasti_apply2.dat` (per-timestep
switch and precedence). All pass on the Professional 02-08-2026.

## Differences with the Professional version

- Truss, spring, interface and beam plasticity are not gated (not measured
  on the Professional); `linear_calculation_apply`, which sets
  `plasti_apply -no` in the Professional, is not implemented.
- `group_materi_plasti_user` / umat laws are not gated (user code).
- The plastic results themselves (switch `-yes`) keep the pre-existing
  GNU/Professional differences of the return algorithms (e.g. Mohr-Coulomb
  sigyy -3.487 vs -3.394 in the test model); the `-no` results are
  identical.
