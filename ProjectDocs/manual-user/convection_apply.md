# convection_apply

## Description

Switch that allows (or forbids) the convection of a material with respect
to the mesh (manual Professional 6.395). If the switch is set to `-yes`
the convective transport terms of the transported fields (velocity,
temperature, stress, history variables, ...) are assembled; if set to
`-no` they are suppressed.

Scope (measured on the Professional 02-08-2026): the switch only governs
the convection of the MATERIAL with respect to the mesh, i.e. the
transport by the material velocity. The heat transport by a prescribed
`group_condif_flow` and by the groundflow velocity (`groundflow_velocity`
dof) is always applied, whatever the value of `convection_apply`.

The record applies for ALL timesteps. The per-timestep form is
`control_convection_apply` (manual Professional 6.113), which applies for
the timestep records with the same index and overrides the global switch.

The record is the Professional input name of the GNU record
`options_convection`, which the element integration already reads
(`general()` in general.cc). Registering the alias makes Professional
input files work unchanged.

## Syntax

```
convection_apply -yes
convection_apply -no
```

The value is a switch: `-yes` or `-no`. No index is used.

Default, if neither `convection_apply` nor `control_convection_apply` is
set: `-no` (manual Professional 6.395; changed 2026-10-06, K03). An
eulerian (`mesh -fixed_in_space`) analysis whose material should be
convected therefore needs an explicit `convection_apply -yes`.

## Example

A 1D heat transfer + flow problem (convection of temperature by the
velocity field, corpus `condif2.dat` / `condif3.dat`):

```
mesh -fixed_in_space
convection_apply -yes
...
group_type 0  -materi -condif
group_condif_density 0  1.0
group_condif_capacity 0  1.0
group_condif_conductivity 0  0.1
```

A stationary Navier-Stokes tube flow that disables the convective term
`rho * v * grad(v)` (corpus `tube1.dat`):

```
mesh -fixed_in_space -fixed_in_space
convection_apply -no
inertia_apply -no
```

## convection_stabilization

Switch that controls the amount of artificial diffusion added to damp the
wiggles of the convective terms (manual Professional 6.396). Values:
`-yes` (minimal artificial diffusion, the default), `-maximal` (maximal
artificial diffusion, plus a bounding clamp for condif temperatures),
`-no` (no stabilization). The record is the Professional input name of
the GNU record `options_stabilization`, consumed by `general()` in
general.cc (peclet-based exponential upwind).

## control_convection_apply

Indexed per-timestep form of `convection_apply` (manual Professional
6.113). It is applied for the timestep records with the same index and
overrides the global switch:

```
control_convection_apply 10 -no
```

## Related records

- `node_convection_apply` (manual 6.880): per-node deactivation of the
  convection contributions. Registered in the GNU as a partial (parses,
  no behaviour) - not covered by this batch.
- `condif_convection_*` / `condif_radiation_*`: different feature
  (convection/radiation boundary conditions of heat conduction), see
  their manuals.

## Tests

`validation-suite/test-2026/convection_apply_default1.dat` (the condif2
model without the record: middle temperature 0.5, pure diffusion, as the
Professional) and `convection_apply_yes1.dat` (same model with `-yes`:
1). The validation-suite inputs `cmat_gate`, `cd_method` and
`materi_direct_wall` carry an explicit `convection_apply -yes` since
their targets were calibrated on the old default.

## Differences with older GNU versions

Before 2026-10-06 the GNU default was `-yes`, and `-no` also switched off
the `group_condif_flow` and groundflow-velocity heat transport. Old
inputs that relied on the default (sfnet `condif2`, `condif3`, `examp1`,
`incnav5`, `strconv1`) now give the Professional results instead of
their legacy targets (BACKLOG B10).

Legacy corpus notes:

Corpus tests that pass with these aliases: `condif2`, `condif3`
(node 2 temperature 0.992028 vs target 1.0 +/- 2e-2; the Professional
gives 1.000000 on the same model) and `tube1` (post-point vely 1066.59
vs target 1000 +/- 100; the Professional gives 1000.49).

`validation_1` parses and runs with the aliases but remains RUNFAIL: the
eulerian velocity self-advection transient of the GNU diverges from the
Professional (per-step post_point_dof trajectory and final .dbs differ;
final state velx +0.0723 vs -0.0391). The eulerian NS transient needs a
dedicated solver/physics sprint - see SEGUIMIENTO-CONVERGENCIA.md.

## Differences with the Professional version

- **explicit `convection_apply -no` next to `materi_displacement`** (GNU
  extension / input dialect; AUDIT 1.1). GNU: Accepted (convection ignored
  with `materi_displacement`). Professional 02-08-2026: Refuses the explicit
  record (`-no` is its default). Status: Input dialect. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
