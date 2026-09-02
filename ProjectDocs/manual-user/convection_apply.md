# convection_apply

## Description

Switch that allows (or forbids) the convection of a material with respect
to the mesh (manual Professional 6.395). If the switch is set to `-yes`
the convective transport terms of the transported fields (velocity,
temperature, stress, history variables, ...) are assembled; if set to
`-no` they are suppressed.

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

The value is a switch: `-yes` (default in the GNU element integration)
or `-no`. No index is used.

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

Corpus tests that pass with these aliases: `condif2`, `condif3`
(node 2 temperature 0.992028 vs target 1.0 +/- 2e-2; the Professional
gives 1.000000 on the same model) and `tube1` (post-point vely 1066.59
vs target 1000 +/- 100; the Professional gives 1000.49).

`validation_1` parses and runs with the aliases but remains RUNFAIL: the
eulerian velocity self-advection transient of the GNU diverges from the
Professional (per-step post_point_dof trajectory and final .dbs differ;
final state velx +0.0723 vs -0.0391). The eulerian NS transient needs a
dedicated solver/physics sprint - see SEGUIMIENTO-CONVERGENCIA.md.
