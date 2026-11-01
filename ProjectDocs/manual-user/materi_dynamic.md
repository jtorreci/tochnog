# materi_dynamic / control_materi_dynamic

Time-integration blending factor of the dynamic (inertia) calculation
(manual Professional 6.800 / 6.141).

## Syntax

```
materi_dynamic factor
control_materi_dynamic index factor
```

`factor` is a number between 0 and 1. Default, when the record is absent:
`factor = 1` (fully implicit).

## Meaning

The record is specifically meant for dynamic calculations. When the
solution is known at time `t` and the new time step solves for `t+dt`, the
stresses at time `t` are used with `(1-factor)` and the stresses at time
`t+dt` with `factor`:

```
sigma_used = (1-factor)*sigma(t) + factor*sigma(t+dt)
```

A factor below 1 makes the scheme less implicit, and thus less numerical
damping occurs (the displacement history keeps its amplitude). Typical use:

```
inertia_apply -yes
materi_dynamic 0.
```

`control_materi_dynamic index factor` is the same record restricted to the
time steps of `control_timestep` with the same index.

## Example

```
( Test dynamic truss )
echo -no
number_of_space_dimensions 1
materi_displacement
materi_velocity
end_initia
inertia_apply -yes
materi_dynamic 0.
node 1  0
node 2  1
element 1  -truss 1 2
group_type 0  -truss
group_truss_area 0  1.0
group_truss_elasti_young 0  1.0
group_truss_density 0  2.0
bounda_dof    10  1 -velx
bounda_force  20  2 -velx
bounda_time   20  1.
control_timestep 10  1.e-1 100.
end_data
```

## Implementation status

Implemented as a momentum blend in the constitutive feedback:

- continuum elements (`materi()`, stress.cc path): the integration-point
  stress used for the internal force is
  `(1-factor)*sigma_old + factor*sigma_new` and the momentum stiffness is
  scaled by `factor`;
- truss elements (`truss()`): the same blend applied to the element axial
  force; the force *state* (`element_truss_force`) always keeps the full
  update.

Known limitation (measured against the Professional binary 25-10-2023): the
GNU dynamics scheme solves velocities with a lumped inertia and assembles
the stiffness implicitly; lowering `factor` reduces the algorithmic damping
but also reduces the stability region towards the explicit limit. The
Professional (displacement-based formulation) stays stable with
`factor = 0` for the same time steps (dynamic1/dynamic8 of the corpus). The
GNU therefore cannot yet reproduce the dissipation-free target values of
the corpus dynamics family: dynamic1/2/5/8 run but fail their targets
(RUNFAIL). This is a solver/scheme blocker (the staggered velocity
formulation of the GNU), not a parsing one.

## Effective factor (F4 T6d, 2026-10-11)

`control_materi_dynamic` with the index of the current time step overrules
`materi_dynamic` (manual 6.141), also for the displacement-primary scheme
used with a factor below 1 (see the developer page
`dynamic_displacement_primary.md`). Before 2026-10-11 that scheme was only
switched on by the global `materi_dynamic`, so an input with only
`control_materi_dynamic 10 0.` ran the old blended path (Pro `dynamic1`
variant: node 2 disx 6.4477 instead of 0.0958890). Measured equal to the
Professional 02-08-2026: only `control_materi_dynamic 10 0.` 0.0958890,
global + control 0. 0.0958890, control 0.5 0.932818, no record 0.995505.
Tests: `validation-suite/test-2026/control_materi_dynamic_only1`,
`control_materi_dynamic_half1`, `truss_dynamic_no_record1`.

## Differences with the Professional version

- None known for the factor resolution: the former C11 difference
  (`control_materi_dynamic` ignored by the displacement-primary gate,
  truss default memory) was resolved on 2026-10-11 (F4 T6d; BACKLOG C11).
  Index: [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
