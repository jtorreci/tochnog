# group_materi_stokes

## Description

Switch per element group that selects Stokes flow for a viscous material
(manual Professional 6.752). With `-yes` the momentum convection term
`rho * v * grad(v)` is left out of the momentum balance of the elements of
that group (creeping flow). With `-no` the term is assembled, so the flow
is a Navier-Stokes flow, as long as convection is active
(`convection_apply -yes` or `control_convection_apply`).

The switch only removes the convection of the velocity (momentum). The
other convected fields of the material, for example `materi_stress` and
the history variables, are still convected when `convection_apply -yes`.

Why it is useful: Stokes flow is a good model for slow, very viscous
flows. Leaving out the non-linear convection term also makes the
analysis easier to converge.

## Syntax

```
group_materi_stokes index switch
```

| parameter | meaning |
|---|---|
| `index` | element group (see `element_group`) |
| `switch` | `-yes` (Stokes flow, no momentum convection) or `-no` (momentum convection when convection is active) |

Default (no record for a group): `-yes`, as in manual 6.752. Before
2026-10-12 the GNU default was `-no`. A eulerian flow that needs the
momentum convection therefore needs an explicit
`group_materi_stokes index -no`.

The record needs `materi_velocity` in the initialisation part.

## Example

Channel flow of the Professional `validation_1` (eulerian mesh,
convection on, no inertia) as a Navier-Stokes flow:

```
materi_velocity
materi_stress
end_initia
mesh -fixed_in_space -fixed_in_space
inertia_apply -no
convection_apply -yes
...
group_type 0  -materi
group_materi_density 0  1.
group_materi_viscosity 0  1.e-2
group_materi_elasti_compressibility 0  1.e-2
group_materi_memory 0  -updated_linear
group_materi_stokes 0 -no
```

Without the last line the same model is a Stokes flow (default `-yes`).
Own tests: `validation-suite/test-2026/stokes_default1.dat`,
`stokes_group1.dat` (two groups, different switches) and
`stokes_group2.dat`.

## Differences with the Professional version

- **index of the record** (deliberate divergence, class A; DIVERGENCES
  `STOKES-IDX`, owner decision 2026-10-12, Professional-suspect). GNU:
  reads the switch per element group, as manual 6.752 says. Professional
  02-08-2026: uses the value of index 0 for every element group and
  ignores the other indices. Measured on the `validation_1` flow: a group
  1 with `-no` that holds the post_point gives velx 0.124105 in the GNU
  and 0.1465713 (the all `-yes` result) in the Professional. With one
  group, or with the same switch in every group, both programs agree.
- **Stokes flow with convected stress** (known difference, class B;
  BACKLOG C72). With `-yes` and `convection_apply -yes`, the Professional
  still convects the integration-point stress, and that changes the
  velocity. The GNU convects only the nodal stress output. Without the
  record, `validation_1` gives velx 0.1465713 in the Professional and
  0.144159 in the GNU. The GNU value is the same as with `convection_apply
  -no` (Professional 0.1441588). The gap grows for a stiffer
  compressibility: with 1e-5, Professional 0.1442 against 0.1386.
- See also [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
