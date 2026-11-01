# inertia_apply

## Description

Global switch that includes or excludes the inertia terms of the principal
dofs: material mass (for `materi_velocity`), heat capacity (for
`condif_temperature`) and water/flow storage capacity (for
`groundflow_pressure`).

- `-yes`: the inertia terms are included — dynamic analyses, transient heat
  conduction and transient ground water flow.
- `-no` (default): the inertia terms are skipped; the corresponding equations
  become quasi-static / steady.

The record is global (`no_index`) and accepts a single switch value that
applies to every principal dof (the Professional also accepts one switch per
principal dof; the GNU reads the first value only). Default is `-no`,
matching the Professional (manual 6.785: "Default, if inertia_apply is not
specified, then each of switch_0, switch_1 etc. is set to -no").

Use `control_inertia_apply` to switch per timestep.

## Usage

```
inertia_apply <switch>
```

`switch` is `-yes` or `-no`. Default is `-no`.

## Parameters

| Parameter | Meaning                                                        |
|-----------|----------------------------------------------------------------|
| `switch`  | `-yes` include the inertia terms, `-no` skip them. Default `-no`. |

## Example

```
inertia_apply -yes
```

A dynamic analysis: the material mass matrix is included and the momentum
equation carries the `M·a` term.

```
groundflow_consolidation_apply -yes
inertia_apply -yes
```

A transient consolidation analysis: both the consolidation coupling and the
groundflow capacity storage (`C·ṗ`, the inertia term of the pressure dof) are
active (the corpus examples `large2`/`large3` set both records).

## Notes

- Lumped masses: every element uses the row-sum lumped mass of its
  integration rule, except the `tet10`, which uses the positive HRZ masses
  (corner V/36, mid-edge 4V/27) since F4 T4 (2026-10-09), and the
  `tria6`, which uses HRZ masses since F4 T6f C66 (2026-10-16): plane
  corner A/19 and mid-side 16A/57, axisymmetric the r-weighted diagonal
  scaled to 2 pi A r_c, for every `group_integration_points` value.
  Example: `validation-suite/test-2026/tria6_mass_hrz1.dat` (one tria6 of
  area 1/2, unit point forces, one step dt 0.01: v = dt/m = 0.38 at a
  corner, 0.07125 at a mid-side node).

- The default changed from `-yes` (GNU legacy) to `-no` (Professional) in the
  u-p consolidation sprint (2026-09-05, commit `0560dcf`). The legacy GNU
  default made every model with a material density run a pseudo-dynamic
  analysis even when the load was quasi-static; the corpus safety/groundflow
  tests (ground14/15/16) run static in the Professional and settle within
  their 1 s window only with the static default.
- Transient analyses that rely on the mass/heat-capacity/storage terms must
  set the switch explicitly. The corpus dynamic (`dynamic*`), transient heat
  (`bicg1`, `phase1`, ...) and consolidation (`large2`, `large3`) examples
  all do so.

## Differences with the Professional version

- **Lumped mass / capacity of the `tet10`** (deliberate divergence,
  project owner 2026-10-09; DIVERGENCES TET10-HRZ). GNU: positive HRZ
  lumping, corner V/36 and mid-edge 4V/27 (the diagonal of the
  consistent mass scaled to the element total). Professional 02-08-2026:
  row-sum lumping of its 5-point rule, corner -V/20 and mid-edge V/5;
  the negative corner masses make dynamics unstable (sudden gravity on a
  48-tet10 column: the Professional solver fails, a GNU row-sum build
  diverges to 5e8). Static results are unaffected.
- **Lumped mass / capacity of the `tria6`** (deliberate divergence,
  project owner 2026-10-16; DIVERGENCES TRIA6-HRZ). GNU: positive HRZ
  lumping of every h-weighted dof (velocity, temperature, pressure,
  density, the h-weighted state variables): plane corner A/19 and
  mid-side 16A/57; axisymmetric the r-weighted diagonal int N_i^2 2 pi r
  dA scaled to 2 pi A r_c; the same masses for every
  `group_integration_points` value. Professional 02-08-2026: row sum,
  corner 0 in plane analysis (a free one-element dynamic step fails with
  a zero pivot), negative corners in axisymmetric analysis (-pi/60 at
  r = 1 for r = 1..2). Static results are unaffected. Tests
  `tria6_mass_hrz1`, `tria6_mass_hrz_axisym1`. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
