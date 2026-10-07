# group_materi_memory — default (no record)

## Description

When a material group has NO `group_materi_memory` record, the GNU now
uses `-updated_linear`, the default of the Professional version: the
stress is the previous stress plus the stress of the linear incremental
strain, without rotation. Until 2026-10-04 the GNU used `-updated`
(rotated stress, rate-of-deformation strains), which gives a nonlinear,
Jaumann-like response in large shear where the Professional stays linear.

An explicit `group_materi_memory` record always wins over the default.

The default can be switched to `-updated` for selected timesteps with
`control_materi_updated_apply index -yes` (manual Professional 6.152;
see [control_materi_updated_apply](control_materi_updated_apply.md)).

## Usage

No input is needed: omit `group_materi_memory`. Example (single quad4,
simple shear gamma = 1, young 1, poisson 0):

```
echo -no
number_of_space_dimensions 2
materi_velocity
materi_velocity_integrated
materi_stress
end_initia
node 1  0 0
node 2  1 0
node 3  0 1
node 4  1 1
element 1  -quad4 1 2 3 4
group_type 0  -materi
group_materi_elasti_young 0  1.0
group_materi_elasti_poisson 0  0.0
bounda_dof 0  -ra 1 2 -ra -velx
bounda_dof 1  -ra 3 4 -ra -velx
bounda_dof 2  -ra 1 2 3 4 -ra -vely
bounda_time 1  0. 1. 1. 1.
control_timestep 0  0.01 1.
control_timestep_iterations 0  20
end_data
```

Result: `sigxy` 0.5, `sigxx` = `sigyy` = 0 (as the Professional). With
`group_materi_memory 0 -updated` the result is sigxx 0.2311, sigxy
0.4207, sigyy -0.2286 in both versions.

## Notes

- On a `-follow_material` mesh the default keeps the CURRENT-geometry
  gradient, as the Professional (uniaxial stretch to twice the length:
  sigxx 0.695653 in both). The reference-geometry gradient that the GNU
  uses for an EXPLICIT `-updated_linear` (see
  [group_materi_memory_updated_linear](group_materi_memory_updated_linear.md))
  is not applied to the default.
- Self-checking inputs: `validation-suite/test-2026/materi_memory_default_shear1.dat`,
  `materi_memory_default_stretch1.dat`,
  `materi_memory_default_updated_apply1.dat`.

## Differences with the Professional version

- With `materi_displacement` kinematics and no memory record the
  Professional behaves as `-total_linear`; the GNU default
  `-updated_linear` still updates the coordinates, so large deformations
  differ (2-quad4 column sheared by 1: node 3 disx GNU 0.472, Pro 0.5).
  Give `group_materi_memory -total_linear` explicitly if that response is
  needed (BACKLOG).
- `control_materi_updated_apply -yes`: the GNU follows the manual (an
  unspecified memory becomes `-updated`); the Professional 02-08-2026
  binary ignores `-yes` and stays linear.
- `-updated_jaumann` (manual 6.685) is not implemented.
- Corpus test `hypo12` (Professional distribution, no memory record)
  misses its sigyy target marginally since this change (off 1.02e-3 with
  tolerance 1e-3).

- **`-updated` large-rotation robustness (`update2`)** (known difference,
  not decided; BACKLOG C39). GNU: "solution seems to diverge" (veliy -86434)
  with p-refinement. Professional 02-08-2026: veliy 12.243 (target 12.24).
  Status: Open, LOW; not investigated (related to AUDIT K19). Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **`group_materi_memory -updated_without_rotation`** (GNU extension / input
  dialect; AUDIT 1.1). GNU: GNU-only memory (legacy). Professional
  02-08-2026: Not in manual 6.686; the Pro's own copies of these tests use
  `-updated_linear` or `-total_linear`. Status: GNU extension; kept. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **`group_materi_memory -total_piola`** (GNU extension / input dialect;
  manual pages only). GNU: Accepted (2014 source; total formulation).
  Professional 02-08-2026: Not listed in manual 6.686; binary behaviour NOT
  measured. Status: GNU extension; Pro behaviour undetermined. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
