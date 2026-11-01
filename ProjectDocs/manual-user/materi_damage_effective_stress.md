# materi_damage_effective_stress (damage stress behaviour)

## Description

With `materi_damage` the stress follows the Professional manual 2.2.5:

```
sigma_damaged = (1 - d) * sigma_undamaged
```

`materi_stress` (the `sig..` dofs) is the DAMAGED stress; the undamaged
(effective) stress is kept internally per integration point and evolves
with the normal material law. This is not a new keyword: it is how every
model with `materi_damage` and `group_materi_damage_mazars` behaves since
2026-09-29.

Before that date the GNU version (and the 2014 sfnet source) multiplied the
already damaged stress by (1 - d) again in every step, so the damage
compounded: a bar pulled to eps = 6e-4 in 60 steps gave 1436.83 instead of
85598.32, and at constant strain the stress decayed to zero. Results of
older GNU runs with damage should be recomputed (`dam` itself was correct).

What it means for a model:

- At constant strain the damaged stress stays constant.
- Plasticity combined with damage: the yield test acts on the damaged
  stress (as in the Professional).
- `materi_damage_apply -no` / `control_materi_damage_apply index -no`: in
  those timesteps the damage data is ignored: `dam` keeps its value and the
  stress increments of the step are added without the (1 - d) factor (as
  in the Professional).

## Usage

No input. Initialize `materi_stress` and `materi_damage` and give a damage
law:

```
materi_stress
materi_damage
materi_strain_elasti
materi_strain_total
end_initia
...
group_materi_damage_mazars 0  2.6e-4 1. 15000. 1.2 649. 1.05
```

## Example

```
echo -no
number_of_space_dimensions 1
materi_velocity
materi_displacement
materi_strain_elasti
materi_strain_total
materi_stress
materi_damage
end_initia
node 1  0
node 2  1
element 1  -bar2 1 2
group_type 0  -materi
group_materi_elasti_young 0  23400.e6
group_materi_damage_mazars 0  2.6e-4 1. 15000. 1.2 649. 1.05
group_materi_memory 0  -total_linear
bounda_dof 0  1 -velx
bounda_time    0  0. 0. 100. 0.
bounda_dof 1  2 -velx
bounda_time    1  0. 1. 100. 1.
control_timestep 0  1.e-5 6.e-4
end_data
```

Result: `dam` 0.9939032534, `sigxx` 85598.32 = (1 - d) * E * eps (same in
the Professional 02-08-2026). Self-checking inputs:
`validation-suite/test-2026/damage_stress1.dat` (load),
`damage_stress_hold1.dat` (load and hold), `damage_stress_3d1.dat` (3D,
compression and shear, load and hold).

## Differences with the Professional version

- **Mazars equivalent strain** (deliberate divergence, class A, project
  decision 2026-10-12; BACKLOG C7). GNU: manual 2.2.5, eps_eq =
  sqrt(sum <eps_i>+^2) of the positive principal strains (hex8
  diag(2e-4,-3e-4,1e-4): dam 0.28610207515). Professional 02-08-2026:
  accumulates them in ascending order as s = 2 s + eps^2 (weights 4, 2, 1;
  the 2014 accumulation defect the GNU fixed on 2026-08-10): dam 0.36409.
  It differs whenever two or three principal strains are positive, also
  through a von Mises + damage path (bar with `vonmises 3e6`: GNU dam
  0.925372, eppxx 1.87863e-4; Professional 0.96946, 1.368e-4). Guard
  tests `test-2026/damage_mazars_eq_two_positive1`,
  `damage_mazars_vonmises1`. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- The triaxiality factor gamma is built from the xx, xy, xz stress
  components in both codes (not objective: a compressive xx with xy or xz
  shear gives an orientation-dependent `dam`). Same in the Professional,
  so not a divergence; noted in BACKLOG C7, no fix planned.
- `group_materi_memory -updated` without `materi_displacement`: the total
  strain differs slightly (5.99823e-4 vs 6.0e-4 in the bar example), so
  `dam` and the stress differ by ~0.2 %. Pre-existing, not related to
  damage.
- Re-enabling damage after timesteps with `control_materi_damage_apply -no`
  and a non-zero damage: not verified against the Professional.
- A restart from a printed database starts from the damaged stress (the
  effective stress is internal and not printed).
