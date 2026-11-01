# group_materi_damage_mazars

## Description

Isotropic damage law of Mazars for concrete-like materials under uniaxial
tension/compression. The scalar damage variable `d` (0 = virgin, 1 = fully
damaged) is a history variable (`materi_damage`), so damage is permanent and
monotonic.

The damage is driven by the equivalent strain `epseq`, computed from the
principal elastic strains and the stress triaxiality `gamma`:

```
epseq = gamma * sqrt(sum of positive principal elastic strains^2)
gamma  = -sqrt(sum(sigma_-^2)) / sum(sigma_-)      (compressive triaxiality)
```

Once `epseq > eps0` the damage grows with the Mazars exponentials for tension
(`at, bt`) and compression (`ac, bc`):

```
dt = 1 - (eps0/epseq)*(1-at) - at*exp(-bt*(epseq-eps0))
dc = 1 - (eps0/epseq)*(1-ac) - ac*exp(-bc*(epseq-eps0))
alpha = epseq^2 / epssiz^2                          (brittleness index)
d  = dt*alpha^beta + dc*(1-alpha)^beta
```

`d` never decreases (`new_damage = max(new_damage, old_damage)`), and the
stress is `(1-d)` times the effective (undamaged) stress, which is kept per
integration point (manual Professional 2.2.5; before 2026-09-29 the damage
compounded every step, see `materi_damage_effective_stress.md`). `gamma` is
evaluated with the effective stress components xx, xy, xz (not the
principal stresses), as in the Professional. With two or three positive
principal strains the damage differs from the Professional, which uses a
non-standard equivalent strain (see "Differences with the Professional
version").

Useful for cracking/failure of concrete, mortar, rock and other quasi-brittle
materials.

Requires `materi_stress`, `materi_damage`, `materi_strain_elasti`,
`materi_strain_total`, and a velocity formulation.

## Usage

```
group_materi_damage_mazars <element_group>  eps0 at bt ac bc beta
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `eps0` | Damage threshold strain (below it the material is elastic). |
| `at` | Post-peak exponential coefficient in tension. |
| `bt` | Post-peak exponential rate in tension. |
| `ac` | Post-peak exponential coefficient in compression. |
| `bc` | Post-peak exponential rate in compression. |
| `beta` | Brittleness exponent weighting tension vs compression damage. |

## Example

```
materi_stress
materi_damage
materi_strain_elasti
materi_strain_total
...
group_materi_elasti_young 0  23400.e6
group_materi_damage_mazars 0  2.6e-4  1.  15000.  1.2  649.  1.05
group_materi_memory 0 -updated_without_rotation
```

This is the `damage1.dat` regression test (uniaxial extension, damage growth).

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

## Integration-point state (C68, 2026-10-12)

`materi_damage` is kept per integration point between steps (Professional manual
4.10), as the hardening parameter `kap`. The node values printed by
`node_dof` are a recovery (averaging) of the integration-point values
with the same weights as the strains; the constitutive update always
uses the integration-point value. Before 2026-10-12 GNU fed the averaged
node value back into the update, which smeared the history across
elements (up to 19 % off the Professional in the reactions of a small
damage test).
