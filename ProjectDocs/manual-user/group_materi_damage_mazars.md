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
stresses are scaled by `(1-d)` in the material update (see
`materi.cc` `materi_damage` block).

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
