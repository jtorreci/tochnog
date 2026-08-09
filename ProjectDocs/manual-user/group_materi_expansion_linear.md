# group_materi_expansion_linear (and group_materi_expansion_volume)

## Description

Thermal expansion of materials under a temperature field (`condif_temperature`).

- `group_materi_expansion_linear alpha` — linear (directional) thermal
  expansion. It contributes temperature strains on the diagonal of the strain
  tensor:

  ```
  inc_temperature_strain[ii] = -alpha * (T_new - T_old)
  new_temperature_strain[ii] = -alpha * T_new
  ```

  (note the minus sign: the strain is subtracted so that a positive
  temperature rise produces the usual positive expansion).

- `group_materi_expansion_volume beta` — volumetric expansion. It modifies the
  material density used for gravity forces:

  ```
  dens = (1 - beta*T) * dens
  ```

Useful for thermal cracking, restrained thermal deformation, and coupled
thermo-mechanical analyses of embankments, pavements, or concrete structures.

Both require the temperature unknown `condif_temperature`; they are read with
`GET_IF_EXISTS` and are inactive when the record is absent.

## Usage

```
group_materi_expansion_linear <element_group>  alpha
group_materi_expansion_volume <element_group>  beta
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `alpha` | Linear thermal expansion coefficient (1/deg). |
| `beta`  | Volumetric thermal expansion coefficient (1/deg), applied to the density. |

## Example

```
condif_temperature
materi_velocity
materi_strain_total
materi_stress
...
group_materi_elasti_young 0  1.0
group_materi_expansion_linear 0  1.0
group_materi_memory 0 -updated_without_rotation
```

This is the `expans1.dat` regression test (linear thermal expansion).
