# group_materi_memory -updated_area

## Description

`-updated_area` is a material memory type of `group_materi_memory`
(manual Professional 6.784, "MATERIAL MEMORY" of the `incremental_driver`
section). It is a SMALL deformation theory, like `-updated_linear`, with
one difference: the area change of the surface where forces are applied IS
taken into account (the forces keep the same stress but act on the CURRENT
area of the deformed specimen).

The three small/large deformation options of the incremental driver are:

| Memory | Deformation theory | Area of the loaded surface |
|--------|--------------------|----------------------------|
| `-updated_linear` | small | NOT updated (initial area) |
| `-updated_area` | small | updated (current area) |
| `-updated` | large | updated (current area) |

The Professional manual documents the semantics only inside the
`incremental_driver` example (section 6.784); the general
`group_materi_memory` section (6.685) lists `-updated`, `-updated_jaumann`,
`-updated_linear`, `-total` and `-total_linear`, but not `-updated_area`.

## Usage

```
group_type            0  -materi
group_materi_memory   0  -updated_area
group_materi_elasti_young    0   100000.
group_materi_elasti_poisson  0   0.23
```

The corpus laboratory tests (oedometric_drained_mohr_coulomb,
triaxial_compression_drained_mohr_coulomb,
triaxial_compression_undrained_mohr_coulomb) use it as:

```
group_type                           0  -materi
group_axisymmetric                   0  -yes
group_materi_memory                  0  -updated_area
```

## Notes

- The kinematic part (strains, rotations) of `-updated_area` is IDENTICAL
  to `-updated_linear`: identity rotation handling and linear engineering
  strains (verified: the GNU runs `-updated_area` and `-updated_linear`
  through the same code path and produces identical results).
- The difference with `-updated_linear` (area update of the loaded
  surface) only becomes observable through the `incremental_driver`
  experiment machinery (force-controlled experiments on a deforming
  specimen). That machinery is NOT implemented in the GNU yet, so today
  both options behave identically in the GNU.
- All four corpus tests that use `-updated_area` (the three above plus
  `incremental_driver_syntax`) remain blocked by the `incremental_driver`
  family, which the GNU does not implement (see SEGUIMIENTO-CONVERGENCIA).
