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

- The rotation handling (no rotation of the stress), the strain measure
  (linear engineering strains) and the compressibility block of
  `-updated_area` are the same as `-updated_linear`. One kinematic
  difference: `-updated_area` builds the strain increment on the CURRENT
  geometry, `-updated_linear` on the reference geometry, so on a mesh that
  follows the material the accumulated strains differ slightly (1D bar
  compressed 0.5%: -0.0050117 vs -0.005).
- The Professional 02-08-2026 accepts `-updated_area` only in the
  `incremental_driver` context: in an ordinary `group_materi_memory` it
  stops with "Error detected for data item : group_materi_memory". The GNU
  accepts it everywhere; outside the incremental driver there is no
  Professional reference.
- 2026-09-29: the stress is no longer rotated by the incremental rotation
  and the compressibility increment uses the incremental strain (same fix
  as `-updated_linear`, audit cluster K01).
- The difference with `-updated_linear` (area update of the loaded
  surface) only becomes observable through the `incremental_driver`
  experiment machinery (force-controlled experiments on a deforming
  specimen). That machinery is NOT implemented in the GNU yet, so today
  both options behave identically in the GNU.
- All four corpus tests that use `-updated_area` (the three above plus
  `incremental_driver_syntax`) remain blocked by the `incremental_driver`
  family, which the GNU does not implement (see SEGUIMIENTO-CONVERGENCIA).
