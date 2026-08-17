# group_materi_plasti_mohr_coul_direct / tension_direct

## Description

`group_materi_plasti_mohr_coul_direct` and
`group_materi_plasti_tension_direct` are **direct stress cut-off** plastic
laws for materials (soils, rock, masonry joints). Unlike the incremental
plasticity models, they do not use plastic strains: they **cap the traction
on a specific plane** with normal vector `n` ("cut off by Tochnog").

- `group_materi_plasti_tension_direct sigy`: caps the NORMAL traction on the
  plane to `sigy` (tensile limit).
- `group_materi_plasti_mohr_coul_direct phi c phi_flow`: caps the SHEAR
  traction on the plane to the Mohr-Coulomb limit
  `max_fric = max(c - sig_n*tan(phi), 0)` where `sig_n` is the normal
  traction (compression increases the limit).

These laws are the material counterpart of the interface laws
(`group_interface_materi_plasti_*_direct`) and are useful to limit stress on
a known plane (e.g. a joint, a bedding plane, a foundation interface) inside
a continuum.

## Uso

Place it in the data part, in the element group definition:

```
group_type 10  -materi
group_materi_elasti_young 10  1000.0
group_materi_elasti_poisson 10  0.0
group_materi_plasti_tension_direct 10  1.0
group_materi_plasti_tension_direct_normal 10  0. 1. 0.
```

or with the Mohr-Coulomb shear limit:

```
group_materi_plasti_mohr_coul_direct 10  0.785398  1.0  0.0
group_materi_plasti_mohr_coul_direct_normal 10  0. 1. 0.
```

## Parámetros

| Record | Parameters | Meaning |
|--------|------------|---------|
| `group_materi_plasti_mohr_coul_direct` | `phi c phi_flow` | Friction angle (rad), cohesion, dilatancy angle (rad, accepted for compatibility, no effect on the cut-off). |
| `group_materi_plasti_mohr_coul_direct_normal` | `nx ny nz` | Plane normal (explicit). |
| `group_materi_plasti_mohr_coul_direct_normal_automatic` | `-yes` | Take the normal from the element normal. |
| `group_materi_plasti_mohr_coul_direct_visco` | `tm` | Visco relaxation time. |
| `group_materi_plasti_mohr_coul_direct_wall` | `phi c phi_flow` | Values used when the element is attached to a wall. |
| `group_materi_plasti_tension_direct` | `sigy` | Tensile limit on the plane normal traction. |
| `group_materi_plasti_tension_direct_normal` | `nx ny nz` | Plane normal (explicit). |
| `group_materi_plasti_tension_direct_normal_automatic` | `-yes` | Take the normal from the element normal. |
| `group_materi_plasti_tension_direct_visco` | `tm` | Visco relaxation time. |
| `group_materi_plasti_tension_direct_wall` | `sigy` | Tensile limit used when the element is attached to a wall. |

The `_normal` and `_normal_automatic` records are optional. Without them the
normal defaults to zero and the laws have no effect; you must provide one of
the two to define the plane.

## Notas

- Requires `materi_stress` in the initialization part.
- The cut-off is applied on the elastic stress BEFORE the incremental
  plastic-yield test, so it combines with the incremental plasticity models.
- `phi_flow` is accepted for interface compatibility but has no effect (no
  plastic flow rule in the direct laws).
- `_visco tm`: viscoplastic relaxation. The stress interpolates between the
  elastic and the fully capped (plastic) response with
  `factor = 1 - exp(-dt/tm)`. `dt << tm` → elastic, `dt >> tm` → fully capped.
- `_wall`: alternative parameters used when the element is attached to a
  wall (its nodes belong to a group listed in `group_materi_plasti_boundary`).

## Estado de implementación

- **Implementado**: `group_materi_plasti_mohr_coul_direct` (+ `_normal`,
  `_normal_automatic`, `_visco`, `_wall`) and `group_materi_plasti_tension_direct`
  (+ `_normal`, `_normal_automatic`, `_visco`, `_wall`). Validated with
  `materi_direct` (uniaxial tension capped to sigy), `materi_direct_mc`
  (pure shear capped to max_fric), `materi_direct_auto` (element-normal
  automatic), `materi_direct_visco` (visco relaxation) and
  `materi_direct_wall` (wall values).
- **Pendiente**: softening via dependency diagrams (`materi_strain_total_shear_kappa`).
