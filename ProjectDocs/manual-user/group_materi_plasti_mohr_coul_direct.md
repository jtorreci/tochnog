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

## Sprint 10 additions: compression_direct + pressure/coord limits

- `group_materi_plasti_compression_direct index sigy`: principal
  stresses lower than sigy are cut off (direct cut-off, no plastic
  strains; manual Professional 6.694). Spectral implementation:
  sigma = V diag(d) V^T, eigenvalues below sigy pulled up.
- `group_materi_plasti_compression_direct_visco index tm`: relaxation
  time; the cut relaxes with factor 1-exp(-dt/tm) (without the record
  the cut is total).
- `group_materi_plasti_pressure_limit index pressure_limit`: neglect the
  DIRECT plasticity laws when the pressure (positive in compression,
  -tr/3) exceeds the limit — free-surface problems (manual 6.688).
- `group_materi_plasti_coord_limit index coord_limit`: neglect them when
  the vertical coordinate exceeds coord_limit (manual 6.689).
- `group_materi_plasti_heat_generation index factor`: Professional name
  (with underscores) of the legacy `group_materi_plasti_heatgeneration`
  (fraction of plastic energy converted to heat; requires
  condif_temperature).

Tests: mdirect_comp (sigyy -10 capped at -5.0 EXACT), mdirect_gate
(pressure_limit disables the cap -> elastic -10; A/B fails without it).

## Sprint 10 lote 2: Drucker-Prager alias + analytic validation

- `group_materi_plasti_druck_prag index phi c phi_flow` — the
  Professional name (with underscore) of the legacy GNU
  `group_materi_plasti_druckprag` (alias in db_number; same physics,
  layout [phi c phi_flow], legacy tests druckpr1/examp23).
- `group_materi_plasti_bounda` / `_factor` — Professional names of
  `group_materi_plasti_boundary` / `_factor`.
- `group_materi_factor index factor` — multiplication factor for the
  material stresses AND stiffness (unit conversion).
- Registered partials: `group_materi_damping_method`,
  `group_materi_density_groundflow` (pre-existing),
  `group_materi_plasti_visco_exponential_limit/_name/_values`.

### Analytic validation (test mdp_shear)

The GNU implements f = sqrt(J2) + 3·alpha·sigma_m − K (the standard DP
form; NOTE sqrt(J2), not sqrt(3·J2)) with the manual's matching
constants. In PURE SHEAR with phi = 0 (alpha = 0, no pressure term, no
dilatancy) the capped stress is exactly

    sigma_xy = K = 6·c·cos(phi) / ( sqrt(3)·(3 − sin(phi)) ) = 2c/sqrt(3)

With c = 30: sigma_xy = 34.641. Measured at dt = 0.05: 34.636 —
0.014% error. (With phi ≠ 0 the associative flow generates normal
stresses in a displacement-imposed shear rig, so sigma_m ≠ 0; use
phi = 0 for the closed-form check.)
