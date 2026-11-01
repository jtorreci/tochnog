# slide_plasti_friction / slide family (slide_geometry law)

## Description

The slide family models a material sliding over a fixed geometry
(`slide_geometry`, manual Professional 6.1042): the slide nodes are
elastically attached to the plane (springs on the TOTAL displacements)
with a Mohr-Coulomb plastic friction cap. This reproduces the
Professional slide1/slide4 tests (block dragged over a plane with
`Fn = kn*penetration` and friction `mu*Fn`).

Records of the family:

```
slide_geometry            index geometry_entity geometry_entity_index
slide_stiffness           index stiffness_n stiffness_t       (6.1046)
slide_plasti_friction     index phi c                         (6.1042)
slide_plasti_tension      index sig_t                         (6.1043)
slide_plasti_residual_stiffness index rn rt
control_slide_plasti_apply     index switch                   (6.371)
control_slide_stiffness_apply  index switch                   (6.372)
```

The node membership is the geometry test on the START coordinates (the
nodes belong to the plane where they were created even after sliding far
along it) plus the explicit `node_slide index slide_number` records
(manual 6.893).

## The slide law (per slide node)

- **Normal spring** (two-sided): `S = -kn*(u.n)` on the node along the
  plane normal n (compression positive). `slide_plasti_tension sig_t`
  only caps the maximum TENSILE force the connection can take
  (`S >= -sig_t`); absent = no tension limit. The node is elastically
  attached to the plane in both directions.
- **Tangential spring**: elastic predictor `-kt*u_t` on the tangential
  displacement, capped by the Mohr-Coulomb limit
  `|F_t| <= c + |Fn|*tan(phi)` (`slide_plasti_friction phi c`, phi in
  RADIANS; the magnitude of the normal force, also when the connection
  is in tension, as in the Professional). Once the cap is exceeded the node slips plastically and the
  tangential force stays at the cap (perfect plasticity; the direction
  opposes the tangential displacement).
- **Axisymmetric problems**: every force/stiffness acts on the whole
  ring: multiplied by `2*pi*r` (r = node radius), triggered by
  `slide_axisymmetric -yes` or by a problem group with
  `group_axisymmetric -yes`. The plain record spelling
  `axisymmetric -yes` / `-no` (manual 6.17, no index) is accepted and
  stored on group 0.
- **Output records** per slide node (Professional .dbs):
  - `node_slide_direction` — (n, t) local frame, t = friction direction
    on the node;
  - `node_slide_f` — plastic yield function value (0 = at the cap);
  - `node_slide_force` — the force the material applies ON the slide
    geometry in the local frame: slot 0 = normal component
    (`kn*un`, compression negative), slot 1 = - tangential force
    magnitude.
- `control_slide_plasti_apply -no` disables the plastic caps (pure
  elastic springs); `control_slide_stiffness_apply -no` disables the
  springs.

## Uso

```
slide_geometry         10  -bottom
slide_plasti_friction  10  0.523598 0.      (phi c: 30 degrees, no cohesion)
slide_stiffness        10  1.e3 1.e3        (kn kt per node)
slide_plasti_residual_stiffness 10  1.e-2 1.e-2   (matrix-only, optional)
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `phi` | Friction angle in RADIANS (coefficient tan(phi)). |
| `c` | Cohesion. |
| `stiffness_n` / `stiffness_t` | Normal / tangential spring stiffness, per node. |
| `sig_t` | Maximum tensile (pull-off) force of the connection. Absent record = no tension limit (since 2026-10-11, F4 T6e; before, an absent record wrongly behaved as `sig_t = 0`, a compression-only gap). |
| `rn` / `rt` | Residual stiffness FRACTION of the elastic stiffness added to the matrix only while the node slips (default 1e-2, like the Professional). |

## Verification / status

- `slide1`, `slide3` and `slide4` of the corpus pass (rc=0); `slide4`
  equals the Professional since 2026-10-11 (F4 T6e). slide1's final state
  matches the Professional .dbs: node_slide_force (-2.1148e-3,
  -1.2210e-3) vs (-2.1129e-3, -1.2199e-3); the top reactions +0.01 /
  -5.7707e-3 (targets).
- `slide2` (axisymmetric) still fails — see the developer manual for
  the diagnosis.
- `slide_friction_cohesion1` (test-2026, F4 T6e): slide4 with cohesion
  `slide_plasti_friction 10 0.523598 0.002`, end time 10: node 1
  tangential `node_slide_force` -2.065482e-3 = Professional 02-08-2026
  (needs the absent-tension default, the `|Fn|` cap and the post-solve
  material state pass).

The legacy penalty behavior (velocity constraint) is kept for slide
geometries WITHOUT `slide_stiffness` (classic inputs such as examp14,
test-2014 `node_slide`): the node gets the reaction
`- slide_penalty * (v.n) * dtime * n` plus the `slide_friction` force on
its VELOCITY equations (fixed 2026-10-11: the terms used to land in other
unknowns, e.g. the total strain `eptxx`, of models with
`materi_strain_total` and `derivatives`).
