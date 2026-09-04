# slide_plasti_friction

## Description

`slide_plasti_friction` (manual Professional 6.1042) specifies the
plastic friction of a slide (`slide_geometry`):

```
slide_plasti_friction index phi c
```

The maximum friction force between the material and the slide surface
equals `c + Fn*tan(phi)`, where `phi` is the friction angle in
radians and `Fn` the normal force.

Related records of the slide family registered in the GNU:

- `slide_geometry index geometry_entity geometry_entity_index`
  (present in the GNU);
- `slide_plasti_tension index sig_t` (6.1043) — maximum tensile force
  of the slide connection;
- `slide_stiffness index stiffness_n stiffness_t` (6.1046) — elastic
  normal/tangential stiffness of the slide connection;
- `slide_plasti_residual_stiffness index rn rt` — residual stiffness
  fraction after plastification (the Professional defaults to 1e-2);
- `control_slide_plasti_apply index switch` (6.371) and
  `control_slide_stiffness_apply index switch` (6.372) — per-control
  gates of the plastic law / the elastic stiffness.

## Uso

```
slide_geometry         10  -bottom
slide_plasti_friction  10  0.523598 0.
slide_stiffness        10  1.e3 1.e3
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `phi` | Friction angle in RADIANS (coefficient tan(phi)). |
| `c` | Cohesion. |
| `stiffness_n` / `stiffness_t` | Normal / tangential stiffness. |

## Verification / status

Records registered and parsed (slide1/slide4 parse cleanly and run);
the CONSUMPTION of the new records by the slide law is PENDING. The
legacy GNU `slide()` (penalty velocity constraint + friction from the
assembled normal force) does not constrain the slide nodes with the
Professional inputs (the block falls through the plane, top reaction
≈ 0) — see the developer manual.
