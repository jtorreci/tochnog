# group_materi_stress_null / _direction / _direction_automatic

## Description

Sets chosen stress components of the material to zero in a local frame
(manual Professional 6.753-6.755). Typical use: a thin layer of normal
continuum elements behaving as a soil-structure interface, keeping only
the normal stress and the shear stresses on the interface plane (manual
9.1).

- `group_materi_stress_null index dir ...`: the local components to set
  to zero: `-x`, `-y`, `-z` (normal) and `-xy`, `-xz`, `-yz` (shear).
- `group_materi_stress_null_direction index` with nine values: the local
  axes x, y, z (three values each). Each axis is normalized; give
  perpendicular axes.
- `group_materi_stress_null_direction_automatic index -yes`: local x is the
  element normal (`element_normal`, inward normal of the mesh boundary),
  y and z are perpendicular to it.
- Without both direction records the global frame is used. The explicit
  direction wins over the automatic one.

The components are set to zero in the final stress, after plasticity, and
the material stiffness gets the same projection.

## Syntax

```
group_materi_stress_null index -x|-y|-z|-xy|-xz|-yz ...
group_materi_stress_null_direction index x0 x1 x2 y0 y1 y2 z0 z1 z2
group_materi_stress_null_direction_automatic index -yes|-no
```

## Example

```
group_materi_stress_null 10 -y -z -yz
group_materi_stress_null_direction_automatic 10 -yes
```

Tests: `validation-suite/test-2026/group_materi_stress_null1`-`5`, `_error1`,
`_direction_error1`.

## Differences with the Professional version

- The GNU projects the material stiffness too, so a step converges in one
  iteration; the Professional keeps the full stiffness and needs more
  equilibrium iterations (with its default iterations the interface layer
  test stops at node disx 1.17808e-4, converged 1.24290e-4 in both codes;
  DIVERGENCES STRESS-NULL-TANGENT).
- A zero axis in `_direction` is an error; the Professional runs.
- Modes without any stiffness left (for example `-y` on a free single
  element) are indeterminate; the GNU keeps them at zero increment, the
  Professional returns another solution.
- The Professional manual example also uses `global_normal -yes`; the GNU
  generates the element normal by itself and does not know that record
  (BACKLOG C14).
- The perpendicular axes y, z of the automatic frame are a GNU choice; the
  result depends on them only when y and z are not nulled together.
