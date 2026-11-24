# slide_axisymmetric

## Implementación

- **Keyword**: `slide_axisymmetric` (INTEGER, length 1, data_class SLIDE),
  registered in `database.cc`. Enum `SLIDE_AXISYMMETRIC` in `tochnog.h` /
  `tochnog-mod.h` (alphabetical order, before `SLIDE_FRICTION`).
- **Ley**: in `slide()` (`slide.cc`), when `slide_axisymmetric == -YES`, the
  friction force `slide_force = friction * Fn` is scaled by `2*pi*r`, with
  `r = NODE[inod][0]` (the radial distance to the symmetry axis, the
  x-coordinate in axisymmetric), consistent with the element assembly that
  uses `2*pi*r` for axisymmetric volumes (`elem.cc`).

## Validación

- `slide_axi`: axisymmetric quad4, left edge (x=r=1) on a `slide_geometry`,
  volumetric load pushing against the surface (Fn != 0), right face pulled
  upwards. With `slide_axisymmetric -yes` the friction (scaled by 2*pi) brakes
  the slide node; the model converges and the slide node velocity is reduced.
- The dynamic discriminator (with/without flag) is documented as
  problematic (the plan 2026-08-13): the Coulomb friction needs Fn != 0 and
  the 2*pi*r factor is validated by inspection of the direct formula.
