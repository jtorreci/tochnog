# element_3d

## Files and functions

- `polynom.cc` — shape functions and integration points for all elements.
- `elem.cc` — generic element loop (uses `ndim` + the shape functions from `polynom.cc`).
- `database.cc` — keyword registration (`-tet4`, `-tet10`, `-hex8`, `-prism6`) and `data_length[ELEMENT]`.
- `check.cc` — validation (`check()` returns 0 for `HEX8`, `HEX27`, `TET4`, `TET10`).
- `initia.cc` — `npointmax` default (=6).

## Implementation details

The 3D element support is built in three ways:

1. **TET4** (`polynom.cc:150`) — explicit shape functions using area
   coordinates `L1,L2,L3,L4`. Minimal (1 point) or maximal (4 points, one per
   node) integration.
2. **TET10** (`polynom.cc:199`) — explicit quadratic tetrahedron shape
   functions, 10 integration points.
3. **HEX8** — uses the general tensor-product formula (`polynom.cc:312-411`).
   `npol=2` for `-hex8` (line 313) → `nnol = nnol_xi*nnol_eta*nnol_zeta = 8`.
   Integration via `integration_gauss`/`integration_lobatto` per direction,
   `npoint_per_dir = npol = 2` for `-maximal` → 2×2×2 = 8 points. The shape
   functions are the tensor product of 1D polynomials, and the Jacobian is
   computed generically with `matrix_inverse` (3×3) in the block at
   `polynom.cc:413+`.
4. **PRISM6** (`polynom.cc:146`) — explicit prism (wedge) shape functions:
   area coordinates for the triangular base × linear interpolation in z.
   `N_i = (1-z)*L_i` for the base nodes (1-3) and `N_i = z*L_i` for the top
   nodes (4-6). 6 integration points (2 z-levels × 3 area points). Derivatives
   are set for the 3 directions (xi, eta from area coords; z from the linear
   interpolation).

`elem.cc` is generic: it does not special-case elements. It reads `ndim`,
calls `polynom` for `h`/`p`, and assembles using the shape function gradients.
No changes to `elem.cc` were needed to add 3D elements.

## External dependencies

- LAPACK (`matrix_inverse` for the 3×3 Jacobian in 3D) — already linked.

## Hardcoded parameters / pending refactorings

- `npointmax = 6` default in `initia.cc:44`. HEX8 needs 8 and PRISM6 needs 6,
  so `number_of_integration_points` must be set in the initia part (≥8 for
  hex8, ≥6 for prism6).
- `HEX27`/`HEX64` are referenced by the tensor formula (`npol=3/4`) but not
  fully registered/validated in all paths.
- PRISM6 integration is fixed at 6 points (no `-minimal`/`-maximal` selection
  yet).
- Consider making the 3D element registration data-driven (a table of nodes,
  integration points, and shape function generators) to reduce duplication.
