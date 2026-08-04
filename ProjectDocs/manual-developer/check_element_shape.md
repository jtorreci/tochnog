# check_element_shape

## Archivos y funciones

- `polynom.cc` — the check lives at the end of the element integration loop,
  after the `volume[ipoint] = weight[ipoint]*detj/...` assignment block
  (`polynom.cc:443-458`). The code is inside the same loop function that fills
  `volume[ipoint]` per integration point.
- `initia.cc:100` — global `double check_element_shape_factor=0.;`
  (default: check off).
- `top.cc:131` — reads the keyword once at startup with
  `db( CHECK_ELEMENT_SHAPE, 0, idum, &check_element_shape_factor, ldum, VERSION_NORMAL, GET_IF_EXISTS )`.
  Note the value arrives through the `dval` argument because the type is
  `DOUBLE_PRECISION`, not `ival`.
- `database.cc:265-268` — keyword registration:
  `strcpy(name[CHECK_ELEMENT_SHAPE],"check_element_shape")`, `type = DOUBLE_PRECISION`,
  `data_length = 1`, `no_index[CHECK_ELEMENT_SHAPE] = 1`.
- Enum `CHECK_ELEMENT_SHAPE` in `tochnog.h:161` / `tochnog-mod.h:154`
  (must stay in sync).

## Detalles de implementación

- The check only runs when `check_element_shape_factor>0.`; a non-positive
  value disables it.
- It uses `volume[ipoint]` as the distortion proxy. `volume[ipoint]` is
  proportional to `detj` (per element type the constants differ: `1/2` for
  triangles, `1/6` for tetrahedra, `2`/`4`/`8` times `detj` for 1D/2D/3D), so
  comparing volume ratios equals comparing determinant ratios within one
  element.
- Distortion per element: `avg( |volume - vol_avg| / vol_avg )` over the
  integration points, guarded by `vol_avg>0.` (a zero average skips the check).
- On a hit it prints to `std::cout`:
  `Warning: element with distorted shape, distortion=... > factor=...`.
  It does not abort the run.
- Note: with the axisymmetric branch (`axisymmetric==-YES`) the volume is
  additionally scaled by `2*PIRAD*radius`, so the distortion metric then
  reflects radius-weighted volumes too.

## Dependencias externas

None. Uses only `scalar_dabs` and locals already in scope in the integration
loop.

## Parameters hardcodeados / refactorizaciones pendientes

- The warning message does not identify WHICH element or integration point is
  distorted — adding element/ipoint context would make the diagnostic usable.
- Output goes to `std::cout`; routing it through `pri()` / `tn.log` would be
  consistent with the rest of Tochnog reporting.
- The per-element-type volume constants (`1/2`, `1/6`, `2`/`4`/`8`) are
  duplicated logic; the distortion metric assumes `volume` is always positive
  (volume is `|detj|`, so it is, but the code relies on that invariant).
