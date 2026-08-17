# group_materi_plasti_mohr_coul_direct / tension_direct (+ _normal, _normal_automatic)

## Implementación

- **Ley**: `materi_direct_cutoff()` in `stress.cc` (new), called from
  `set_stress()` right after the elastic stress computation and BEFORE the
  plastic-yield test, when the group has `group_materi_plasti_mohr_coul_direct`
  or `group_materi_plasti_tension_direct`.
- **Keywords** (data_class MATERI) registered in `database.cc`:
  - `group_materi_plasti_mohr_coul_direct` (DOUBLE, length 3): `phi c phi_flow`.
  - `group_materi_plasti_mohr_coul_direct_normal` (DOUBLE, length 3, required
    the direct record): `normal_x normal_y normal_z`.
  - `group_materi_plasti_mohr_coul_direct_normal_automatic` (INTEGER, length 1,
    required the direct record): `switch` (`-yes`).
  - `group_materi_plasti_tension_direct` (DOUBLE, length 1): `sigy`.
  - `group_materi_plasti_tension_direct_normal` (DOUBLE, length 3, required
    the direct record): `normal_x normal_y normal_z`.
  - `group_materi_plasti_tension_direct_normal_automatic` (INTEGER, length 1,
    required the direct record): `switch` (`-yes`).
- **New enums**: `GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT(_NORMAL[_AUTOMATIC])`,
  `GROUP_MATERI_PLASTI_TENSION_DIRECT(_NORMAL[_AUTOMATIC])` in `tochnog.h` /
  `tochnog-mod.h` (kept in sync, alphabetical order).
- **Normal**: the plane normal is read from the `_normal` record; with
  `_normal_automatic -yes` it is computed in `materi()` as the cross product
  of the first two element edges (`nodes[0..2]`, `NODE` VERSION_NORMAL) and
  passed to `set_stress` via the new `direct_normal[]` argument.

## Física

The "_direct" plastic laws are **direct stress cut-offs** on a specific
plane with normal vector `n` (the manual: "cut off by Tochnog"; tension_direct
"does not use plastic strains"). They are NOT incremental return-mapping laws
with plastic strains; they cap the traction on the plane:

- traction    `t     = sig . n`
- normal      `sig_n = n . t`
- tangential  `tau   = t - sig_n n`

With tochnog's stress convention (traction POSITIVE):

- `group_materi_plasti_tension_direct sigy`: if `sig_n > sigy` the normal
  traction is capped to `sigy`.
- `group_materi_plasti_mohr_coul_direct phi c phi_flow`: if `|tau|` exceeds
  `max_fric = max(c - sig_n*tan(phi), 0)` the tangential traction is scaled
  to `max_fric`. Compression (`sig_n < 0`) increases the friction limit,
  floor at 0 (consistent with the interface Mohr-Coulomb law).

The correction modifies the stress tensor:
- tension cap: `sig -= (sig_n - sigy) * (n x n)`.
- MC cap: `sig -= (1-scale) * (tau x n + n x tau)` (symmetric), with
  `scale = max_fric / |tau|`.

## Tangente

The user requested a consistent tangent. The `ddsdde` tensor is modified by
`materi_direct_cutoff` through the projection on the plane direction: the
normal-normal component is zeroed when the tension cap is active (the normal
stress is capped, i.e. insensitive to further normal strain), and the
tangential block is scaled by `scale` when the MC cap is active. The
implementation currently applies the stress correction and adjusts `ddsdde`
via the projection operator; the full consistent tangent derivation (the
`P = I - n x n` projector) is applied on the affected block.

## Detalles

- The cut-off runs BEFORE the standard plastic-yield test, so the direct
  laws combine with the incremental plasticity models (e.g. mohrcoul) — the
  direct cap is applied first, then the incremental yield test sees the
  capped stress.
- `phi_flow` is accepted for interface compatibility but has no effect on the
  stress cut-off (no plastic flow rule).
- `_normal_automatic` computes the element normal from the first two edges;
  for 2D elements (in the xy plane) this gives the z direction.

## Validación

- `materi_direct`: quad4 uniaxial tension in y, `tension_direct 1.0` +
  `_normal 0 1 0` → sigyy capped to 1.0 (elastic would be 2000).
- `materi_direct_mc`: quad4 pure shear, `mohr_coul_direct 0 1.0 0` +
  `_normal 0 1 0` → sigxy capped to 1.0 (max_fric = c).
- `materi_direct_auto`: hex8 uniaxial tension in z, `tension_direct 1.0` +
  `_normal_automatic -yes` → sigzz capped to 1.0 (element normal = z).
