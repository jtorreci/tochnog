# groundflow_phreatic_level_static — developer

Mechanics-only static groundflow (manual Professional 6.574): with
`groundflow_phreatic_level_static -yes` and no `groundflow_pressure`
unknown, the element groups with `-groundflow` in `group_type` get the
static pore pressure of the phreatic level and the wet/dry density.
Audit cluster K22 (`AUDIT-ALL-VARIABLES-2026.md`), fixed 2026-09-29.

## Files and functions

- `database.cc` — record `GROUNDFLOW_PHREATICLEVEL_STATIC` (INTEGER,
  no index, `-yes`/`-no`). Registered earlier (see
  `force_element_edge_water.md`).
- `calcul.cc` — `groundflow_static_pressure_active()`: 1 when the single
  switch or any `groundflow_phreatic_level_multiple_static` is `-yes`. It
  also enables the post_calcul pressure monitors without a pressure
  unknown (see `groundflow_multilevel_mechanics.md`, dev/excv).
- `groundfl.cc`
  - `groundflow_phreatic_coord()` — without a pressure unknown and with a
    static level, returns `total = static = rho_w*g*(z_level - z)`,
    clamped to `groundflow_pressure_atmospheric`, plus
    `groundflow_addtopressure`. `level_source` is 1 when a level covers the
    point.
  - `groundflow_static_pressure_group( element_group )` (new) — 1 when
    there is no pressure unknown, `groundflow_static_pressure_active()` is
    true, and the `GROUP_TYPE` of the group holds `-GROUNDFLOW`.
- `materi.cc`
  - `materi()`, block "forces on nodes": new `else if
    ( groundflow_static_pressure_group( gr ) )` branch next to the
    `groundflow_pressure` branch. It calls `groundflow_phreatic_coord( -1,
    coord_ip, ... )` at the integration point. When `level_source==1` it
    adds `groundflow_pressure_factor * total_pressure` to the diagonal of
    `total_new_sig` (the stress used for the internal force). The stress
    dofs keep the effective stress (`new_sig`).
  - `get_materi_density()` — the wet/dry rule of the pressure-dof branch
    (average of the nodal `groundflow_phreatic_coord` totals, wet if < 0)
    now also runs when `groundflow_static_pressure_group()` is true. The
    old `NODE_PHREATICLEVEL` branch stays for the other cases.
- `tochnog.h` / `tochnog-mod.h` — prototypes of
  `groundflow_static_pressure_active()` (was local to `calcul.cc`) and
  `groundflow_static_pressure_group()`, same place in both headers.

## Implementation details

- **Why both parts are needed.** Before the fix the `NODE_PHREATICLEVEL`
  branch chose the density. That record is only written by
  `groundflow_phreaticlevel_bounda`, which needs the pressure unknown, so
  without it every element was wet. And the pore pressure was only added
  inside `if ( groundflow_pressure )`. The column of `ground18` carried
  the full wet weight as effective stress: bottom sigyy -1995 instead of
  -1247.5. The density fix alone gives -1745 (the total stress), the
  pressure fix alone -1497.5.
- **Tangent.** The static pressure does not depend on any unknown, so it
  is a load: no coupling term in the element matrix (the
  `groundflow_pressure` coupling block of `materi()` stays gated on the
  pressure unknown).
- **Group gate.** Measured on the Professional 02-08-2026: `ground18` with
  `group_type -materi` only gives -1995 (wet density everywhere, no pore
  pressure). So both the density choice and the pore pressure need
  `-groundflow` in the group.
- **Multiple levels.** `groundflow_static_pressure_group()` also covers
  `groundflow_phreatic_level_multiple_static`, but `check.cc` still
  demands the pressure unknown for `groundflow_phreatic_level_multiple`, so
  that path is not reachable without it today.
- **Output.** Unchanged: `calcul.cc` already evaluated `-total_pressure`
  etc. from the same `groundflow_phreatic_coord()`.

## Differences with the Professional (measured, 02-08-2026)

The Professional applies the static pore pressure to the skeleton only
when `post_calcul -groundflow_pressure -total_pressure` is present.
`-static_pressure` does not trigger it. Measured variants:

| input | Pro | GNU after |
|---|---|---|
| `ground18` (has the post_calcul) | -1247.5 | -1247.5 |
| `ground18` without post_calcul | -1745 | -1247.5 |
| `ground18` with `-static_pressure` instead | -1745 | -1247.5 |
| `ground18` without density_groundflow, without post_calcul | 0 | 497.5 |
| `ground18` without density_groundflow, with post_calcul | 497.5 | 497.5 |
| `excavate1` (no post_calcul), node 13 vely | -0.00210987 | -0.00253891 |
| `excavate1` + post_calcul, node 13 vely | -0.00253891 | -0.00253891 |

The GNU follows the manual (the pore pressure is part of the physics, not
of the output request). Consequence: the Professional target of
`excavate1` (-0.00210987, produced without the pore pressure) now fails in
the GNU. See BACKLOG C10.

- **`post_calcul -groundflow_pressure -total_pressure` with a pressure dof**
  (GNU extension / input dialect; AUDIT 1.1; BACKLOG C10). GNU: Not
  required. Professional 02-08-2026: Required when a pressure dof exists
  (and couples the static pore pressure, C10). Status: Input dialect. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).

## Verification

- `ground18` (Pro test/other) vs Pro: every node_dof and node_dof_calcul
  equal (sigyy 1e-9, total pressure exact); identical to its
  pressure-dof twin `ground4`.
- `validation-suite/test-2026/phreatic_static1.dat` (level inside the
  column) and `phreatic_static_submerged1.dat` (level 10 above the top,
  `force_edge_water` on the top): targets = analytic discrete values,
  equal to the Pro; rc=1 on the pre-fix binary (-1995 and -2095).
- Regression: every `.dbs` of validation-suite, sfnet and test-2026
  identical to the pre-fix binary (no corpus test uses the switch
  without a pressure unknown).

## Pending

- Decide whether to follow the Professional coupling to the post_calcul
  record (BACKLOG C10).
