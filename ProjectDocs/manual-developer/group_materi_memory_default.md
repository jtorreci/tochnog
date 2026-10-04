# group_materi_memory — default (no record)

## Implementation

- **general.cc, `materi_memory_default()`**: single source of the default
  memory. Returns `-UPDATED` when `control_materi_updated_apply` exists
  and its record at the current `ICONTROL` is `-YES`; otherwise
  `-UPDATED_LINEAR`. Declared in `tochnog.h` and `tochnog-mod.h`.
- Every reader of `GROUP_MATERI_MEMORY` assigns the default right before
  its `db(GROUP_MATERI_MEMORY, ..., GET_IF_EXISTS)` (an explicit record
  overwrites it):
  - `materi.cc` `materi()` and `set_deften_etc()`;
  - `stress.cc` `set_stress()` (the `-no` gate of
    `control_materi_updated_apply` stays there, after the read);
  - `elem.cc` `elem()`;
  - `tendon.cc` `tendons()`; its assert now accepts `-UPDATED_LINEAR`
    (same no-rotation tendon direction as `-UPDATED_WITHOUT_ROTATION`).
- **elem.cc, reference gradient**: the `NODE_START_REFINED`
  reference-geometry gradient (`ref_active`) is loaded only when
  `memory==-UPDATED_LINEAR` AND the group has an EXPLICIT record
  (`db_active_index(GROUP_MATERI_MEMORY, element_group, VERSION_NORMAL)`).
  The default keeps the current-geometry gradient, which matches the
  Professional on `-follow_material` meshes (stretch probe 0.695653 in
  both; the reference gradient gives 1).
- Previously the default `-UPDATED` was hard-coded as the initial value
  of `memory` at those five places.

## Verification (2026-10-04)

| case | GNU before | GNU after | Professional 02-08-2026 |
|---|---|---|---|
| shear gamma 1, no record (sxx / sxy / syy) | 0.2311 / 0.4207 / -0.2286 | 0 / 0.5 / 0 | 0 / 0.5 / 0 |
| same + `control_materi_updated_apply 0 -yes` | 0.2311 / 0.4207 / -0.2286 | 0.2311 / 0.4207 / -0.2286 | 0 / 0.5 / 0 |
| `-follow_material` stretch x2, no record (sxx) | 0.695653 | 0.695653 | 0.695653 |

Tests: `validation-suite/test-2026/materi_memory_default_shear1`
(rc=1 before, rc=0 after), `materi_memory_default_stretch1` and
`materi_memory_default_updated_apply1` (regression guards: the old
default was already `-updated`; the latter fails with the new binary if
the `-yes` line is removed). Registered in `scripts/build_safe.sh`.

Corpus A/B (60 s, timeouts rerun at 300 s; load-induced timeouts rerun
on both binaries): validation-suite 233/1, test-2026 44 (+3 new),
sfnet 187 PASS unchanged; Pro `other`: `mohr_coul_direct6` and
`mohr_coul_direct8` RUNFAIL -> PASS, `hypo12` PASS -> RUNFAIL (marginal,
sigyy off 1.02e-3 vs tolerance 1e-3; Professional distribution file,
target not edited); Pro `validation` unchanged.

## Hardcoded parameters / pending

- `-updated_jaumann` (manual 6.685) is not implemented.
- With `materi_displacement` the Professional default behaves as
  `-total_linear`. A variant mapping the default to `-TOTAL_LINEAR` in
  that case (A2) was rejected: it breaks `mohr_coul_direct7` and sfnet
  `temp5`. The remaining difference comes from the coordinate update in
  `elem.cc` (around line 569) and is tied to the mesh default (BACKLOG C8).
- `control_materi_updated_apply -yes` follows the manual; the
  Professional binary ignores it (deliberate divergence).
