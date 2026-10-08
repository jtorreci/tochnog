# group_materi_memory — default (no record)

## Implementation

- **general.cc, `materi_memory_default()`**: single source of the default
  memory. Returns `-UPDATED` when `control_materi_updated_apply` exists
  and its record at the current `ICONTROL` is `-YES` (checked first);
  otherwise `-TOTAL_LINEAR` when `materi_displacement` is initialised
  (manual Professional 6.797, F3 2026-10-07, C18) and `-UPDATED_LINEAR`
  else. Declared in `tochnog.h` and `tochnog-mod.h`.
- Every reader of `GROUP_MATERI_MEMORY` assigns the default right before
  its `db(GROUP_MATERI_MEMORY, ..., GET_IF_EXISTS)` (an explicit record
  overwrites it):
  - `materi.cc` `materi()` and `set_deften_etc()`;
  - `stress.cc` `set_stress()` (the `-no` gate of
    `control_materi_updated_apply` stays there, after the read);
  - `elem.cc` `elem()`;
  - `tendon.cc` `tendons()`; its assert accepts `-UPDATED_LINEAR` and,
    since F3, `-TOTAL_LINEAR` (same no-rotation tendon direction as
    `-UPDATED_WITHOUT_ROTATION`).
- `elem.cc` skips the coordinate update of `materi_displacement` for
  `-TOTAL_LINEAR` (`if ( materi_displacement && memory!=-TOTAL_LINEAR )`),
  so the `-TOTAL_LINEAR` default removes the C18 column difference.
- The former elem.cc reference gradient (`ref_active`, explicit
  `-updated_linear` only) was removed in F3; see
  `group_materi_memory_updated_linear.md` section 7.
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

- F3 (2026-10-07, design A): the `materi_displacement` -> `-TOTAL_LINEAR`
  default (variant A2, rejected on 2026-10-04 alone) was adopted together
  with the `-fixed_in_space` mesh default. Known consequences: Pro
  `mohr_coul_direct7` (`-total_linear` with a velocity `mpc_node_number`
  tie) now fails (separate GNU defect, BACKLOG C-item); sfnet `temp5` (and
  seven other sfnet targets) were calibrated on the old GNU defaults
  (BACKLOG B-item).
- Test: `validation-suite/test-2026/materi_memory_default_displacement1`
  (column node 3 disx 0.5 = Pro; old 0.472222).
- `control_materi_updated_apply -yes` follows the manual; the
  Professional binary ignores it (deliberate divergence).

## Differences with the Professional version

- **`control_materi_updated_apply -yes`** (deliberate divergence, project
  decision 2026-10-04; BACKLOG C21; SEGUIMIENTO Descarte; test-2026
  `materi_memory_default_updated_apply1`). GNU: Follows manual 6.152: an
  unspecified memory becomes `-updated` (simple shear 0.231 / 0.421 /
  -0.229). Professional 02-08-2026: Ignores `-yes` (stays 0 / 0.5 / 0 for
  every kinematics and index). Status: Deliberate divergence, project owner:
  follow the manual. Also Pro-suspect. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **default memory with `materi_displacement`** (C18): RESOLVED in F3
  (2026-10-07), `-total_linear` as in the Professional.
- **Pro `hypo12` target after the default-memory change** (known difference,
  not decided; BACKLOG C22). GNU: sigyy misses the target by 1.02e-3
  (tolerance 1e-3). Professional 02-08-2026: Meets its own target. Status:
  Open, LOW; revisit with the Masin kernel calibration. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **`-updated` large-rotation robustness (`update2`)** (known difference,
  not decided; BACKLOG C39). GNU: "solution seems to diverge" (veliy -86434)
  with p-refinement. Professional 02-08-2026: veliy 12.243 (target 12.24).
  Status: Open, LOW; not investigated (related to AUDIT K19). Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **`group_materi_memory -updated_without_rotation`** (GNU extension / input
  dialect; AUDIT 1.1). GNU: GNU-only memory (legacy). Professional
  02-08-2026: Not in manual 6.686; the Pro's own copies of these tests use
  `-updated_linear` or `-total_linear`. Status: GNU extension; kept. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **`group_materi_memory -total_piola`** (GNU extension / input dialect;
  manual pages only). GNU: Accepted (2014 source; total formulation).
  Professional 02-08-2026: Not listed in manual 6.686; binary behaviour NOT
  measured. Status: GNU extension; Pro behaviour undetermined. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
