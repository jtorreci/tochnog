# groundflow_total_pressure_limit (developer)

## Files and functions

- `tochnog.h` / `tochnog-mod.h` — enum `GROUNDFLOW_TOTAL_PRESSURE_LIMIT`
  (after `GROUNDFLOW_SEEPAGE_NODE`; both headers in sync) and declaration
  of `double groundflow_total_pressure_limit( void )`.
- `database.cc` — keyword registration: `DOUBLE_PRECISION`,
  `data_length = 1`, `no_index = 1`, `data_class = GROUNDFLOW`.
- `check.cc` — requires the `groundflow_pressure` unknown.
- `groundfl.cc`:
  - `groundflow_total_pressure_limit()` — single source of truth of the
    limit: the record value, or `GROUNDFLOW_TOTAL_PRESSURE_LIMIT_DEFAULT`
    (0., manual 6.587) when the record is absent.
  - `groundflow_phreatic_coord()` — caps the calculated `total_pressure`
    at the limit (both cap points: after the phreatic-level branches and
    after the `post_calcul_static_pressure_height` branch). Every consumer
    of the total pressure goes through this function: `post_calcul
    -groundflow_pressure -total_pressure` (calcul.cc), the pore pressure on
    the skeleton (materi.cc), the wet/dry density, `-total_pressure`
    boundary conversions (bounda.cc), area.cc and group.cc.
    `groundflow_pressure_atmospheric` (GNU-only record, default 0) now
    caps the STATIC pressure only.
  - `groundflow()` — dry rule: when `groundflow_velocity` and
    `groundflow_pressure` are unknowns, the integration-point total
    pressure is interpolated from the nodal total pressures
    (`groundflow_phreatic_coord()` on `NODE` VERSION_NORMAL coordinates and
    `NODE_DOF` VERSION_NEW, i.e. the capped nodal values, as the
    Professional `element_total_pressure` output shows). When it is
    `>= limit - EPS_TOTAL_PRESSURE_LIMIT*(1+|limit|)` the velocity source
    `pe*grad(pres)` of that integration point is replaced by 0, so the
    projected `gvel` dof is 0 there. The pressure equation itself
    (permeability and consolidation terms) is not touched.
- `top.cc` — no post-solve hook any more (the former
  `groundflow_total_pressure_limit_apply()` clamp of the `pres` dof was
  removed).

## Hardcoded parameters

- `GROUNDFLOW_TOTAL_PRESSURE_LIMIT_DEFAULT` = 0. (manual 6.587).
- `EPS_TOTAL_PRESSURE_LIMIT` = 1.e-10 (relative tolerance of the
  "limit reached" test; the interpolation of capped values that all equal
  the limit can differ from it by rounding).

## Semantics (measured on the Professional 02-08-2026, F4 T1)

Probe: one quad4, prescribed solid velocity with divergence 2,
`groundflow_consolidation_apply -yes`, `inertia_apply -yes`, pres 0
prescribed at the bottom, free top nodes (uncapped pres 1.6667), plus the
same with the opposite velocity (pres -1.6667), each with no record and
limits 0 / 0.5 / 100 / -0.5:

- `pres` dof: never cut (1.6667 or -1.6667 in every case).
- nodal total pressure: `min(pres + static, limit)` at every node,
  prescribed ones included (limit -0.5: -0.5 at the prescribed bottom).
- integration-point total pressure (`element_total_pressure`): the
  interpolation of the capped nodal values (limit 0.5: 0.10566 / 0.39434).
- groundflow velocity: 0 when the integration-point total pressure is
  `>=` the limit (no record, 0, -0.5 on the positive field), otherwise
  `k grad(pres)` from the uncapped `pres` (limit 0.5: 0.16667, not the
  gradient of the capped total pressure).
- consolidation: never skipped (the free-node `pres` develops with the
  default limit 0), although manual 2.4.3 says it is skipped when the
  limit is 0 and the total pressure is 0. The GNU follows the binary; the
  former GNU dry rule (limit 0 + pres 0 -> `materidivergence = -NO`) was
  removed, because with the default limit 0 it would switch the
  consolidation off in every model that starts at pres 0.

All 8 probes (positive/negative field, no record / 0 / 0.5 / 100 / -0.5)
and `ground1` (no record / 1.e10) give the Professional values with the
GNU after F4 T1.

## Verification

- `test-2026/groundflow_total_pressure_limit_default1`, `_high1`,
  `_negative1` (pres, total pressure via `-node_dof_calcul`, gvely).
- `test-2026/groundflow_velocity_dry1`, `_dry2` (`ground1`, K10).
- `test-2014/groundflow_total_pressure_limit` and `_dry` retargeted to the
  Professional values (pres 1.6667 in both; gvely 0.16667 / 0). Old
  targets: pres 0.5 / 0 (GNU dof clamp and dry consolidation skip).

## Pending refactorings

- `groundflow_phreatic_coord()` still adds `groundflow_addtopressure`
  twice when both cap blocks run (pre-existing; the record is GNU-only and
  unused by every corpus).
- The dry test calls `groundflow_phreatic_coord()` for every node of the
  element at every integration point; caching the nodal total pressures
  per element would avoid the repeated phreatic-level lookups.

## Differences with the Professional version

None known (former DIVERGENCES TPL and K10, closed in F4 T1 by following
manual 6.587 and the Professional after the project-owner decision of
2026-10-07).

## Related

- `groundflow_consolidation_apply` family (switches of the consolidation
  term; independent of this limit).
- `groundflow_pressure_atmospheric` (GNU-only cap of the static pressure).
- [groundflow_velocity_sign](groundflow_velocity_sign.md).

- **Deliberate divergence from the manual (decision 2026-10-09, the "ley seca"):** the dry rule of manual 2.4.3 (limit 0 and pres 0 skip the consolidation term) is not applied, as in the Professional 02-08-2026. See `ProjectDocs/DIVERGENCES-PROFESSIONAL.md` (DRY-RULE).
