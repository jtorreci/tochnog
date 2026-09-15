# groundflow_total_pressure_limit (developer)

## Files and functions

- `tochnog.h` / `tochnog-mod.h` — enum `GROUNDFLOW_TOTAL_PRESSURE_LIMIT`
  (after `GROUNDFLOW_SEEPAGE_NODE`; both headers in sync, 977 enums) and
  declaration of `groundflow_total_pressure_limit_apply()`.
- `database.cc` (after the seepage registrations) — keyword registration:
  `DOUBLE_PRECISION`, `data_length = 1`, `no_index = 1`,
  `data_class = GROUNDFLOW`.
- `check.cc` — requires the `groundflow_pressure` unknown.
- `groundfl.cc`:
  - `groundflow_total_pressure_limit_apply()` — per-node clamp on
    `NODE_DOF[pres]` (VERSION_NEW). Skips nodes bounded in `pres`
    (Dirichlet); cuts `pres > limit` down to `limit`. Early-returns unless
    the record exists (record-gated: without the record nothing is
    clamped; difference with Professional, which defaults the limit to 0).
  - `groundflow()` — dry check: when the record exists, `limit == 0` and
    the interpolated `new_unknowns[pres_indx]` is 0 (within TINY), the
    element is considered dry and `materidivergence` is forced to `-NO`
    (the consolidation term is skipped for that element). This composes
    with the `groundflow_consolidation_apply` switch family (any -no also
    forces `-NO`).
- `top.cc` — `groundflow_total_pressure_limit_apply()` is called right
  after `solve( options_solver )` inside the equilibrium iteration loop,
  so the clamped values are what the next iteration and the end-of-step
  `VERSION_NEW -> VERSION_NORMAL` copy see.

## Semantics

- Clamp: solved pressures above the limit are cut to the limit. Prescribed
  (bounded) pressure dofs are never cut.
- Dry/consolidation interaction (manual 2.4.3): `limit 0` + pressure 0 →
  "no water" → skip the material divergence term for the element. After a
  clamp to 0 the value is exactly 0.0, so the dry check sees it cleanly.
- The clamp only runs post-solve; the assembled equations of the CURRENT
  iteration still see the unclamped increment. The iteration loop makes
  this self-consistent (next iteration assembles with clamped values).

## Verification

Tests `groundflow_total_pressure_limit` and `_dry`
(validation-suite/test-2014):

- Same mesh as `groundflow_consolidate_off`: prescribed velocity field
  with divergence 2, pres=0 Dirichlet at the bottom, free top nodes.
- **The pressure only develops when the consolidation coupling and the
  capacity storage term are active**, which is not the default since
  `0560dcf` aligned both switches with the Professional (manual 6.557:
  "Default switch is -no"). The tests therefore carry the explicit
  records `groundflow_consolidation_apply -yes` and `inertia_apply -yes`
  (the same convention documented in `groundflow_consolidation_apply.md`
  and used by the corpus consolidation examples `large2`/`large3`). With
  both records the uncapped free-node pressure is `1.6667` (the ~1.67 of
  the manual), i.e. `C·dp/dt = div(v_solid)` at `dt=1`, `C=1`.
  Without them the free nodes stay at 0 and the limit record has nothing
  to clamp — this was the stale reconstruction (see
  `CORPUS-RECONSTRUCTION.md`, fourth pass).
- `groundflow_total_pressure_limit`: `limit 0.5` → free nodes cut to 0.5
  (targets 0.5), Dirichlet nodes unaffected.
- `groundflow_total_pressure_limit_dry`: `limit 0.` → elements are dry,
  the consolidation term is skipped and the pressure stays 0 (targets 0),
  even with consolidation + inertia active.

## Related

- `groundflow_consolidation_apply` family (per-group/timestep switches of
  the same term).
- `groundflow_pressure_atmospheric` (cap on the STATIC pressure derived
  from phreatic levels, GNU-inherited; this keyword caps the SOLVED
  pressure instead).

## Pending refactorings

- None specific. Note the general gotcha (fixed 2026-08-24 in
  interface.cc): `db(..., GET_IF_EXISTS)` does not write `dval` when the
  record is missing — every buffer passed to GET_IF_EXISTS must be
  initialized before the call.
