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
    clamped; the Professional defaults the limit to 0 and caps the
    calculated total pressure instead, see below).
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

## Differences with the Professional version

- **What is cut and the default** (known difference, pending project-owner
  decision; DIVERGENCES TPL; measured 2026-10-07). Professional 02-08-2026
  (manual 6.587, "Default the limit is set to 0"): reads the record, caps
  the calculated total pressure (`node_dof_calcul` of `post_calcul
  -groundflow_pressure -total_pressure`) at every node including prescribed
  ones, never cuts `NODE_DOF[pres]`, default limit 0; the groundflow
  velocity is 0 where the total pressure reaches the limit (this explains
  AUDIT K10 / `ground1`: limit 1.e10 gives -0.05 = GNU). GNU: record-gated
  clamp of `NODE_DOF[pres]` on free nodes in
  `groundflow_total_pressure_limit_apply()`; dry consolidation skip in
  `groundflow()`; the post-processed total pressure comes from
  `groundflow_phreatic_coord()` and is capped by `pressure_atmospheric`
  (default 0), not by this record. Probe numbers (free-node pres GNU / Pro,
  top total pressure GNU / Pro): no record 1.6667 / 1.6667, 0 / 0;
  limit 0: 0 / 1.6667, 0 / 0; limit 0.5: 0.5 / 1.6667, 0 / 0.5;
  limit 100: 1.6667 / 1.6667, 0 / 1.6667; limit -0.5: -0.5 / 1.6667,
  -0.5 / -0.5 (Pro also -0.5 at the prescribed node). Probe inputs: the
  `groundflow_total_pressure_limit` validation test in Professional
  spelling (`bounda_dof`, `mesh`, `-updated_linear`, no positional
  `node_dof`, plus `post_calcul -groundflow_pressure -total_pressure`,
  without which the Professional stops with "Error: you need to include
  post_calcul -groundflow_pressure -total_pressure"). The earlier
  statement that the Professional ignores the record (AUDIT K09/K16) was
  wrong. Index: [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
