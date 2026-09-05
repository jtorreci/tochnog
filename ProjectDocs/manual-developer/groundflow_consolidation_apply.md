# groundflow_consolidation_apply

## Files and functions

- `groundfl.cc` — in `groundflow()` (lines 49-67). The switch is read after the
  existing `GROUP_GROUNDFLOW_MATERIDIVERGENCE` and applied on top of it:
  ```c
  db( GROUNDFLOW_CONSOLIDATION_APPLY, 0, &groundflow_consolidation_apply,
    ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  ...
  if ( groundflow_consolidation_apply==-NO ) materidivergence = -NO;
  ```
  `materidivergence==-NO` skips the divergence block (lines 83-100), which is
  the consolidation coupling between `materi_velocity` divergence and the
  pressure equation.
- `database.cc` — keyword registration (after `GROUNDFLOW_ADDTOPRESSURE`):
  type `INTEGER`, `data_length = 1`, `no_index = 1`,
  `data_class = GROUNDFLOW`.
- `check.cc` — requires `groundflow_pressure` to be initialized.
- Enum `GROUNDFLOW_CONSOLIDATION_APPLY` in `tochnog.h` (alphabetical, between
  `GROUNDFLOW_ADDTOPRESSURE` and `GROUNDFLOW_DENSITY`).

## Implementation details

- `no_index`: the database accessor is called with the fixed index `0`.
- Read with `GET_IF_EXISTS`: when the record is absent, the variable keeps its
  default `-NO` and nothing changes. The default is `-NO` (aligned with the
  Professional manual 6.556; the GNU legacy initializer was `-YES`, changed in
  the u-p consolidation sprint, commit `0560dcf`).
- Precedence (from lowest to highest): `GROUP_GROUNDFLOW_MATERIDIVERGENCE`,
  `GROUP_GROUNDFLOW_CONSOLIDATION_APPLY` (group), `GROUNDFLOW_CONSOLIDATION_APPLY`
  (global), `CONTROL_GROUNDFLOW_CONSOLIDATION_APPLY` (timestep),
  `OPTIONS_SKIP_GROUNDFLOW_MATERIDIVERGENCE` /
  `CONTROL_OPTIONS_SKIP_GROUNDFLOW_MATERIDIVERGENCE` (legacy skip switches).
  Any `-NO`/`-YES`-skip forces `materidivergence = -NO`; any explicit `-yes`
  record activates the coupling unless a higher-precedence record turns it off.

## Why the default changed (measured)

- Ground14/15/16 of the corpus (safety piping/lifting of an artesian column,
  manual Pro 2.4.1/6.919) do NOT request consolidation. The Professional
  binary (25-10-2023) with the record absent reaches the drained steady state
  inside the 1 s window: sigma'(0) = -10.000 exact at t = 0.5 s (its default
  `groundflow_consolidation_apply` is `-no`, so no skeleton-volume source in
  the pressure equation and no u-p transient).
- The GNU legacy default `-yes` ran the true consolidation transient of the
  model (excess pore pressure dissipating over ~10^2 s: sigma'(0) = -9.53 at
  t = 1 s) and the corpus targets failed.
- With the coupling OFF, the GNU reproduces the Professional drained state to
  8 digits in the same two 0.5 s steps: ground14/15/16 rc = 0.
- Cross-check with the coupling FORCED on both binaries (`-yes`): the GNU and
  the Professional transients agree within a few % (sigma'(0) -9.53 vs -9.29
  at t = 1 s), so the GNU coupling itself was never the problem — only the
  default.

## External dependencies

- Core `db()` accessor only; no external library.
- Globals `materidivergence`, `materi_velocity`, `pres_indx` (in `tochnog.h`).

## Hardcoded parameters / pending refactorings

- The three switches (group/global/timestep) plus the two legacy skip switches
  all funnel into the single `materidivergence` flag; the precedence chain is
  simple but implicit. A small helper would make the intent clearer.
