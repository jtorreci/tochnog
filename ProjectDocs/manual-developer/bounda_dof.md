# bounda_dof

## Files and functions

- `bounda.cc` — within `bounda()`: `bounda_dof` is implemented as a functional
  alias of `bounda_unknown`.
  - Local `max_bounda_dof=0` (line 29), added to the `max_bounda` calculation:
    `db_max_index( BOUNDA_DOF, max_bounda_dof, ... )` (line 64) and
    `if ( max_bounda_dof>max_bounda ) max_bounda = max_bounda_dof;` (line 83).
  - The iboun loop treats both records as "unknown" (Dirichlet):
    `unknown = db_active_index( BOUNDA_UNKNOWN, iboun, ... ) ||
               db_active_index( BOUNDA_DOF, iboun, ... );` (lines 86–87).
  - Value read (lines 183–188): uses `BOUNDA_DOF` when it is active, otherwise
    falls back to `BOUNDA_UNKNOWN`.
- `database.cc` — keyword registration (lines 188–192): name `bounda_dof`,
  `type = INTEGER`, `data_length = MBOUNDA`, `fixed_length = 0`,
  `data_class = BOUNDA`.
- `tochnog.h` — enum `BOUNDA_DOF` (line 142), before `BOUNDA_UNKNOWN` (line 154).
- `tochnog-mod.h` — mirror enum `BOUNDA_DOF` (line 135), must stay in sync.

## Implementation details

- `bounda_dof` is a pure alias: the record is read into the same `val` array
  and processed by the identical code path as `bounda_unknown`, so it supports
  the same selection modes — node ranges (`-range`), `-all`, `-node_set` and
  geometries (`-geometry_line`, `-geometry_point`, ...) — through the shared
  `val[0]` dispatch in the node loop (lines 288–305).
- The `unknown` flag drives both the value read and the dof-label lookup; only
  one of `bounda_dof`/`bounda_unknown` or `bounda_force` may be active per
  boundary index.
- The dofs prescribed are the primary dofs (e.g. `-velx`, `-temp`, `-pres`)
  resolved through `dof_label`/`array_member`, identical to `bounda_unknown`.
- Values over time come from the matching `bounda_time` record, as for
  `bounda_unknown`.

## External dependencies

None. Reuses the existing `bounda_unknown` code path and the `geometry()`
routine for geometry-based selection.

## Hardcoded parameters / pending refactorings

- The dof lookup error `db_error( BOUNDA_UNKNOWN, iboun );` fires even when
  the input used `bounda_dof` (messages may confuse users of the
  professional keyword).
- The alias duplicates the `bounda_unknown` registration pattern in
  `database.cc`; the two records could share a single registration helper to
  avoid drift between `BOUNDA_DOF` and `BOUNDA_UNKNOWN`.
- `data_required` is unset for `BOUNDA_DOF`, so `bounda_dof` works standalone
  (like `bounda_unknown`).

## bounda_dof -topres (total pore pressure prescription)

- `db_number("topres")` resolves to the `GROUNDFLOW_PRESSURE` keyword enum
  (NOT to the dynamic dof label that `-pres` resolves to), so a `-topres`
  token stored in the BOUNDA_DOF record does not match `dof_label[]`.
- `bounda.cc` `bounda()`: the per-dof resolution loop special-cases
  `val[iu] == -GROUNDFLOW_PRESSURE` when `groundflow_pressure` is active:
  it maps the token onto the pres dof (`iuknwn = pres_indx`) and sets the
  per-record flag `topres_bounda` (reset per iboun next to `rotate`).
- The load application (`else` branch of the value setting) converts the
  prescribed total pressure `load` into the head value per node by INVERTING
  `groundflow_phreatic_coord()`:
  - `found` (phreatic level / static height covers the node):
    `new_node_dof = load - static_pressure - addtopressure`;
  - otherwise (no level): `new_node_dof = load + dens*g*z - addtopressure`
    (the no-level total is `pres_dof - rho*g*z`).
  The node coordinate is the current one (NODE + materi displacement).
- GOTCHA: the inversion uses the SAME static the machinery adds, so it is
  exact for whatever covers the node — including a
  `post_calcul_static_pressure_height` region (ground13: bottom/top rows lie
  inside the region height_ref=123, so the prescribed -20/-10 become the
  uniform head +1210, and `-to_pres` = -20/-10 / `-dy_pres` = 1210 match the
  Professional .dbs digit by digit).
- Explicit bounds always win over automatic defaults: `bounda()` zeroes
  `NODE_BOUNDED`, then `groundflow_phreatic_apply()` applies the phreatic
  conditions and finally the bounda_dof records overwrite their nodes.
