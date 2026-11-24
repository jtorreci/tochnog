# groundflow_pressure_factor

## Files and functions

- `materi.cc` — inside the `if ( materi_stress )` block that computes the nodal
  forces (lines 369-380). The factor is read at lines 372-374 and applied at
  line 378:
  ```c
  double gpf = 1.;
  db( GROUNDFLOW_PRESSURE_FACTOR, 0, idum, &gpf, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  new_pres = new_unknowns[pres_indx];
  if ( groundflow_phreatic_coord( -1, coord_ip, new_unknowns,
    total_pressure, static_pressure, location ) ) new_pres = total_pressure;
  new_pres *= gpf;
  for ( idim=0; idim<MDIM; idim++ ) total_new_sig[idim*MDIM+idim] += new_pres;
  ```
- `database.cc:2058-2062` — keyword registration:
  `strcpy(name[GROUNDFLOW_PRESSURE_FACTOR],"groundflow_pressure_factor")`,
  type `DOUBLE_PRECISION`, `data_length = 1`, `no_index = 1`,
  `data_class = GROUNDFLOW`.
- Enum `GROUNDFLOW_PRESSURE_FACTOR` in `tochnog.h:485` and mirror enum in
  `tochnog-mod.h:478` (must stay in sync).

## Implementation details

- The record is `no_index`; the database accessor is called with the fixed
  index `0`.
- It is read with `GET_IF_EXISTS`, so it is inactive when the record is absent
  (`gpf` stays `1`, the default).
- `new_pres` is the pore pressure at the integration point: taken from the
  unknown `new_unknowns[pres_indx]`, or from `total_pressure` when
  `groundflow_phreatic_coord()` returns true (phreatic-level based pressure).
- The factor is applied BEFORE the pore pressure is added to the total stress:
  `new_pres *= gpf` precedes
  `total_new_sig[idim*MDIM+idim] += new_pres`. The effective stress in
  `new_sig` is not modified.
- The coupling matrix entry for the pressure dof (line 434-438, the
  `if ( groundflow_pressure )` branch) is NOT scaled by the factor; only the
  right-hand-side (nodal force) term is affected.

## External dependencies

- Core database `db()` accessor only; no external library.
- `groundflow_phreatic_coord()` (`groundfl.cc`) — provides the phreatic-based
  `total_pressure`/`static_pressure` when applicable.
- Globals `groundflow_pressure`, `pres_indx`, `total_new_sig`,
  `new_unknowns`, `coord_ip` (declared in `tochnog.h`).

## Hardcoded parameters / pending refactorings

- Default `gpf = 1` is hardcoded before the `db()` read; no `data_required`
  dependency is declared in `database.cc`, so an invalid/absent record cannot
  be detected at input-check time.
- The factor only affects the stress/nodal-force contribution; the pressure
  dof coupling term and the `groundflow_pressure_atmospheric` clamp are not
  scaled, which may be surprising for large factors in unsaturated flow.
- Like other `GROUNDFLOW` records, the enum must be mirrored in
  `tochnog-mod.h`; a single shared enum source would remove the risk of drift.
