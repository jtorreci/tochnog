# group_groundflow_total_pressure_tension

## Files and functions

- `materi.cc` — in `materi()` inside the `if ( groundflow_pressure )` block
  (lines 420-452). After `new_pres` is obtained (unknown or phreatic-based),
  the correction is applied before the pressure factor:
  ```c
  if ( db_active_index( GROUP_GROUNDFLOW_TOTAL_PRESSURE_TENSION, gr,
      VERSION_NORMAL ) ) {
    double gptt[2], epp_princ[3], epp_t[MDIM*MDIM], dens_water=0.;
    db( GROUP_GROUNDFLOW_TOTAL_PRESSURE_TENSION, gr, idum, gptt, ldum,
      VERSION_NORMAL, GET );
    db( GROUNDFLOW_DENSITY, 0, idum, &dens_water, ldum, VERSION_NORMAL,
      GET_IF_EXISTS );
    // epp_t from new_unknowns[epp_indx+stress_indx(idim,jdim)*nder]
    matrix_eigenvalues( epp_t, epp_princ );
    tmp = max(epp_princ[0..2]);
    if ( tmp > gptt[0] ) {
      double static_pres =
        force_gravity[ndim-1] * dens_water * ( gptt[1] - coord_ip[ndim-1] );
      if ( scalar_dabs(static_pres) > scalar_dabs(new_pres) ) new_pres = static_pres;
    }
  }
  new_pres *= gpf;
  for ( idim=0; idim<MDIM; idim++ ) total_new_sig[idim*MDIM+idim] += new_pres;
  ```
- `database.cc` — keyword registration (after `GROUP_GROUNDFLOW_POROSITY`):
  type `DOUBLE_PRECISION`, `data_length = 2`, `data_class = GROUNDFLOW`,
  `data_required = GROUP_TYPE`.
- `check.cc` — requires `groundflow_pressure` and `materi_strain_plasti`.
- Enum `GROUP_GROUNDFLOW_TOTAL_PRESSURE_TENSION` in `tochnog.h`.

## Implementation details

- Uses the water density `GROUNDFLOW_DENSITY` (NOT the material density) for
  the static pressure — the material density would be 0 in groundflow-only
  groups.
- `matrix_eigenvalues()` (math.cc) returns the 3 eigenvalues unsorted; the
  maximum is taken explicitly.
- The correction replaces `new_pres` BEFORE `new_pres *= gpf`, so the pressure
  factor still applies on top.
- The correction affects the total stress used for the nodal forces; the
  stored effective stress (`new_sig`) is unchanged.

## External dependencies

- `matrix_eigenvalues()` (math.cc), `scalar_dabs()`, `force_gravity_calculate()`.
- Globals `materi_strain_plasti`, `epp_indx`, `pres_indx`, `coord_ip`, `ndim`.

## Hardcoded parameters / pending refactorings

- `plastic_tension_minimum` is compared against the LARGEST eigenvalue (tensile
  crack); the manual's wording matches.
