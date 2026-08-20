# group_groundflow_nonsaturated_vangenuchten

## Files and functions

- `groundda.cc` — in `groundflow_data()` (lines 79-154). When
  `GROUP_GROUNDFLOW_NONSATURATED_VANGENUCHTEN` is active and the
  `groundflow_nonsaturated_apply`/`control_groundflow_nonsaturated_apply`
  switches are `-yes`, the law is evaluated at the integration point:
  ```c
  head = -new_unknowns[pres_indx] / ( dens * gravity );
  tmp  = 1. + pow( vg[2]*fabs(head), vg[4] );
  S    = vg[0] + (vg[1]-vg[0]) * pow( tmp, (1.-vg[4])/vg[4] );
  Se   = (S - vg[0]) / (vg[1]-vg[0]);
  krel = pow( Se, vg[3] ) * pow( 1. - pow( 1. - pow( Se, vg[4]/(vg[4]-1.) ),
    (vg[4]-1.)/vg[4] ), 2. );
  ...
  for ( idim=0; idim<ndim; idim++ ) pe[idim] *= krel;
  dS = (vg[1]-vg[0]) * ((1.-vg[4])/vg[4]) * pow( tmp, (1.-vg[4])/vg[4]-1. )
       * vg[4] * pow( vg[2], vg[4] ) * pow( fabs(head), vg[4]-1. );
  if ( head<0. ) dS = -dS;
  C += por * dS;
  ```
  `vg[0..4] = {Sresidu, Ssat, ga, gl, gn}`; `por` from
  `GROUP_GROUNDFLOW_POROSITY`; `dens` from `GROUNDFLOW_DENSITY`; `gravity` is
  the sum of |components| of `force_gravity_calculate()`.
- `database.cc` — keyword registration: type `DOUBLE_PRECISION`,
  `data_length = 5`, `data_class = GROUNDFLOW`,
  `data_required = GROUP_TYPE`.
- `check.cc` — requires `groundflow_pressure` and `groundflow_saturation`.
- Enum `GROUP_GROUNDFLOW_NONSATURATED_VANGENUCHTEN` in `tochnog.h`.

## Implementation details

- The law modifies `pe[]` (permeability, `krel * ksat,i`) and `C` (capacity,
  `csat + n*dS/dhead`), which are then used by `general()` for the diffusion
  and inertia terms of the `GROUNDFLOW_PRESSURE` equation.
- The saturation is stored in the `groundflow_saturation` dof
  (`node_dof[gsat_indx]`), which is initialised in the initialisation part but
  not solved (no principal dof).
- `gravity` is computed as the sum of the absolute components of
  `force_gravity`, so a 2D gravity `(0,-1)` gives |g|=1 and `head = -pres`.
- Edge cases: `S`, `Se` and `krel` are clamped to `[0,1]`/`[0,...]`;
  `eps_permeability` from `group_groundflow_nonsaturated_eps_permeability`
  floors `krel`; the head division is skipped when `dens*gravity==0`.

## External dependencies

- Core `db()` accessor, `get_group_data()`, `force_gravity_calculate()`.
- Globals `groundflow_saturation`, `gsat_indx`, `pres_indx`, `ndim`.

## Hardcoded parameters / pending refactorings

- The van Genuchten `dS/dhead` derivative is implemented analytically; the
  manual derives the non-saturated capacity as `c = csat + n*dS/dphi_p`.
- `gravity` as sum of |components| is a simplification; for rotated gravity
  vectors the vertical component would be more correct.
