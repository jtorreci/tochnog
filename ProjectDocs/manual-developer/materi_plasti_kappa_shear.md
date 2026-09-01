# materi_plasti_kappa_shear

## Implementation

- **Enum**: `MATERI_PLASTI_KAPPA_SHEAR` in tochnog.h and tochnog-mod.h
  (in sync; right after `MATERI_PLASTI_KAPPA`).
- **Globals** (initia.cc): flag `materi_plasti_kappa_shear` and index
  `kapsh_indx` (pattern of `kap_indx`).
- **Parser** (input.cc): `materi_plasti_kappa_shear` registers a SCALAR
  dof `-MATERI_PLASTI_KAPPA_SHEAR` with `n = 1` unknown:
  ```c
  materi_plasti_kappa_shear = 1;
  kapsh_indx = unknown_indx;
  n = 1;
  array_set( &dof_type[kapsh_indx], -MATERI_PLASTI_KAPPA_SHEAR, n*nder );
  array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
  ```
- **Registration** (database.cc): name[] of the initia + basename
  `kapsh` in the DOF_LABEL table (next to `kap`).
- **Unknown handling** (general.cc): `unknown_belongs_to_type = 1`,
  `inertia = 1.`, `conv_part = 1.` (same pattern as `materi_plasti_kappa`).
- **Accumulation** (stress.cc, `set_stress`): in both the incremental
  and the converged plastic branches,
  ```c
  tmp = array_inproduct( epp_dev, epp_dev, MDIM*MDIM );   // deviatoric epp
  new_kapsh = old_kapsh + sqrt( 0.5 * tmp );
  ```
  with `epp_dev = inc_epp - mean(inc_epp)*I` (deviatoric part of the
  plastic strain increment).
- **Element RHS** (materi.cc, `element_rhside`): the kapsh dof is
  integrated with the standard transient pattern
  `volume * h[inol] * (new_kapsh - old_kapsh) / dtime` (same as kap).
- **Hardening** (plasti.cc, `plasti_rule`): the Mohr-Coulomb
  hardening-softening laws read the hardening variable from
  `new_unknowns[kapsh_indx]` when `materi_plasti_kappa_shear` is set,
  falling back to `kap` otherwise.
- **set_stress signature**: `old_kapsh`/`new_kapsh` were added to
  `set_stress` (tochnog.h / tochnog-mod.h in sync); the caller
  (materi.cc) passes the dof values (0 when the dof is absent).

## Verification

- dam_building: parses and runs (17-dof reset including `-kapsh`); the
  full layered model needs ~10 minutes, so it stays RUNFAIL in the
  corpus (45 s budget) but the kapsh dof and its reset are exercised.
- slope_classical_numerical / slope_nonlocal_refine: parse and run
  (the kappa_shear initia is consumed).
