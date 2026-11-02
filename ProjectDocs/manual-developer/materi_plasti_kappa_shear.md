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
- **Element RHS** (materi.cc, `element_rhside`): since C68 (F4 T6f,
  2026-10-12) `volume * w * (new_kapsh - old_kapsh) / dtime` with
  `w = strain_recovery_weight()`, outside the h guard, and the IP value
  is written to ELEMENT_INTPNT_DOF (same as kap; see "Integration-point state
  (C68)" below).
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

## Integration-point state (C68, F4 T6f, 2026-10-12)

`materi_plasti_kappa_shear` is an integration-point history state (Professional manual
4.25), like the hardening parameter kappa (`materi_plasti_kappa.md`).

- `elem.cc` `elem()`: the ELEMENT_INTPNT_DOF read and write masks include
  `(materi_plasti_kappa_shear && i==kapsh_indx)`, so `old_kapsh` in `materi.cc` is the value remembered at the
  integration point, not the h-lumped nodal value interpolated back.
- `materi.cc` `materi()`: one block right after the kappa block (after
  the `stress_recovery_post_pass` skip, outside the `|h| > EPS_H`
  guard) assembles the nodal (output) numerator
  `volume * w * (new - old) / dtime` with `w = strain_recovery_weight()`
  and writes the IP value to `new_unknowns[<indx> + k*nder]` when
  `options_element_dof -yes`. The old h-weighted block was removed.
- `general.cc` `general()`: the lumped denominator of `MATERI_PLASTI_KAPPA_SHEAR` uses the
  same weight.

Verification: `test-2026/kappa_shear_ip_state1` (four tria3, von Mises
hardening on `-kapsh`, targets = Professional 02-08-2026): previous binary
rc 1 (1.3e-2 relative off), fix rc 0 (1.5e-4). The Mohr-Coulomb
hardening-softening use of kapsh keeps an open 1.2x difference on a
homogeneous control (BACKLOG C68c, not an IP-state issue).

## Differences with the Professional version

- **Mohr-Coulomb hardening-softening on kapsh** (C68c, class B, cause unknown): 1.2x between GNU and the Professional on a homogeneous control. Index: [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
