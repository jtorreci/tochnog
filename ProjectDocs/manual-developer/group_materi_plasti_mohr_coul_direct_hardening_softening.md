# group_materi_plasti_mohr_coul_direct_hardening_softening

## Implementation

- **Enum**: `GROUP_MATERI_PLASTI_MOHR_COUL_DIRECT_HARDENING_SOFTENING`
  in tochnog.h and tochnog-mod.h (in sync; right after
  `GROUP_MATERI_PLASTI_MOHR_COUL_HARDENING_SOFTENING`).
- **Registration** (database.cc): DOUBLE_PRECISION, length 7,
  MATERI class, requires GROUP_TYPE.
- **Consumption** (plasti.cc, `plasti_rule`): the hardening-softening
  branch now accepts the non-direct record OR the direct record
  (`get_group_data` of the direct record sets the local
  `mc_hs_direct` flag). When the direct variant is used, the angles
  (indices 0, 2, 3, 5 = phi_0, phiflow_0, phi_1, phiflow_1) are
  converted with `plasti_data[i] *= PIRAD/180.` BEFORE the interpolation
  block (the cohesion entries are not angles).
- **Hardening variable**: the branch reads `new_unknowns[kapsh_indx]`
  when `materi_plasti_kappa_shear` is declared, else
  `new_unknowns[kap_indx]` (manual 6.731: the hardening variable is the
  shear plastic strain kappa_shear).

## Verification

- dam_building: parses and runs with the record
  `0. 16.5 0. 0. 12. 0. 0.18` (angles 0 deg, cohesion 16.5 -> 12 kPa).
  The full model exceeds the corpus time budget (~10 minutes vs 45 s),
  so the rc=0 check is not reachable in the corpus; the record, the
  degree conversion and the kappa_shear hardening are exercised.

## Pending

- The softening branch is exercised only by dam_building in the corpus;
  a dedicated small test (analytical calibration like the non-direct
  variant, SEGUIMIENTO sprint 10 lote 4) is pending.
