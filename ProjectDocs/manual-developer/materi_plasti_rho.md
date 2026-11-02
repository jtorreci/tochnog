# materi_plasti_rho — developer notes

## Where

- `input.cc`: `materi_plasti_rho = 1`, `rho_indx`, 6 components
  (`dof_type = -MATERI_PLASTI_RHO`, `n*nder` slots).
- `materi.cc` `materi()`: `old_rho[idim*MDIM+jdim] =
  old_unknowns[rho_indx + stress_indx(idim,jdim)*nder]`; passed to
  `set_stress()`.
- `stress.cc` `set_stress()`: `inc_rho = plasti_kinematic_hardening *
  inc_epp`, `new_rho = old_rho + inc_rho`; the yield test uses
  `test_sig = new_sig - new_rho`.
- `general.cc` `general()`: `inertia = 1`, `conv_part = 1` (kappa
  pattern) and, since C68, the lumped denominator with
  `strain_recovery_weight()`.

## Integration-point state (C68, F4 T6f, 2026-10-12)

Before C68 the ELEMENT_INTPNT_DOF masks of `elem.cc` omitted the rho block, so
`old_rho` was the h-lumped nodal value interpolated back to the point
(the kappa defect of c3e6d2a). Now:

- `elem.cc` `elem()`: read and write masks include
  `materi_plasti_rho && rho_indx <= i < rho_indx + 6*nder`.
- `materi.cc` `materi()`: the C68 block (right after the kappa block,
  outside the `|h| > EPS_H` guard) assembles
  `volume * w * (new_rho - old_rho) / dtime` per component
  (`k = stress_indx(idim,jdim)`, row `rho_indx/nder + k`) with
  `w = strain_recovery_weight()`, and writes
  `new_unknowns[rho_indx + k*nder] = new_rho` for ELEMENT_INTPNT_DOF. The old
  h-weighted block was removed.

Verification: `test-2026/rho_ip_state1` (four tria3, membrane von Mises
+ kinematic hardening, targets = Professional 02-08-2026): previous
binary rc 1 (reactions 1.1e-1 relative off), fix rc 0 (7e-6); nodal
rhoyy within 2e-3.

## Hardcoded / pending

- Pre-existing, not changed in C68 (BACKLOG C68d): in `materi.cc` (old-
  state read loop) `old_rho[indx]`, `old_epi[indx]` and `old_epp[indx]`
  use `indx` after the Maxwell loop has overwritten it with
  `m*MDIM*MDIM + idim*MDIM + jdim`. With one chain this equals the
  intended index; with two or more chains (and rho, epi or epp
  initialised) the reads go to the wrong slots, beyond the 9-entry
  arrays. Untested combination.
