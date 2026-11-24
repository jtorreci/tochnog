# materi_plasti_kappa — developer notes

## Where

- `materi.cc` `materi()`: `old_kappa = old_unknowns[kap_indx]`; the
  constitutive update (`set_stress`, `plasti.cc`) returns `new_kappa`.
  The kappa block of the nodal assembly (F4 T6f WU2, 2026-10-11) sits
  next to the strain recovery, outside the `|h| > EPS_H` guard:
  numerator `volume * w * (new_kappa - old_kappa) / dtime` with
  `w = strain_recovery_weight()` (the Professional per-element
  least-squares weight F/volume; h for SRI elements), and
  `new_unknowns[kap_indx] = new_kappa` when `options_element_dof -yes`
  (the default) so the integration-point value reaches ELEMENT_INTPNT_DOF.
- `elem.cc` `elem()`: the ELEMENT_INTPNT_DOF read mask (restore of the
  remembered integration-point dofs into `old_unknowns`) and the write
  mask (`element_dof_npoint = new_unknowns`) include
  `(materi_plasti_kappa && i==kap_indx)`.
- `general.cc` `general()`: the lumped denominator and inertia term of
  `MATERI_PLASTI_KAPPA` use `strain_recovery_weight()` (as the strain
  tensors, K11), so numerator and denominator match.

## Defect fixed (F4 T6f WU2)

The ELEMENT_INTPNT_DOF masks omitted `kap_indx`, so `old_kappa` was the h-lumped
NODAL kappa interpolated back to the point (`old_unknowns` holds the
nodal interpolation for every dof outside the masks). The hardening then
read a smeared, lagging value: `validation_3` displacements up to 1.25 %
(of the field maximum) off the Professional, nodal kap 11.6 % off. Fix A
(masks + IP write) brings the solution to the Professional
(displacements within 0.06 %, post_node_result vely -129.353 vs
-129.302, HEAD -128.835); fix B (the recovery weight of the nodal
output) brings the nodal kap within 2.7 % (fix A alone: 16 %, the
h-weighted nodal average of IP values).

Verification:
- `test-2026/kappa_ip_state1` (two quad4, von Mises + kap hardening):
  HEAD rc 1 (node 1 reaction -11.85835 vs Pro -11.91214), fix rc 0
  (-11.91338); nodal kap node 1 4.555e-3 vs Pro 4.635e-3 (HEAD 3.463e-3).
  Fix A alone passes the reactions but not the nodal kap.
- `test-2014/validation_3_kappa`: Pro targets post_node_result vely,
  disx/disy, kap node 1; HEAD rc 1, fix rc 0.

## Hardcoded / pending

- `kapsh_indx`, `cap1_indx`, the `rho` block, `dam_indx` and the Maxwell
  `mstres_indx` block had the same defect; fixed with the same pattern in
  C68 (F4 T6f, 2026-10-12): one block right after this kappa block (see
  `nodal_strain_recovery.md`).
- The quad4 nodal value keeps the GNU bilinear operator (C42/C64).
