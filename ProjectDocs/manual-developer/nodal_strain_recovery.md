# nodal_strain_recovery — developer notes (K11, 2026-10-06)

## Where

- `miscel.cc` `strain_recovery_weight()`: `sri_stress_recovery_weight(0,
  ...)` (the Professional per-element least-squares weight F/volume),
  or `h` for SRI elements. Declared in `tochnog.h` and `tochnog-mod.h`.
- `materi.cc` `materi()`: the `materi_strain_elasti`, `_plasti` and
  `_total` numerators (`volume * w * inc_ep? / dtime`) moved out of the
  `h > EPS_H` block into one block right after the stress recovery
  (after the `stress_recovery_post_pass` skip), with `w =
  strain_recovery_weight()` and the guard `|w| > EPS_H || |h| > EPS_H`
  (the least-squares weight is non-zero where h vanishes, as for the
  stress). The plastic heat generation stays in the old block.
- `general.cc` `general()`: the lumped denominator and inertia term of
  `MATERI_STRAIN_TOTAL/_ELASTI/_PLASTI` use the same weight, so the
  incremental update converges to `old + sum(F inc)/sum(F)` per element
  and the arithmetic mean over the elements, as the stress.

## Verification

- `tslv_part` (validation-suite, converted, Professional 02-08-2026):
  every strain within 4e-11 of the Professional (before: clamped-edge
  eptyy -2.35e-6 vs 0, eptxy a factor 3 low).
- Pro `tension_direct1`: the eptyy residual 0.667 vs 0.5 at the shared
  nodes is gone (every variable = Professional; with
  `materi_strain_plasti` also eppyy).
- Test `test-2026/strain_recovery_clamped1`: RED (rc 1) on the old
  binary, GREEN new, the Professional passes it.
- Corpus A/B (k13 -> k11 build, audit runner): 46 tests change, only
  strain dofs and their derivative slots (`tept`/`xept`/`yept`, `epp` in
  `mhardsoil_plast`, `examp16`, `ho_mech2`); no classification change.
  Against the Professional on the converted inputs it runs (26 of 44):
  14 tests move to the Professional values (msf_*, qsri_*, tslv_*,
  gforce7q4, mesh_act_grav, cmat_gate, element_normal2d/3d_column1 eptzz,
  examp16, force3 eptzz); the rest differ only in round-off of zero
  fields (1e-15 / 1e-9) or keep a pre-existing difference of another
  origin (beam2d family, examp17).

## Hardcoded / pending

- kappa (`materi_plasti_kappa`) uses the same weight since F4 T6f WU2
  (2026-10-11), with its integration-point value kept in ELEMENT_DOF
  (`materi_plasti_kappa.md`); validation_3 nodal kap within 2.7 % of the
  Professional (quad4 operator, C42/C64).
- C68 (F4 T6f, 2026-10-12): kappa_shear (`kapsh`), the cap1 pc, the
  damage `dam`, the kinematic hardening `rho` (6 components) and the
  Maxwell stresses `msig` (6 per chain) follow kappa: ELEMENT_DOF masks
  in `elem.cc` (both), one `materi.cc` block right after the kappa block
  (outside the h guard, recovery weight numerator, IP write
  `new_unknowns[<indx>+k*nder]`), and the `general.cc` denominator.
  Verified on four-tria3 probes against the Professional (tests
  `test-2026/{damage,kappa_shear,rho,maxwell}_ip_state1`); cap1 by the
  manual only (BACKLOG C68b).
- The other integration-point dofs (f, history variables, softening
  variables, model-specific plastic strains `epphs*` ...) keep the h
  weight; not measured against the Professional.
- K02 (one-iteration lag of the non-principal dofs when the stress
  recovery is not absolute, e.g. `materi_strain_total` + one equilibrium
  iteration) is not addressed: see BACKLOG C35.
