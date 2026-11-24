# group_materi_plasti_hypo_niemunis_visco (developer)

Also `group_materi_plasti_hypo_niemunis_visco_ocr` and
`control_materi_plasti_hypo_niemunis_visco_ocr_apply`.

## Where

- `database.cc`: `GROUP_MATERI_PLASTI_HYPO_NIEMUNIS_VISCO` (DOUBLE, 9,
  MATERI) and `_OCR` (DOUBLE, 1, MATERI); enums after
  `GROUP_MATERI_PLASTI_HYPO_WOLFFERSDORFF_NIEMUNIS` in `tochnog.h` /
  `tochnog-mod.h`. The control record existed already. `check.cc` requires
  `materi_stress`, the hypo history and `materi_strain_intergranular`.
- `hypoplas.cc` `hypoplasticity()`: dispatch after wolffersdorff /
  lowangles, `hypo_type` 2, data checks (db_error, b2 < 0 error), no
  pressure dependent void ratio. Initial void ratio from the OCR when
  `step_first_of_control` (set in `top.cc` at the first step of a control
  index) and the control switch of the current index is -yes:
  `niemunis_visco_peplus()` on the old stress.
- `hypo.c`: `niemunis_visco_shape()` (That, That*, F, pe+),
  `niemunis_visco_peplus()`, `niemunis_visco_rate()` (L, flow rule, OCR,
  creep rate, M with intergranular strain, stress update, `his[5]` OCR,
  stiffness size). `sigma_()` calls it for `ihypotype == 2` after the void
  ratio update and jumps to the mobilised friction angle.
  The `hypo_()` driver is shared: substepping, midpoint intergranular strain
  evolution with the total D (the thesis does not use D - Dvis).

## Implementation details

- Explicit Euler per substep with L, M and Dvis at the start of the substep.
- Loading test of M: h:D > 0 (thesis); the von Wolffersdorff branch uses
  >= 0 for the L:hh term.
- pe+: betaR within 1e-8 of 1 uses the limit p (1 + eta^2).
- Tension (trT >= 0) flags a substep error (substep halving).
- The manual flow rule sign is garbled; the form above follows Niemunis 2003
  eq. 4.73 (isotropic stress gives a compressive creep).
- The b2 term is a Tochnog extension (not in the thesis).
- The old `masin_niemunis_visco_umat` (masin_visco.c) is the Masin clay
  route and is not used here.

## Verification

Professional other/hypo6 and 20 Professional probes (isotropic and
anisotropic relaxation, betaR 0.5 / 1.5, phi 25, nu 0.3, Dr 10, Iv 0.1,
oedometric compression and extension with intergranular strain, OCR
application at index 20 and 30): agreement within 5e-5 relative. An
independent Python implementation of the equations (explicit Euler) equals
the GNU to 1e-8 for relaxation. The isotropic relaxation equals the Euler
recursion of dp/dt = -sqrt3 Dr/kappa p (p/pe)^(1/Iv).

## Differences with the Professional version

- b2 < 0 rejected (class D NIEMUNIS-VISCO-B2).
- `hyhis3` content (class B C120-NIEMUNIS-VISCO-HYHIS3).

## Pending

- No consistent (implicit) tangent as in Niemunis 2003 eq. 4.128-4.134.
