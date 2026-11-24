# Bonded-soil generalised Cam clay — developer

Records: `group_materi_plasti_generalised_non_associate_cam_clay_for_bonded_soils` (6.700),
`materi_plasti_generalised_non_associate_cam_clay_for_bonded_soils_history` (4.21),
`materi_strain_plasti_generalised_non_associate_cam_clay_for_bonded_soils` (4.38). User page of the
same name has the model.

## Where

- Enums (`tochnog.h`, `tochnog-mod.h`): `GROUP_MATERI_PLASTI_GENCAM_BONDED`,
  `MATERI_PLASTI_GENCAM_BONDED_HISTORY`, `MATERI_STRAIN_PLASTI_GENCAM_BONDED` (short enum names; the
  record names are 71 to 73 characters, within `MCHAR` 100). Flags
  `materi_plasti_gencam_bonded_history`, `materi_strain_plasti_gencam_bonded`, index
  `gencamepp_indx` (defaults in `initia.cc`).
- `database.cc`: the group record (DOUBLE, 11 values, class MATERI), the two initialisation names,
  and the dof labels `gencamhis0/1` and `eppgencamxx..zz` in `db_initialize`.
- `input.cc`: the history joins the shared `materi_history_variables` branch (2 variables, dof type
  `MATERI_PLASTI_GENCAM_BONDED_HISTORY` at `hisv_indx`) and is exclusive with every other history
  declaration; the plastic strain sets `gencamepp_indx`; the history type is an element dof.
- `general.cc` (inertia, recovery, conv_part), `so.cc` (no fast solver), `check.cc` (the group needs
  `materi_stress`, the history and `materi_strain_plasti`).
- `materi.cc`: the per-model plastic strain arrays grew from 9 to 10 entries (`pm_count`).
- `plasti.cc`: helpers `gencam_eta_limit`, `gencam_integral`, `gencam_eta_of_x`, `gencam_yield`,
  `gencam_potential`, and a model block in `plasti_rule()` after the Cam clay block.

## Design notes

- **Closed form of I(eta).** With the real roots r1, r2 of Q, partial fractions give
  I = (r1 ln(1-eta/r1) - r2 ln(1-eta/r2)) / (A (r1-r2)), A = 1-mu. The roots come from
  qr = -(b + sgn(b) sqrt(D))/2, r1 = qr/A, r2 = c/qr, with `log1p`, so the form stays accurate when
  mu -> 1 (r1 -> infinity) and alpha -> 0 (r2 -> 0). Complex roots use the arctan form; alpha = 0 and
  mu = 1 have their own forms. An earlier (1/(2A)) ln Q form lost 4 digits at mu = 1 - 1e-7.
- **Yield function.** `p_c (ln(p*/p_c) + I(eta))` for p* >= 0.05 p_c, and `q - p* eta(p*/p_c)` near
  and beyond the apex, with eta(x) the inverse of I by safeguarded Newton. The ln form alone tends to
  -infinity at the apex, so isotropic tension jumped back into compression. The q form alone made the
  stress.cc secant iteration oscillate at the critical state, and the |deps_v^p| noise then destroyed
  the bonding (p_m 48.8 -> 47.1 over 0.3 strain). Both forms have the same sign and the same zero
  set, so the switch at 0.05 p_c does not move the surface.
- **Potential.** `ln(p*) + I_g(eta)` (times p_c): its level sets are the potential surfaces of every
  size, so its gradient is the flow direction (the size p_g never has to be found). At and beyond the
  apex the smooth `q^2/(2 C p_c) - C p*` gives purely volumetric flow under isotropic tension; the
  cone `q - C p*` gave spurious shear flow at the vertex.
- **Hardening.** p_s and p_m are updated inside `GET_YIELD_RULE` from the current `inc_epp`
  (exponential, exact for constant rates) and written to `new_hisv` only there, as the Cam clay block
  does; `GET_FLOW_RULE` reads the old values.

## Gotchas

- The flow direction is evaluated once, at the elastic trial stress (stress.cc), for every model.
  The error is first order in the step (critical stress ratio 1.125 at strain steps of 5e-4, 1.195 at
  5e-5, M_g = 1.2).
- Do not write `new_hisv` in the `GET_FLOW_RULE` branch: `plasti_rule()` runs every model block for
  each of the 18 central-difference calls.
