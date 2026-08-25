# group_materi_elasti_hardsoil

## Implementación

- **Ley elástica**: new block in `set_stress()` in `stress.cc`, right
  after the `GROUP_MATERI_ELASTI_YOUNG_POWER` block. It reads the group
  record `GROUP_MATERI_ELASTI_HARDSOIL` (DOUBLE, fixed length 7) and
  builds C/Cmem with `C_matrix()`. Like young_power, the block CLEARS
  C/Cmem first (C_matrix ACCUMULATES into its target, see the lot 6
  fix): the hardsoil law IS the Young modulus.
- **Keywords** (data_class MATERI) registered in `database.cc`:
  - `group_materi_elasti_hardsoil` (DOUBLE, length 7, required
    GROUP_TYPE): `Eref_50 sigmaref_50 nu50 m Eref_ur sigmaref_ur nuur`.
- **New enums** in `tochnog.h` / `tochnog-mod.h` (kept in sync, same
  order):
  - `GROUP_MATERI_ELASTI_HARDSOIL` (after `GROUP_MATERI_ELASTI_K0`);
  - `GROUP_MATERI_PLASTI_HARDSOIL` (after `GROUP_MATERI_PLASTI_GURSON`);
  - `MATERI_PLASTI_HARDSOIL_HISTORY` (after
    `MATERI_PLASTI_CAP1_HISTORY`);
  - `MATERI_STRAIN_PLASTI_HARDSOIL` (after `MATERI_STRAIN_PLASTI`);
  - `ELEMENT_INTPNT_MATERI_PLASTI_HARDSOIL_GAMMAP_INITIAL` (after
    `ELEMENT_INTERFACE_FORCE_TANG2`);
  - `CONTROL_MATERI_PLASTI_HARDSOIL_GAMMAP_INITIAL` already existed
    (Sprint 9 stub: enum + registration in `database.cc` were present;
    only the behavior was missing).
- **check.cc**: `GROUP_MATERI_ELASTI_HARDSOIL` requires
  `materi_stress`, `materi_velocity` and `materi_plasti_hardsoil_history`
  (the loading/unloading switch needs the maximum |p| history).

## Física

The manual (theory "Hardening-Soil model"):

```
E50 = Eref_50 * ((sig3 + c*cot(phi))/(sigmaref_50 + c*cot(phi)))^m
Eur = Eref_ur * ((sig3 + c*cot(phi))/(sigmaref_ur + c*cot(phi)))^m
```

**Sign mapping (documented decision)**: the manual orders the principal
stresses `sig3 > sig2 > sig1` with `sig1` the LARGEST compressive stress.
In this code (tension positive, compression negative) the largest
compression is the SMALLEST algebraic eigenvalue, so the manual's `sig3`
(the least compressive = the confining stress in triaxial) is the
**LARGEST algebraic eigenvalue** of `new_sig` (`matrix_eigenvalues`).
The base is `sig3_manual + c*cot(phi) = -(largest eigenvalue) +
c*cot(phi)`. NOTE: the task brief said "smallest eigenvalue (more
negative = more compression)" — that is the manual's `sig1` (the AXIAL
stress) and would make the stiffness depend on the axial load, NOT the
HS model; the tests verify the correct (largest-eigenvalue) mapping
analytically. The brief's first sentence ("sig3 is the least compressive
= largest algebraic") is the correct one and was followed.

**Cohesion term**: `c*cot(phi)` is read from the plastic group
`GROUP_MATERI_PLASTI_HARDSOIL`; without it `c = 0` and the term is
skipped (never evaluate `cot(phi)` for `c = 0`). `phi` must be in
`(0, pi/2)` (db_error otherwise, same as cap1).

**Loading/unloading switch**: the maximum |p| history dof
(`materi_plasti_hardsoil_history` = the shared `sph` dof, updated in
`dof.cc` as the running max of |p| over the stress dofs). The decision
uses the step pressure ESTIMATE `p_est = p_old + dp` with
`dp = -mean(C:inc_ept)` evaluated with a TRIAL C built with the
first-loading modulus (same pattern as
`group_materi_elasti_stress_pressure_history_factor`, lot 6): this
catches the FIRST unloading step. If `|p_est| < old_unknowns[sph_indx]`
the material is unloading/reloading -> Eur/nuur, otherwise first loading
-> E50/nu50. `matrix_a4b` is NOT in-place safe (scratch `work`).
IEEE `-0.0` normalized (`p == 0. -> p = 0.`).

**Clamps (documented decision)**: (1) base `<= 0` (very tensile sig3
beyond the cohesion) -> `E = Eref` (the power law is not evaluated on a
non-positive base; avoids zero/negative stiffness); (2) `m == 0` ->
`E = Eref` (avoids `0^0` and makes the law exactly linear);
(3) `Eref_50 <= 0 || Eref_ur <= 0` or a non-positive reference base ->
db_error.

**Coupling with plasticity (documented decision)**: the elastic C is
evaluated at the START-of-step stress (one step lag, like all
pressure-dependent laws of this code); the plastic law evaluates its own
E50/Eur at the CURRENT trial stress (see the plasti page). The manual
couples them; the one-step lag is the tochnog convention.

## Gotchas

- `new_sig` at the top of `set_stress` is the start-of-step stress
  (copied from the VERSION_NORMAL dofs), NOT the current iteration
  stress.
- With `m > 0`, `c = 0` and `sig3 ~ 0` the power law is EXTREMELY
  sensitive (a numerical `sig3` residual of 1e-2 collapses `E` by two
  orders of magnitude: `E = Eref*(1e-4)^0.5 ~ 1% Eref`): use cohesion
  and/or confinement, or `m = 0` (documented; the tests use `c = 10`
  which keeps the base at `c*cot(phi) = 17.32`).
