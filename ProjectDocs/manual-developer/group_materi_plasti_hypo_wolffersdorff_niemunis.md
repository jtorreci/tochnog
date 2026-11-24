# group_materi_plasti_hypo_wolffersdorff_niemunis (developer)

## Where

- `database.cc`: `GROUP_MATERI_PLASTI_HYPO_WOLFFERSDORFF_NIEMUNIS` (DOUBLE,
  2, MATERI); enums in `tochnog.h` / `tochnog-mod.h`; `check.cc` requires
  `materi_stress` and the hypo history, as the base record.
- `hypoplas.cc` `hypoplasticity()`: reads nu z into the thread-local
  `hypo_niemunis_on`, `hypo_niemunis_nu`, `hypo_niemunis_z` (defined in
  `hypo.c`) and checks -1 < nu < 0.5.
- `hypo.c` `sigma_()`, only for `ihypotype == 0` (von Wolffersdorff):
  1. after `fb`: `f_d_bar` from eq. 4.221 with `M_e = 1`,
     `M_T = -n e_d/h_s (-tr T/h_s)^{n-1}` (eq. 4.211), and the signed
     interpolation 4.223 (`f_d < 1` only); it replaces the reference
     `f_d` and its `e < 1.0001 e_d` clamp.
  2. after the second linear term: `lmat += fac b2 (delta_ik delta_jl -
     1/3 delta_ij delta_kl)`, `b2` from eq. 4.183.
  3. after `nmat`: `N += b2 dev(Lhat^-1 : N)` with
     `Lhat^-1 = (I - a^2 That That / (F^2 + a^2 That:That)) / F^2`
     (Sherman-Morrison). This keeps the flow rule m of eq. 4.138
     unchanged, as the thesis requires below eq. 4.183; adding b2 to L
     alone gave sigyy -64.5 on `other/hypo11` (Professional -0.977).

## Verification

Professional 02-08-2026 on `other/hypo11` and variants (sigyy at t 0.1):
nu 0.2 -0.9769159 / -0.9768868, nu 0.365 (b2 ~ 0) -1.0172854 / -1.0173444,
nu 0.3822 -1.0455057 / -1.0454883, nu 0.49 -0.3668340 / -0.3667830.
Dense e0 0.48 (f_d_bar active): -7.413458 / -7.418368 at t 0.03; both
collapse near t 0.04 and end 0.4 % apart. e0 0.50 with z 1 follows the
Professional to 0.1 % up to the collapse at t 0.09. Variants tried and
rejected against these histories: e_d0 instead of e_d in M_T, p instead
of -tr T.

## Pending

- The thesis derives `f_d_bar` for the isotropic direction D_bar (eq.
  4.221, as the manual); the general D_bar of eq. 4.216 is not used, as
  in the Professional.
