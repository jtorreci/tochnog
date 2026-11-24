# group_materi_elasti_borja_tamagnini (developer)

## Where

- `database.cc`: `GROUP_MATERI_ELASTI_BORJA_TAMAGNINI` (DOUBLE, 4, MATERI);
  enums in `tochnog.h` / `tochnog-mod.h`; `check.cc` requires
  `materi_stress` and `materi_strain_elasti`.
- `hyperela.cc` `hyperelasticity()`: analytic branch before the
  energy-difference laws. Stress from the elastic strain `epe` (eq. 3.12,
  `P0 = -pr`, `eps_v0 = 0`) and the exact tangent
  `C = (k_t - 2mu/3) I1 + 2 mu I - 2 alpha pr exp(Omega)/khat (1 e + e 1)`,
  `k_t = pr exp(Omega)/khat (1 + 3 alpha eps_s^2/(2 khat))`, into `Chyper`.
  `set_stress()` calls it for the old and the new elastic strain and adds
  the difference, as for the other hyperelastic laws; the tangent enters
  `ddsdde` through `Chyper`.

## Implementation details

- Small-strain invariants (tochnog tensor strains); the paper uses
  logarithmic principal stretches, which coincide for small strains.
- khat <= 0 is an error.
- Not combined with `group_materi_hyper_stiffness` (the tangent is always
  analytic).

## Verification

Analytic (homogeneous hex8): isotropic compression -190.8553692 for 4 steps
and 1 step; with initial stress -pr, compression and shear: -202.8639229 /
5.0171074; pure shear 1.2 / -0.1. The Professional 02-08-2026 differs (see
the user page and DIVERGENCES BORJA-TAMAGNINI).

## Pending

- Membrane stress state not checked.
