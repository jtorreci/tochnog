# group_materi_elasti_c / group_materi_elasti_c_direction (developer)

## Where

- `database.cc`: `GROUP_MATERI_ELASTI_C` (DOUBLE, length MDIM^4 = 81) and
  `GROUP_MATERI_ELASTI_C_DIRECTION` (DOUBLE, length MDIM^2 = 9), MATERI,
  `GROUP_TYPE`; enums appended in `tochnog.h` / `tochnog-mod.h`.
- `check.cc`: both need `materi_stress` and `materi_velocity`.
- `stress.cc` `set_stress()`: block after the last elastic law (TSKH /
  smallstrain), because the young_power, hardsoil and other branches clear
  C / Cmem before filling them. The record value index ((i*3+j)*3+k)*3+l
  is C'_ijkl; with the normalized axes dir_p (rows of c_dir)
  C_ijkl += dir_pi dir_qj dir_rk dir_sl C'_pqrs, added to both C and Cmem.

## Implementation details

- Membrane (plane stress): refused with
  `Error: group_materi_elasti_c is not available with group_materi_membrane -yes.`
  (no plane-stress condensation; the Professional also refuses it).
- Zero axis: `Error: group_materi_elasti_c_direction has a zero direction.`
- The axes are not checked for orthogonality.

## Verification

Professional 02-08-2026 probes on the unit square (plane strain,
total_linear, oedometer-like and shear): isotropic C equal to young 1000 /
poisson 0.3 (sigyy -10.98901099, sigxy 3.846153846), added to young /
poisson (-21.97802198), orthotropic C_xxxx 3000 rotated so local x = global y
(sigyy -27.52747253, the hand value lambda eps_xx + 3000 eps_yy), axes of
length 2 (same). GNU equal to 1e-9.

## Pending

- Plane-stress condensation for membranes.
- Element-specific values follow get_group_data (as every group record).
