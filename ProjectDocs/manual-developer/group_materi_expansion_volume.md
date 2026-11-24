# group_materi_expansion_volume — developer notes

## Where

- `materi.cc` `materi()`: `dens = (1 - beta*T) * dens` before the gravity
  force (`h * dens * force_gravity`), with `T = new_unknowns[temp_indx]`.
- `general.cc` `general()`, `type==-MATERI` branch (K04 incnav3,
  2026-10-06): the same reduction right after `get_materi_density()`, so
  the lumped inertia of `MATERI_VELOCITY` (`inertia = dens`) and the
  convective part (`conv_part = dens`) use rho = rho0*(1 - beta*T).

## Evidence

Pro test `incnav3` (quad4, T = 0.5, beta = 1, rho0 = 1, gravity +1,
compressibility 1, one step of 0.01): `tvely(2)` GNU before -2.4976, after
-4.990415335 = Professional 02-08-2026 (-4.990415335463), also on the
converted sfnet `incnav3`.

## Pending / hardcoded

- `group_materi_expansion_linear` is read next to it in `materi.cc`; no
  change.
- The reduction is evaluated per integration point from the interpolated
  temperature; no clamp for beta*T >= 1 (negative density).
