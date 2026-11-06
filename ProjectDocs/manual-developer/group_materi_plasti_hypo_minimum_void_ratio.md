# group_materi_plasti_hypo_minimum_void_ratio

## Files and functions

- `database.cc`: `GROUP_MATERI_PLASTI_HYPO_MINIMUM_VOID_RATIO` (F5 T3 C1,
  DOUBLE, length 1, class MATERI, requires `GROUP_TYPE`). Enum appended in
  `tochnog.h` / `tochnog-mod.h`.
- `hypoplas.cc`, `hypoplasticity()`: reads the record for the group into
  the thread-local `hypo_minimum_void_ratio` and sets
  `hypo_minimum_void_ratio_on`.
- `hypo.c`, `sigma_()`: after `e = eold + edot*dt`, `e` is raised to the
  minimum (von Wolffersdorff and lowangles, every substep, before the
  stress uses it).
- `hypoplas.cc`, Masin block: `new_hisv[e_slot]` is raised to the minimum
  after the kernel returns (the kernel itself is not changed).

## Verification

Professional 02-08-2026 probes: `other/hypo1` + 0.599 sigyy -848.5617,
hyhis0 0.606516 (GNU -848.5630, 0.606528); `other/hypo7` + 1.1207 sigyy
-224.8340, hyhis0 1.124789 (GNU -224.8345, 1.124790).

## Pending

- The SANISAND route ignores the minimum.
- In the Masin route the kernel integrates the whole step with the
  unclamped void ratio; only the stored value is reset.
