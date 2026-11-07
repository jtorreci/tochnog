# group_materi_plasti_hypo_void_ratio_linear (developer)

## Where

- `database.cc`: `GROUP_MATERI_PLASTI_HYPO_VOID_RATIO_LINEAR` (INTEGER, 1,
  MATERI); enums in `tochnog.h` / `tochnog-mod.h`. The former `db_number`
  alias to `GROUP_MATERI_PLASTI_HYPO_PRESSUREDEPENDENTVOIDRATIO` is removed
  (BACKLOG C118).
- `hypoplas.cc` `hypoplasticity()`: reads the record (default `-YES`) into
  the thread-local `hypo_void_ratio_exact` (defined in `hypo.c`, the same
  pattern as `hypo_minimum_void_ratio_on`).
- `hypo.c` `sigma_()`: after the linear update, when `hypo_void_ratio_exact`
  the void ratio is `(1 + e_old) det(I + dstr) - 1`, `dstr` the strain
  increment of the substep; the minimum void ratio clamp follows.

## Verification

Professional `other/hypo11` (biaxial, Karlsruhe sand) without
`_wolffersdorff_niemunis`: `-yes` sigyy -1.0198992, e 0.6588397
(Professional -1.019954, 0.6588469), `-no` -1.0199336, 0.6588243
(-1.019924, 0.6588206). With the old alias the GNU gave -1.15668.

## Pending

- `control_materi_plasti_hypo_pressure_dependent_void_ratio` (6.147), the
  Professional switch of the pressure dependent correction, is registered
  but not read; the GNU switch is the group record
  `group_materi_plasti_hypo_pressuredependentvoidratio`.
- The Masin kernels do not read the record.
