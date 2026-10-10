# db() buffer contract and the scratch-buffer audit (F4 T13, BACKLOG C49)

Internal contract, no user-visible keyword. It explains how to size the
buffers given to the database getter and how to rerun the mechanical audit
(`scripts/audit_db_buffers.py`) that checks every call site.

## The contract

`db( idat, index, ival, dval, length, version, task )` (`database.cc`):

| task | what it does with the buffers | `length` |
|------|-------------------------------|----------|
| `GET` | writes the **stored** length of the record into `ival` (INTEGER records) or `dval` (DOUBLE records) | output only: the value passed in is ignored and overwritten |
| `GET_IF_EXISTS` | as `GET` when the index is active; returns 0 and writes nothing otherwise | output only |
| `GET_AND_CHECK` | writes `length` values, `db_error` when the stored length differs | input |
| `PUT` | reads `length` values (must be `<= data_length`) | input |

The stored length of a record is at most `data_length[idat]`
(`db_data_length()`); for a fixed-length record (`fixed_length = 1`, the
default) it is always exactly `data_length`. A variable-length record holds
what the input file or the code put there, so a short buffer overflows as
soon as an input gives more values than the code uses.

Rules:

1. The buffer written by a read (`ival` for INTEGER, `dval` for DOUBLE)
   must hold `data_length[idat]` values. The buffer of the other type is
   never touched and may be the usual `idum[1]` / `ddum[1]` dummy.
2. A scalar (`&x`) is only valid for a record of `data_length` 1.
3. Capping `length` before a `GET` does **not** limit what is written
   (`interface.cc` capped `length_stiff` to 3 and still received up to
   `DATA_ITEM_SIZE` values into `ddum3[3]`).
4. When only the first values are needed, read into a buffer of the record
   length (`DATA_ITEM_SIZE` for the generic variable records) and copy the
   values used; the stack cost is irrelevant (190 doubles = 1.5 kB).
5. For a record chosen at run time (`db( data_item_name, ...)`), dispatch on
   `db_type()` and give the matching buffer `DATA_ITEM_SIZE` values.
6. The same applies to `get_group_data( idat, gr, element, unknowns,
   values, n, task )` (`group.cc`): `values` receives up to
   `data_length[idat]` doubles.
7. `array_*( ..., n )` helpers (`array_distance( a, b, work, n )`,
   `array_move`, `array_set`, ...) touch `n` values of every array argument,
   including the work array of `array_distance`.
8. Other helpers with output arrays have their own size contract; the
   audit knows `geometry( inod, co, entity, in, factor, normal, penetration,
   projection, ... )` (`normal` is always cleared over `MDIM` values,
   `projection` gets `ndim` values, `co` is read when `inod < 0`: all three
   need `MDIM`) and the static `msf_element_rule()` of `calcul_force.cc`
   (`nper[3]`). Add new ones to `CALLEE_CONTRACT` in the script.

Why it matters: the overflows are silent at `-O1` (the stack layout happens
to put a dead variable after the buffer) and break at `-O2`
(`print_beam_force_moment` cut distance 0, `creset_iface` / Pro
`interface7` / `interface10` / `mesh_cut_2` database errors, F4 T11).

## The audit script

```bash
scripts/audit_db_buffers.py            # hits + summary; exit 1 when a hit exists
scripts/audit_db_buffers.py --all      # also REVIEW and unresolved calls
scripts/audit_db_buffers.py --tsv out.tsv   # every scanned call with its verdict
```

What it does (pure text analysis, no compiler):

1. Constants: `#define` of `tochnog.h` plus the `#define` of each source
   file. Run-time lengths are evaluated at their upper bound: `ndim` =
   `MDIM`, `nuknwn` = `MUKNWN`, `npuknwn` = `MPUKNWN`, `nprinc` = `MPRINC`,
   `npointmax` = `MPOINT`, `materi_maxwell_stress` = `MMAXWELL-1`
   (`input.cc` asserts `< MMAXWELL`), `nder` = 3.
2. Record table: every `type[X]`, `data_length[X]` (last assignment wins,
   as at run time) and `fixed_length[X]` of `database.cc`.
3. Calls: every `db(...)` with a reading task, every `get_group_data(...)`,
   every `array_*(...)` with a constant length and every call of a
   function listed in `CALLEE_CONTRACT` (`geometry`, `msf_element_rule`),
   in every `.cc` / `.c`.
4. Buffer size: nearest declaration before the call in the same function,
   then the function parameters, then file scope (globals and struct
   members), then the `tochnog.h` globals. `name[N]` gives N (products for
   multi-dimensional arrays), `&scalar` / `&s.member` 1, `&a[k]` / `a+k`
   N-k.
5. Verdicts: `HIT` (buffer < length: a real or latent overflow), `ok`,
   `REVIEW` (run-time record with a buffer below `DATA_ITEM_SIZE`: usually
   the dummy of the other type behind a `db_type()` check),
   `unresolved:*` (pointer, parameter, variable offset), `skip:*`
   (write tasks; run-time record with large buffers).

Limits: a pointer buffer (`get_new_dbl()`, `db_dbl()`) or an array
parameter is not sized (the caller owns the size); a variable offset
(`&coord[inol*ndim]`) is not evaluated; records chosen at run time are only
flagged for review.

## Result of the F4 T13 audit (2026-10-10)

Scanned 3450 calls in HEAD `dca225a`: 1972 `db()`, 129
`get_group_data()`, 1293 `array_*()` with a constant length, 51
`geometry()`, 5 `msf_element_rule()`.

- 45 hits at 34 call sites, all fixed; after the fix the script reports 0
  hits (exit 0). None was a false positive once the script resolved struct
  members, prototype parameters and the Maxwell bound (the first versions
  flagged `plasti_i.cc` `sig[MDIM*MDIM]` against a prototype `sig[6]` of
  `tochnog.h` and `viscoela.cc` `maxwell_chain` against an unbounded
  `materi_maxwell_stress`).

| site | record / call | buffer -> fix | class |
|------|---------------|---------------|-------|
| `bounda.cc` `bounda()` | `bounda_factor`, `bounda_factor_parabolic_x` (variable, `DATA_ITEM_SIZE`) | `[4]` / `[3]` -> `[DATA_ITEM_SIZE]` | latent (input with more values) |
| `bounda.cc` `bounda()` | `groundflow_phreaticlevel` (variable) | `&wl` -> `phreatic[DATA_ITEM_SIZE]`, `wl = [0]` | real when the level has more than one value |
| `calcul_force.cc` `msf_calculate_node_2d/3d()` | `post_calcul_materi_stress_force_plot_switch` | `[3]` / `[4]` -> `[DATA_ITEM_SIZE]` | latent |
| `calcul_force.cc` `msf_element_contribution_3d()` | `..._direction_exclude/include` | `dir[MDIM]` -> `[DATA_ITEM_SIZE]` | latent |
| `data.cc` `data()`, `interface.cc` `interface_element()`, `print_interface_stress.cc` | `group_interface_materi_elasti_stiffness` (`DATA_ITEM_SIZE`) | `[3]` -> `[DATA_ITEM_SIZE]` | latent (the `length_stiff` cap had no effect) |
| `delete.cc` `mesh_cut()` | `control_mesh_cut_node_force` (up to `MDIM` integers) | `idum[1]` -> `ival[MDIM]`, one `GET_IF_EXISTS` | real (`mesh_cut_2` at `-O2`) |
| `general.cc` `general()`, `truss.cc` `truss()`, `top.cc` `post_solve_state_pass_gated()` | `control_options_inertia` (multi-valued, one switch per principal dof in the Professional) | `&options_inertia` -> `inertia_switches[DATA_ITEM_SIZE]`, value `[0]` | real with the Professional multi-value form |
| `geometry.cc` `geometry()` | `array_subtract( xyint, centre, ..., ndim )` | `xyint[2]` -> `[MDIM]` zero-initialised (read overflow in 3D) | latent |
| `materi.cc` `strain_settlement_creep()` | `mesh_activate_gravity_time` (fixed, 2) | `ddum[1]` -> `gravity_time[2]` | real |
| `new_mesh.cc` `new_mesh()` | `control_mesh_new_mesh_region` (`6+nuknwn`) | `&scalar` -> `region[6+MUKNWN]`, value `[0]` | latent |
| `print_beam_force_moment.cc` (2 sites) | `array_distance( p1, p2, ddum, 3 )` | work `ddum[1]` -> `work[3]` | real (`bmom*` cut distance 0 at `-O2`) |
| `print_interface_stress.cc` | `control_print_interface_stress_2d_coordinates` (fixed, 4) | `ddum[1]` -> `coordinates[4]` | real (gcc `-Warray-bounds` at `-O2`) |
| `print_interface_stress.cc` | `element_interface_force_norm`, `_strain_normal`, `_force_tang`, `_force_tang2` (per integration point, up to 4) | scalars -> `values[DATA_ITEM_SIZE]`, value `[0]` | real (`creset_iface`, `interface7`, `interface10` at `-O2`; ASan) |
| `top.cc` `top()` | `control_timestep_iterations_automatic` (up to 3) | `[2]` -> `[3]` (the code read `[2]`) | real with the 3-value form |
| `top.cc` `step_start()` | `control_mesh_copy` (variable) | `control_copy[MDIM]` -> `[DATA_ITEM_SIZE]` | latent |
| `condif.cc` `condif()`, `data.cc` `control_reset_interface()`, `groundfl.cc` `groundflow_phreatic_level_multiple_find_element()`, `mpc.cc` `mpc_element_group_generate()`, `print_interface_stress.cc`, `print_node.cc` `print_node()`, `unknown.cc` `unknown_reset()` | `geometry()` with `ddum[1]` as `co` / `normal` / `projection` | `ddum[1]` -> `ddum[MDIM]` | real: `normal` is cleared over 3 values on every call (`creset_iface` ASan, `-O2` database error) |
| `calcul_force.cc` `msf_element_contribution_2d()` | `msf_element_rule()` writes `nper[0..2]` | `nper[2]` -> `nper[3]` | real (ASan in every `post_calcul -materi_stress -force` 2D test: `msf_*`, `qsri_*`, `gforce7`) |

- 38 REVIEW calls (record chosen at run time): reviewed by hand. They pass
  the matching large buffer after a `db_type()` check or read a record of
  a known type (node coordinates, `-yes/-no` switches, the normal records).
  The sites that assume a DOUBLE record without checking the type
  (`post.cc` `post_integrate()` / `parallel_post_node()`, `print_el.cc`,
  `repeat.cc`, `data.cc` arithmetic and spatial resets, `calcul.cc`)
  would overflow `idum[1]` only if the input names an INTEGER record where a
  DOUBLE one is expected: BACKLOG C54.
- 273 unresolved: 252 pointer buffers, 15 variable offsets (all
  `&coord[inol*ndim]` into `MNOL*MDIM` arrays), 6 array parameters.

## Verification (F4 T13)

- ASan + UBSan build (`-fsanitize=address,undefined -fno-omit-frame-pointer
  -g -O1`, `ASAN_OPTIONS=hard_rss_limit_mb=4000:detect_leaks=0`, no
  `ulimit -v`): build_safe 28 failing runs -> 1; corpus (omen, 6 jobs,
  360 s) 45 tests with a report at 9 sites -> 2. Left: BACKLOG C55
  (`print_unknownsrhside()` heap read, wrong index) and C56 (parser token
  longer than `MCHAR`). No UBSan report.
- `-O1`: build_safe 187 OK / 0 FALLO; corpus `.dbs` identical to HEAD.
- `-O2` (now the default, `scripts/build_safe.sh`): every finished `.dbs`
  identical to `-O1`; the F4 T11 `-O2` failures run.

## Pending refactorings

- Several variable records are declared with `data_length =
  DATA_ITEM_SIZE` although the code uses 1-4 values
  (`control_options_inertia`, the `bounda_factor*` and `plot_switch`
  records); a tighter declared length would turn an over-long input into a
  "Length too small" error instead of silently ignored values.
- `GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS` is assigned twice in
  `database.cc` (3, then `DATA_ITEM_SIZE`).
- A type check for run-time record names (BACKLOG C54).
