# processors_maximum — developer

## Implementation

- Record `PROCESSORS_MAXIMUM` (INTEGER, data_length 1, `no_index`;
  `database.cc` after `OPTIONS_PROCESSORS`, enum after `OPTIONS_PROCESSORS`
  in `tochnog.h` / `tochnog-mod.h`).
- `top()` (`top.cc`), where `OPTIONS_PROCESSORS` is resolved: when no
  `OPTIONS_PROCESSORS` record exists (the `processors` alias and the
  `TOCHNOG_OPTIONS_PROCESSORS` environment variable, stored as the record
  by `miscel.cc`, both create it), `processors_maximum -yes` sets
  `number_processors = available_processors()`. The existing clamps
  follow: `MTHREAD` (64) and one thread with `group_user_umat`; the value
  is stored back in `OPTIONS_PROCESSORS`.
- `available_processors()` (`top.cc`): `CPU_COUNT` of
  `sched_getaffinity()`, fallback `sysconf(_SC_NPROCESSORS_ONLN)`, at
  least 1.

## Hardcoded parameters

- GNU default `-no` (Professional `-yes`), comment in `top()`.

## Threading defects found

- `force_point_calculate()` ran inside `parallel_new_dof_before()`, which
  `parallel_sys_routine()` runs once per thread: the point load was added
  `nthread` times (disx 8 instead of 2 on 4 threads) and raced on
  `NODE_RHSIDE`. It now runs serially in `top()` before both
  `parallel_sys_routine( &parallel_new_dof_before )` calls. Test
  `force_point_threads1` (`options_processors 4`).
- Known single-thread-only path: the `node_support_edge_normal_force`
  accumulation (`area.cc`, warning).

## Pending (F6, BACKLOG C92)

- `parallel_sys_routine()` creates and joins `nthread` pthreads per call
  and `element_loop()` allocates `nthread` copies of the nodal arrays per
  sweep: high system time on small models. A persistent thread pool and
  per-thread buffers allocated once would allow the Professional default.
- `processors_used` output record not written.
