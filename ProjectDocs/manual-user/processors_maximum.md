# processors_maximum

## Description

`processors_maximum switch` (manual Professional 6.1005): with `-yes` the
number of threads (the `processors` record) is set to every CPU the run
may use. It is not used when `processors` (or `options_processors`, or
the `TOCHNOG_OPTIONS_PROCESSORS` environment variable) is given.

## Usage

```
processors_maximum switch
```

- `switch`: `-yes` or `-no`. GNU default `-no` (one thread), see below.

The number of CPUs is the affinity mask of the process (so `taskset` and
container CPU limits are honoured), capped at 64 threads. With
`group_user_umat` the run always uses one thread. The threads used are
written to the `.dbs` as `options_processors`.

## Example

```
processors_maximum -yes
```

Tests: `validation-suite/test-2026/processors_maximum1` (absent: 1
thread), `processors_maximum2` (`-yes`: every CPU), `processors_maximum3`
(`processors 2` wins over `-yes`).

## When to use it

Threads pay off on large meshes only: a 22 500-element plastic quad4 model
runs 1.7x faster on 4 threads, and no faster on 16 (the solver is serial).
On small models they slow the run down (80 elements, 300 steps: 1.3x
slower on 4 threads, 5.5x on 16), because the GNU starts and joins the
threads in every parallel region. Use `processors 4` for large models.

## Differences with the Professional version

- Default: the Professional default is `-yes` (it writes `processors_used
  16` on a 16-core machine); the GNU default stays `-no` until the thread
  overhead is removed in F6 (project owner decision 2026-11-02;
  DIVERGENCES `PROC-MAX`, class B; BACKLOG C92).
- The Professional writes the threads used as `processors_used`; the GNU
  writes `options_processors`.
