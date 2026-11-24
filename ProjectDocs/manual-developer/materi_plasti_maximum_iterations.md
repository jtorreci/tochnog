# materi_plasti_maximum_iterations

## Files and functions

- `database.cc`: two records, registered next to each other.
  - `GROUP_MATERI_PLASTI_MAXIMUM_ITERATIONS` (GNU, per group): INTEGER,
    length 1, class MATERI, requires `GROUP_TYPE`.
  - `MATERI_PLASTI_MAXIMUM_ITERATIONS` (F5 T3 B, manual Professional
    6.804): INTEGER, length 1, class MATERI, `no_index`.
- The enum `MATERI_PLASTI_MAXIMUM_ITERATIONS` is appended in `tochnog.h` /
  `tochnog-mod.h` (F5 T3 B block).
- `stress.cc`, routine `stress()`: the block that sets `max_plasti_iter`.

## Implementation details

- Defaults: `materi_plasti_f_nonlocal && viscoplasti` → 2,
  `viscoplasti` → 10, otherwise `MAX_PLASTI_ITER` (20000, manual 6.804,
  owner 2026-11-06, BACKLOG C112; `MAX_ITER` 1000 before). `MAX_ITER`
  (1000) still bounds the membrane plane-stress loop.
- Override: `db( GROUP_MATERI_PLASTI_MAXIMUM_ITERATIONS, gr, ...,
  GET_IF_EXISTS ) || db( MATERI_PLASTI_MAXIMUM_ITERATIONS, 0, ...,
  GET_IF_EXISTS )`, so the group record wins over the global one. A value
  below 1 prints an error naming both records and exits.
- Loop guard: `while ( !plasti_found && total_plasti_iter<max_plasti_iter )`.
  Reaching the cap is silent.

## Verification

Professional 02-08-2026 probes on the `plasti_apply1` square without
`plasti_apply -no`: default sigyy -3.39353 (GNU equal), cap 1 -21.9773
(GNU -21.978), cap 0 accepted as 1. Suite: `materi_plasti_maximum_iterations1`,
`_2` (precedence), `_error1`, `_default1` (default 20000). Corpus impact
of the default: BACKLOG C112.

## Pending

- The viscoplastic defaults 2 and 10 are hardcoded numbers.
- The error message does not carry the group number.
