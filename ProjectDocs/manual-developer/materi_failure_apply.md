# materi_failure_apply — developer

## Implementation

- **Enum**: `MATERI_FAILURE_APPLY` in `tochnog.h` / `tochnog-mod.h` (in
  sync), right after `CONTROL_MATERI_FAILURE_APPLY`.
- **Record** (`database.cc`): INTEGER, length 1, `no_index 1`, data_class
  CONTROL. `control_materi_failure_apply` stays indexed (manual 6.143).
- **Consumer** (`failure.cc`, `failure()`): `failure_off =
  materi_apply_gate_off( MATERI_FAILURE_APPLY, CONTROL_MATERI_FAILURE_APPLY
  )` (helper in `general.cc`, see `materi_damage_apply.md`); the detection
  `if ( !failure_off && tmp>threshold && !ELEMENT_DELETE_TIMES )` is
  skipped, so no new `ELEMENT_DELETE_TIMES` record is created. Elements
  already failing keep being deleted according to their stored schedule
  (design choice: "failure data ignored" = the criteria are not evaluated;
  could not be measured, see below).
- **Removed**: the old `CONTROL_MATERI_FAILURE_APPLY` test in `damage()`.
  Before this change `control_materi_failure_apply -no` froze damage and
  did not touch the failure criteria, the opposite of the Professional.

## Verification (Professional binary 02-08-2026 as oracle)

No Professional distribution test uses the failure switches. Own model:
two stacked quad4 (element 1 group 0 with `group_materi_failure_rupture 0
0.005 0.2`, element 2 group 1 without), pulled in x (velx 0.02 constant,
t = 0.5), `sigxx` at nodes 1 and 5.

| Case | Pro | GNU |
|------|-----|-----|
| (a) failure active | stops: "Program error detected for data item: element_rhside_delete" (any delete_time) | 8.07692 |
| (b) `materi_failure_apply -no` | 13.4615 | 13.4615 |
| (c) failure record removed | 13.4615 | 13.4615 |
| `control_materi_failure_apply -no` (Pro: unindexed) | 13.4615 | 13.4615 (GNU: index 0) |
| `materi_damage_apply -no` | program error (failure still active) | 8.07692 |
| global `-no` + control `-yes` | program error (failure active) | 8.07692 |
| global `-yes` + control `-no` | 13.4615 | 13.4615 |

Findings: `materi_failure_apply -no` equals removing the failure data;
the control record overrules the global one; damage and failure switches
are independent. The Professional binary registers
`control_materi_failure_apply` without index (indexed form = parse error),
contrary to its manual; the GNU keeps the manual's indexed form.

Self-checking inputs `validation-suite/test-2026/materi_failure_apply1.dat`
(passes on both binaries) and `materi_failure_apply2.dat` (GNU only).

## Pending

- The Professional's element deletion crashed in every small model tried
  (rupture, several delete_time values), so case (a) and the behaviour of
  elements that failed before a later `-no` timestep are unverified.
