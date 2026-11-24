# bounda_found

## Files and functions

- `database.cc` — keyword registration (lines 178–181): name `bounda_found`,
  `type = INTEGER`, `data_length = 1`, `data_class = BOUNDA`.
- `bounda.cc` — within `bounda()`: the internal `found` logic already exists.
  The local variable `found` is reset per boundary (`found = 0`, line 200) and
  set to 1 in the node loop (lines 287–320) whenever a node matches the boundary
  selection (`-range`, `-all`, geometry, `-node_set`, or a direct node); the
  condition `if ( found )` (line 321) then applies the bounda to that node.
- `tochnog.h` — enum `BOUNDA_FOUND` (line 140).
- `tochnog-mod.h` — mirror enum `BOUNDA_FOUND` (line 133), must stay in sync.

## Implementation details

- The flag can be read/written with `db( BOUNDA_FOUND, iboun, ... )`, with one
  integer value per boundary: `-YES` when the boundary was applied to at least
  one node, `-NO` otherwise.
- The actual "found" detection is already in place: a boundary is considered
  found when `found` becomes 1 inside the node loop, i.e. at least one active
  node matched the boundary selection and the bounda was applied to it.
- The keyword is registered with `data_required` unset (unlike
  `BOUNDA_CONSTANT`, which requires `BOUNDA_UNKNOWN`), so it works alongside
  both `bounda_unknown` and `bounda_force` records.

## External dependencies

None. Internal database API only.

## Hardcoded parameters / pending refactorings

- PENDING: the flag is never written to the database. `bounda_found` is
  registered (record exists and the input keyword is accepted) but no code
  stores `-YES`/`-NO` when a boundary is used or unused, and nothing prints the
  result. The write should happen after the node loop in `bounda()`, using the
  accumulated `found` per `iboun`.
- The `found` variable is reused for two different purposes (time-interpolation
  match in the time branch and node-match in the node loop); do not confuse the
  two scopes when adding the write.
- Report/print behaviour of the `-yes` flag is not yet defined (no output
  mechanism wired to the record).
