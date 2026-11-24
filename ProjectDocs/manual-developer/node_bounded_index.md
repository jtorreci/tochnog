# node_bounded_index (developer)

## Where

- `bounda.cc`, `bounda()`: `db_set_int( NODE_BOUNDED_INDEX )` next to
  `NODE_BOUNDED`, every slot set to `-NOTHING`; where a `bounda_dof` /
  `bounda_unknown` record bounds a dof (`node_bounded[ipuknwn] = 1`) the
  slot gets `iboun`. Rebuilt at every `bounda()` call, so the record holds
  the state of the last step.
- `database.cc`: INTEGER, length `npuknwn`, `print_only`, class NODE.
- Enum appended before `LAST_DUMMY` in `tochnog.h` / `tochnog-mod.h`.

## Notes

- Dofs bounded by other mechanisms (groundflow phreatic / seepage, mpc,
  bounda force-type branches) set `node_bounded` but have no `bounda_dof`
  index: they stay `-nothing`, as the manual describes ("the index of the
  bounda_dof record").
- Test `node_bounded_index1`.
