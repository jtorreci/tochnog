# bounda_constant

## Files and functions

- `bounda.cc` — within `bounda()`: state variable `bounda_constant` (line 45),
  read of the keyword (lines 152–154, `GET_IF_EXISTS`), application in the
  time branch of the node loop (lines 444–447).
- `tochnog.h` — enum `BOUNDA_CONSTANT` (line 140).
- `tochnog-mod.h` — mirror enum `BOUNDA_CONSTANT` (line 133), must stay in sync.
- `database.cc` — keyword registration (lines 184–188).

## Implementation details

- Read: `db( BOUNDA_CONSTANT, iboun, &bounda_constant, ... )` returns an
  integer; active when equal to `-YES`.
- Application (time branch, `else` of the rotation/sine block): if
  `bounda_constant==-YES`, `iuknwn>=vel_indx`, `!force` and
  `node_dof[iuknwn]!=0.`, the new value is copied from the previous step
  (`VERSION_NORMAL`): `new_node_dof[iuknwn] = node_dof[iuknwn]`. Otherwise the
  normal `factor * load` path is used. The zero-check means a dof that has never
  been active is not locked at zero.
- The previous-step value is used, so the dof is "frozen" at whatever value the
  last computed step had (first step still applies `bounda_time`).
- `database.cc`: `type = INTEGER`, `data_length = 1`, `data_class = BOUNDA`,
  `data_required = BOUNDA_UNKNOWN`.

## External dependencies

None. Internal database API only.

## Hardcoded parameters / pending refactorings

- Activation flag is a plain integer compared against `-YES`; consistent with
  the rest of `bounda.cc` but not enforced as a boolean type.
- The `node_dof[iuknwn]!=0.` condition silently ignores nodes whose previous
  value is exactly zero (they keep following `bounda_time`); intended but
  undocumented, review if zero-valued constant dofs are needed.
- Only applies to the `unknown` branch (velocities/pressures/temperatures);
  applied forces (`bounda_force`) are excluded via `!force`.
