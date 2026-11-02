# groundflow_apply

## Description

`groundflow_apply switch` (manual Professional 6.555, data part, no
index): with `-no` the groundflow equation is skipped and all groundflow
data is ignored, for all timesteps. Default `-yes`.

Useful to run a coupled soil model drained (no pore pressure) without
editing the groundflow records out of the input.

Measured on the Professional 02-08-2026 and done the same way in the GNU:

- every groundflow dof (`pres`, `gvel*`, saturation, pressure gradient)
  is 0 at every node, also with an initial `node_dof` value;
- the groundflow dofs are removed from `bounda_dof` (a prescribed pressure
  is ignored) and `control_reset_dof` (a pressure reset is ignored);
- the solid sees no pore pressure: the consolidation column of
  `k24_consol_3d1` gives the drained elastic answer.

## Usage

```
groundflow_apply -no
```

Tests: `validation-suite/test-2026/groundflow_apply1` (one quad4,
prescribed pressure 5 and initial pressure 7 ignored: pres 0 everywhere),
`groundflow_apply2` (3D consolidation column, disx 0.1, sigxx -1, pres 0
= Professional).

## Differences with the Professional version

None measured. The GNU also deletes every `groundflow_*`,
`group_groundflow_*`, `group_interface_groundflow_*` and
`control_groundflow_*` record ("all groundflow data is ignored"), so a
phreatic level adds no static pore pressure either; the Professional
probes had no phreatic level, so that part follows the manual text only.
