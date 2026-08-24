# force_edge / force_volume aliases + control_materi gates (developer)

## Prefix translation — the critical placement

First attempt translated the keyword in input.cc right before
`idat = db_number(str)`. BROKEN: the end-of-values detection of
variable-length records also calls db_number (input.cc stops reading
values when db_number(str)>=0). With the Professional name that check
returned -1, the parser tried to read `force_volume` as a double and
died with "Problem reading : bounda_time". FIX: the translation lives
INSIDE db_number (database.cc) so every caller sees it. Rule: keyword
aliases must resolve at the db_number level, never at a single call
site.

## Files

- tochnog.h / tochnog-mod.h: 24 new enums (edge/water/volume variants +
  13 CONTROL_MATERI_*), 1080 in sync.
- database.cc: 24 registrations; db_number prefix translation.
- area.cc: force_edge_companion helper (element/group/side/node/
  element_node per family; -1 = not available), element-level
  restriction block + node-level checks and node_factor application in
  the three force branches; water_factor in the water branch.
- force.cc: _element/_element_group restriction in
  force_element_volume.
- general.cc: control_materi_gate_off(item) — reads ICONTROL-indexed
  switch, 1 when -no.
- viscosit.cc / damage.cc / stress.cc: the five wired gates (viscosity,
  damage+failure, updated [-no direction only], plasti_tension,
  plasti_visco). The updated_apply gate sits at the single
  GROUP_MATERI_MEMORY read of set_stress (canonical point).

## Verification (suite 91/91)

fedge_alias (pure Professional syntax, sigxx 5.0), fedge_restrict
(_element off-geometry -> disx 0; needs derivatives+displacement+
total_linear), fvol_elem (group restriction: 0.5/0.0), cmat_gate
(tension cutoff ignored -> linear -88.4 vs capped ~1; A/B fails without
the gate).
