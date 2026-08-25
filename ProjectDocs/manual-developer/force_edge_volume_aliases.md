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

## force_edge_projected (Sprint 9 lote 2)

Master as area.cc type 9 (MTYPES 10) with the full companion set in
force_edge_companion (incl. _node_factor). The projection: ph/pv from
the linear field at the NODE coordinates; vd normalized (fallback
(0,-1,0)); hd = tunnel x vd (2D: (-vd_y, vd_x, 0)); ndothd/ndotvd/
tdothd/tdotvd explicit inner products (the first draft accumulated
n2/t2 inside a loop with a dimensionally WRONG scalar — always compute
t·σ·n as ph(t.hd)(n.hd) + pv(t.vd)(n.vd)). Tangent from the outward
normal (2D rotation). The load sign follows the force_edge_normal
convention (positive pushes along +n, i.e. outward into the void).

Test gotcha (bit twice now): quad4 Z-convention — node 3 is TOP-LEFT,
node 4 TOP-RIGHT. A crossed mesh gives sigxx=0 mysteries and
"-ra 2 3" velx constrains bottom-right + top-LEFT (crossed). Verified
with force_edge_normal on the same mesh first (known-good family) to
separate test bugs from feature bugs.

Verification: fproj_tunnel (sigxx=10 ph exact on the vertical wall,
sigyy=20 pv on the horizontal wall; 2:1 ratio discriminates the
projection from a uniform pressure). Suite 92/92.

## Lote 3 (data family + solver aliases)

- data_activate/data_delete: time-gated loops in data() BEFORE the
  control_data_* blocks (same destructive patterns). EPS_TIME gate.
- data_ignore: input.cc, right after idat resolution — the record loop
  condition `while(strcmp(str,"end_data"))` re-tests str, so consuming
  until db_number>=0 and `continue` handles both the next keyword and
  end_data. GOTCHAS: the ignored item enum is stored NEGATIVE
  (match -idat), and the data_ignore record must precede the records it
  ignores (it is evaluated while reading).
- db_number explicit aliases: control_solver→control_options_solver,
  control_solver_bicg_error→control_options_solver_bicg_error,
  axisymmetric→group_axisymmetric, bounda_print_mesh_dof*→print_mesh_dof*.
- control_solver_bicg_stop -no: wired in so_bicg.cc (warn+continue
  instead of exit; GET_IF_EXISTS only — safe in the solver context).
- print_mesh_dof: one-shot dump in data() (static flag), file
  print_mesh_dof.dat; _geometry filter via geometry(); _values
  registered but unused by the dump (documented).
- DATA has no data_class — the data_* records use CONTROL.
