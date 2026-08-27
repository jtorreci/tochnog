# control_print_dof_smooth_dof

## Implementación

- Read INSIDE `print_dof()` in `print_hi.cc`: the smoothing is a
  modifier of the `control_print_dof` output, NOT a print itself (it is
  NOT dispatched from `top.cc` and is NOT gated by
  `control_print_frequency_allowed` — the gate only wraps
  `control_print_dof`, and when the gate blocks the print the smoothing
  is irrelevant).
- `print_dof_smooth_apply( icontrol, max_node, nuknwn, smooth_field )`
  (new static helper in `print_hi.cc`): builds the smoothed field and
  returns 1 when active. `print_dof` writes
  `smooth_field[inod*nuknwn_+indx]` instead of `node_dof[indx]` when
  the record is active.
- Registered in `database.cc`: `CONTROL_PRINT_DOF_SMOOTH_DOF` (INTEGER,
  variable length, data_required = CONTROL_PRINT_DOF) and
  `CONTROL_PRINT_DOF_SMOOTH_N` (INTEGER, length 1, data_required =
  CONTROL_PRINT_DOF_SMOOTH_DOF). Enums in `tochnog.h` /
  `tochnog-mod.h` in sync.

## Diseño / decisiones

- **Smoothing formula** (documented decision): each pass replaces the
  value of a node by the AVERAGE OF ITS NEIGHBOUR NODES — nodes
  connected to it by at least one element — WITHOUT including the node
  itself. A node with no neighbours (isolated) keeps its value; an end
  node of a chain averages over its single neighbour. This matches the
  manual wording ("the average of its neighbouring nodes") and makes a
  linear field a fixed point of one pass (interior nodes unchanged),
  while repeated passes converge to the mean value of the field.
- **Adjacency**: built from the ELEMENT records (VERSION_NORMAL; the
  mesh does not change within a step). Each element contributes all its
  node PAIRS (undirected). Stored as head/to/next_edge lists indexed by
  POSITION (the compacted node order of print_dof's renumbering): the
  element node numbers are original, mapped with orig2pos.
- **Position mapping**: `pos2orig`/`orig2pos` built from the active
  NODE indices in ascending order (the same order that
  `renumbering(VERSION_PRINT)` produces). CRITICAL: the node_dof VALUES
  are read by POSITION (`db_dbl(NODE_DOF, pos, VERSION_PRINT)`), NOT by
  the original node number — VERSION_PRINT is indexed by the compacted
  numbers (GOTCHA discovered in test calibration: reading by original
  id shifts the field by one node).
- **Components**: `-all` -> every `nuknwn` component; otherwise the
  listed labels are matched against `dof_label` with `array_member`
  (unmatched labels are skipped silently, same leniency as
  `print_history`). Per-component sums (GOTCHA: an initial
  implementation accumulated ONE shared sum over all components,
  corrupting every component; fixed with a per-component `sum_comp`
  buffer).
- **N passes**: the number comes from `CONTROL_PRINT_DOF_SMOOTH_N`
  (default 10, `nsmooth<1` -> `db_error`). Each pass uses the PREVIOUS
  pass state (two buffers, `array_move` copy per pass).
- The node id column of `control_print_dof_id` is not smoothed (it is a
  mapping, not a dof value).

## Detalles

- The smoothing runs on the same data the write loop prints
  (VERSION_PRINT after `db_version_copy`/`renumbering`), so the printed
  values are the same ones smoothed.
- `smooth_field` is `(max_node+1)*nuknwn` doubles; freed when the
  record is not active (`smooth_field=NULL` -> raw values).

## Pendiente

- None. (Pre-existing GNU quirk, NOT touched: in 1D, `materi_stress`
  registers 6 matrix unknowns so `control_print_dof` in 1D writes 37
  component blocks, most of them with out-of-range reads; the smoothing
  tests use a 1D model WITHOUT `materi_stress` so only the velocity
  block is printed.)
