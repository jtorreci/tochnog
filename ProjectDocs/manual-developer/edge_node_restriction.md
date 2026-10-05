# edge_node_restriction — developer notes

## Where

`area.cc`:

- `side_node_restriction_ok( master, ind, element, nodes, sides, iside,
  nnol_side )` (new static helper, after `force_edge_is_master()`):
  returns 0 when the `_node` companion of `master` with index `ind` is
  active and a node of side `iside` is not listed, or when the
  `_element_node` companion is active and the element differs or a local
  side node is not listed.
- `area()`, side loop: after the geometry / `type_area` test of every side
  node, `if ( ok && !side_node_restriction_ok(...) ) ok = 0;` — the side
  is skipped as a whole, exactly like a side outside the geometry.

Companion lookup: `flux_edge_companion(master, 6/7)` (condif heat,
groundflow flux), `force_edge_companion(master, 3/4)` (force_edge,
_normal, _water, _projected; -1 when absent), `conv_rad_companion(master,
3/4)` (convection/radiation edge normal; the legacy `condif_convection` /
`condif_radiation` masters have no node companions), and
`support_edge_companion(3/4)`.

## Why

The per-family load code checked `array_member(nodes, inod)` for each
side node separately (`use_it`), so a side with only one listed node still
received that node's share. The Professional 02-08-2026 loads a side only
when all its nodes are listed (probes on a 10 x 10 quad4, nodes 1 3 vs 1
vs 1 3 4, for heat flux, `_element_node`, groundflow flux, convection,
`force_edge`, `force_edge_normal`, `support_edge_normal`; all equal the
Professional after the fix). Pro test `condif7`: temp 20 -> 10.

## Notes / pending

- The per-node `use_it` checks in the load code are kept; for a selected
  side they are always true now (redundant, harmless). Removing them is a
  possible clean-up.
- `_element_node_factor` (flux) and `_node_factor` (force) still index the
  factors by the ELEMENT local node number (`j==inol`); the manual says
  "factor0 is the factor for the first node on the side". Not verified
  against the Professional (no corpus example).
- Hardcoded: none.
