# element_intpnt_h / _iso_coord / _method / _npoint (developer)

## Where

- `polynom.cc`, `element_intpnt_set()`: called at the end of `top.cc`
  `step_start()` (serial, after every mesh operation of the step and
  before the parallel element loop, where db allocation is not allowed;
  also on the `step_start( NO, ...)` path of a control index without
  time steps). For each active element:
  - interface group (`GROUP_INTERFACE` active):
    `element_intpnt_interface()` builds the layout of
    `interface_element()` (`interface.cc`): one point per facing pair
    (i, i+ns1), h = 1/2 on both nodes; 2D ns1 = 2..4 and 3D ns1 = 4 / 9
    Lobatto points (`integration_lobatto()`), method `-LOBATTO`; 3D
    ns1 = 3 (prism6) the triangle corners, method 0. Other interface
    shapes (bar2 with ns1 = 1, ...) get no records.
  - otherwise only the shapes `pol()` integrates, with the matching
    dimension and node count (`element_intpnt_supported()`: bar2/3/4 in 1D;
    tria3/6, quad4/9/16 in 2D; tet4/10, prism6/15, hex8/27 in 3D): one
    `pol()` call on the `NODE` coordinates (VERSION_NORMAL, old = new), then
    `npoint`, `h`, and the two new thread-local outputs of `pol()`.
- `polynom.cc`, `pol()`: `pol_method_used` (the final
  `integration_method` of the tensor branch, after the SRI and C71/C80
  switches; 0 for the triangle, tetrahedron and prism branches) and
  `pol_iso_used[MPOINT*MDIM]` (the point coordinates, point-major), both
  `static thread_local`, written by every `pol()` call and read only by
  `element_intpnt_set()` right after its own call. Chosen over a capture
  pointer argument: no signature change for the other `pol()` callers
  (`elem.cc`, `calcul_force.cc`) and no shared state across threads.
- `database.cc`: `ELEMENT_INTPNT_H` and `ELEMENT_INTPNT_ISO_COORD`
  DOUBLE_PRECISION, `DATA_ITEM_SIZE`, variable length; `ELEMENT_INTPNT_METHOD`
  and `ELEMENT_INTPNT_NPOINT` INTEGER, length 1; all `print_only`, class
  ELEMENT. Enums appended before `LAST_DUMMY` in `tochnog.h` and
  `tochnog-mod.h`; prototype `element_intpnt_set()` in both headers.

## Implementation details

- No effect on the results: `pol()` resets its thread-local recovery
  state (`stress_recovery_*`, `lumped_hrz_active`) at every call and the
  element loop calls `pol()` again before any reader; the element shape
  check (`check_element_shape_factor`) is set to 0 during the serial call
  (same pattern as `pef_element_internal_forces_pol()`, `calcul_force.cc`),
  so no warning is printed twice.
- Iso coordinates per branch: tria3 `-minimal` (1/3,1/3), `-maximal` the
  nodes (1,0),(0,1),(0,0); tria6 (L1,L2); tet10 (L1,L2,L3); tet4 (L1,L3,L2)
  because the tet4 branch puts L3 on node 2 and L2 on node 3 (the
  Professional convention is coordinate k = h of corner k; the `-maximal`
  points are listed in the Professional order, point k+1 with 1/2 on
  corner k, since 2026-11-02: TET4-MAX-ORDER resolved); prism6
  (L1, L2, 2 zeta - 1) because the GNU prism6 reference cell has zeta in
  [0,1] and the Professional prints [-1,1]; prism15 (L1,L2,zeta); tensor
  elements (xi, eta, zeta) up to ndim, xi fastest.
- Measured on the Professional 02-08-2026 (omen probes, every element type
  of the tests): the records exist in every run with no switch; truss,
  beam and spring elements get none; tria3 / tria6 / tet4 / tet10 / prism6
  and the prism6 interface print method 0; `-reduced` prints `-reduced`.
- The axisymmetric C71 / C83 / C80 switches use the radius of the
  undeformed `NODE` coordinates (the element loop uses the current ones);
  equal unless a large-displacement mesh moves a node onto the axis.

## Hardcoded parameters

- `DATA_ITEM_SIZE` (190): `element_intpnt_h` is skipped (and an older
  record deleted) when npoint x nnol is larger (hex27, prism15, quad16
  `-maximal`); the Professional writes it. Writing it would need a longer
  record, but generic readers (`print_data_versus_data`, filters) use
  `DATA_ITEM_SIZE` buffers.
- The interface layouts (ns1 ranges, the prism6 corner order) are the
  measured Professional ones, not derived from a rule table.

## Pending

- `element_intpnt_h` for the rules above 190 values (DIVERGENCES
  INTPNT-H-LENGTH, class B).
- Elements created after the last `step_start()` (refinement at the end
  of a step) have no records until the next step.
- Tests: `element_intpnt_quad4_1`, `_tria3_1`, `_interface_quad4_1`,
  `_1d1`, `_2d1`, `_3d1`, `_interface2d1`, `_interface3d1`
  (`scripts/build_safe.sh`, `intpnt_check`).
