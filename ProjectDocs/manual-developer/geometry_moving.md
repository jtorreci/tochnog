# geometry_moving

`geometry_moving`, `geometry_moving_parameter`, `geometry_moving_operat`,
`geometry_moving_operat_parameter`, `geometry_moving_operat_time`,
`geometry_moving_n` (manual Professional 6.533-6.538) and
`control_geometry_moving` (6.135); F5 T4 U8. None was registered before.

## Files and functions

- `geometry_moving.cc` (new module, in `makefile`):
  - `geometry_moving_initialise(first_of_control)`: called from
    `step_start()` (top.cc) on the first step of a control index with an
    active `control_geometry_moving` (value `-initialise`, anything else is a
    `db_error`). Reads and checks every `geometry_moving` index
    (`gm_read`), removes the previous `element_geometry_moving*` records and,
    per entity, writes a scratch geometry record of the entity at its start
    position (`gm_put_scratch`) at an index no `geometry_*` record uses
    (`gm_scratch_index`), so that no per-index option of a user geometry
    (eps_iso, node type, filters) applies. For every element:
    `gm_element_samples` (sample points and weights) and
    `gm_element_entity` (time loop). The scratch record is deleted at the end.
  - `gm_element_entity`: for the `ntime+1` moments
    `t_i = t_start + i*(t_end-t_start)/ntime` the sample points are moved
    back over the translation (`x - v*(t_i - t_start)`) and tested with
    `geometry(-1, x, {entity, scratch}, ..., PROJECT_INSIDE)`; a bounding
    box prefilter (entity box at the start plus a margin of radius + 1 % of
    its size) skips the moments the element cannot touch. Per point: first
    moment inside (`first`), last moment inside before an outside one
    (`leave_index`). Entry time = `t_{first-1}` (`t_0` when inside at `t_0`),
    leave time = `t_{leave_index}`: the Professional 02-08-2026 convention,
    measured on its `geometry_moving1` (element 6: nodes enter at 1.063,
    1.080, 1.132 -> in_time 1.05 1.125; leave at 1.183, 1.253 -> out_time
    1.17 1.245, sample step 0.015). The part is the weighted fraction of the
    points ever inside.
  - `gm_element_samples`: tensor-product elements (`-bar2/3/4`, 2-node
    `-truss/-beam/-trussbeam`, `-quad4/9/16`, `-hex8/27`, nodes numbered with
    the first direction fastest) use the midpoints `-1 + (2i+1)/nspace` of
    `nspace` cells per isoparametric direction, the Lagrange shape functions
    of `interpolation_polynomial()` (polynom.cc) and the weight
    `sqrt(det(J^T J))` (also for a line element in 2D/3D); simplices
    (`-tria3/6`, `-tet4/10`, corner nodes 0 1 2 / 0 2 5 / 0 1 2 3 / 0 2 5 9 as
    `element_volume_set`) the collapsed (Duffy) map of the unit square / cube
    cell midpoints, weights `1-v` / `(1-v)(1-w)^2` (straight-sided only);
    every other element (prism, interface, spring, contact spring) its nodes
    with equal weights (the Professional rule).
  - Several entities: an element keeps the entity with the largest part
    (the lower index on a tie; `ELEMENT_GEOMETRY_MOVING_PART` read back).
  - `geometry_moving_apply(time_current)`: every `step_start()` (after the
    initialisation, before `delete_geom`). For an element with
    `ELEMENT_GEOMETRY_MOVING_PART` and `time_current > in_time[0]`: ramp
    `min(1, (t - in0)/(in1 - in0))` (1 when `in1 == in0`), factor
    `1 - part*ramp`; above `1.e-6` it is PUT as `ELEMENT_DELETE_FACTOR` (the
    elem.cc frozen released-stress path), otherwise `delete_element()` and,
    once, `mesh_has_changed()` (removes the unused nodes). The linear ramp
    was measured on the Professional (element 6: 0.9 0.7 0.5 0.3 at 1.06 1.08
    1.10 1.12, 0.25 from 1.125; element 5 0.4167 / 0.0833 at 1.04 / 1.06, gone
    at 1.08). Before the PUT it allocates `ELEMENT_RHSIDE_DELETE` up to the
    highest element (`db_allocate -MINIMAL`): elem() stores it inside the
    parallel element loop; without the allocation `geometry_moving_ramp1`
    and `_band2` aborted with "Data is allocated in a parallel loop" (the
    Professional has the same abort for its example with threads). The
    `delete_geom` path does not abort on the same meshes (delete3, a 2 x 2
    probe with `control_mesh_delete_geometry_factor` at a later index); why
    it escapes was not isolated.
- `top.cc` `step_start()`: the two calls before `delete_geom()`.
- `database.cc`: the six input records and `control_geometry_moving`
  (class CONTROL: not a geometry entity for the `-geometry_* index`
  references), the four element records (class ELEMENT, printed, removed by
  `delete_element`), the words `initialise` and `translate`.
- `tochnog.h` / `tochnog-mod.h`: the enums (before `XY`) and the two
  prototypes.

## Hardcoded parameters

- `GM_NTIME_DEFAULT 100`, `GM_NSPACE_DEFAULT 5` (manual 6.538).
- `GM_EPS_FACTOR 1.e-6` (as `EPS_ELEMENT_DELETE_FACTOR`, delete.cc),
  `GM_EPS_TIME 1.e-12`, `GM_BOX_MARGIN 1.e-2` (prefilter only).
- Triangle and quadrilateral scratch records get the tolerance `EPS_COORD`
  (a filled area in 2D; in 3D a surface excavates practically nothing).

## Verification

- Professional 02-08-2026, `test/other/geometry_moving1.dat` with
  `processors_maximum -no` (with threads it aborts): the records, the time
  conventions and the ramp above; the part is the fraction of nodes ever
  inside (`nspace` 2, 5, 200 give identical records; a point passing between
  the nodes excavates nothing) -> DIVERGENCES GEOMETRY-MOVING-SAMPLING.
- Analytic: `geometry_moving_band1` (part 0.4 and all entry / leave times by
  hand), `_band2` (deletion, nodes), `_simplex1` (0.7496 vs 0.75),
  `_hex1` (0.5 exact, tet 0.8758 vs 0.875 with the `point_el` tolerance),
  `_ramp1` / `_ramp2` (delete3 stresses 0.6 / 0.76), `_error1`-`_error5`.
- `geometry_moving1` (copy of the Professional example) is not in the
  build_safe suite: 2000 steps take more than 10 minutes (the bicg solver
  breaks down on the zero-stiffness nodes of the partly excavated elements
  and every iteration falls back to SuperLU, BACKLOG C130). The full run
  (35 min) gives sigyy -0.348912581336 against the Professional
  -0.3268 +- 0.02, a consequence of GEOMETRY-MOVING-SAMPLING; owner
  decision 2026-11-09: class D, the copy targets the GNU value (1.e-6).

## Pending refactorings

- The sampling rules (cell midpoints, Duffy map) are local to the module;
  a shared iso -> global mapping helper for all element types would remove
  the node fallback for prisms.
- Curved simplices (`-tria6`, `-tet10`) are sampled on their corners.
- Only `-translate` exists (manual); a rotation operator would need the
  scratch record to move instead of the points.

## Differences with the Professional version

- GEOMETRY-MOVING-SAMPLING (class D): `gm_element_samples` follows the
  manual sample points; the Professional uses the nodes (the node fallback
  of the module is its rule).
- GEOMETRY-MOVING-PARAMETERS (class D): `gm_read` takes the manual lengths
  (`nvertex*ndim` + radius, `ndim` velocities).
- GEOMETRY-MOVING-THREADS (Professional defect): the `ELEMENT_RHSIDE_DELETE`
  allocation in `geometry_moving_apply`.
- GEOMETRY-MOVING-CONTROL-VALUE (class C): `-initialise` only.
