# print_dof_component_loop

Fix of the dof-component loop of `control_print_dof` (BACKLOG C23, bug
fix, F0 T6, 2026-10-05).

## Implementation

- File `print_hi.cc`, function `print_dof()`, the write loop over the
  dofs.
- The node dof record holds `nuknwn = npuknwn*nder` values: primary
  unknown `ipuknwn` starts at `iuknwn = ipuknwn*nder` (value, then its
  derivatives). `DOF_SCAL_VEC_MAT` has `nuknwn` entries and `input.cc`
  sets it for EVERY value and derivative entry of a dof (`array_set(...,
  n*nder)`), so a vector of `ndim` components spans `ndim*nder` entries
  and the stress matrix `6*nder`.
- Old loop: `for ( ipuknwn=0; ipuknwn<nuknwn_; ipuknwn++ )`, testing
  `dof_scal_vec_mat[ipuknwn]` and reading `node_dof[ipuknwn*nder_ + ...]`.
  Two defects:
  1. every entry of a vector/matrix dof started a new block, so the
     components were printed `nval*nder` times, shifted by one entry
     each time (into the next dof);
  2. `ipuknwn*nder_` (plus `k*nder_` or `stress_indx*nder_`) went past
     the end of the `NODE_DOF` record: the values of the next nodes and,
     for the last nodes, heap memory past the database array (garbage
     denormals such as `6.9e-310` / `2.7e-312`, different every run).
     With `derivatives` (`nder = ndim+2`) the overrun is large
     (`dofid_no`: index up to ~210 for a 48-value record); without
     derivatives it is small (the last matrix entries) but present.
- Fix: loop over the PRIMARY unknowns as `print_vt.cc` does:
  `npuknwn_ = nuknwn_/nder_`, `for ( ipuknwn=0; ipuknwn<npuknwn_;
  ipuknwn+=nval )`, `iuknwn = ipuknwn*nder_`, type read from
  `dof_scal_vec_mat[iuknwn]`. Each vector/matrix dof is printed once,
  per component (`iuknwn + k*nder_`, `iuknwn + stress_indx(kk,ll)*nder_`).
  This is the layout documented in the user manual
  (`control_print_dof.md`: one line per node per dof component).
- Not a `get_new_dbl()` initialization defect (the expected pattern of
  `a32ea74`): the garbage came from an out-of-bounds read of the
  database array, so zero-initializing a buffer would not have fixed it.

## Verification

- `dofid_no` (test-2014, 2D, `derivatives`, 12 components x 4 nodes):
  `dof.21` 768 lines (20 with `e-3xx` garbage, run-dependent) -> 48
  lines, two runs byte-identical, no `e-3xx`.
- `build_safe.sh`: the `dof.21` check now also asserts 48 lines and no
  `e-3xx`; the `freq_timeint` count is 40 lines per print (10
  components x 4 nodes; was 176 with the repeated blocks).
- Corpus A/B old vs new: see `SEGUIMIENTO-CONVERGENCIA.md` (F0 T6 row).

## Notes

- The smoothed path (`control_print_dof_smooth_dof`) uses the same
  `indx`, so it is fixed by the same loop change.
- No new records, enums or hardcoded parameters.

## Pending

- The same `ipuknwn<nuknwn_` loop pattern exists in `print_fr.cc`
  (`print_frd`, result blocks) and `print_gm.cc` (`print_gmsh`, node
  data): BACKLOG C24, not fixed here because the change alters the
  block structure of those formats and their checks.
