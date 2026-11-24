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

## FRD and gmsh printers (BACKLOG C24, F1 T1, 2026-10-05)

- Files `print_fr.cc` (`print_frd()`, result blocks) and `print_gm.cc`
  (`print_gmsh()`, `$NodeData` / `$ElementData` / `$ElementNodeData`
  blocks) had the same `ipuknwn<nuknwn_` loop. Both now loop over the
  primary unknowns (`npuknwn_ = nuknwn_/nder_`, step `nval`, `iuknwn =
  ipuknwn*nder_`, type, label and dof type read at `iuknwn`).
- `print_frd()`: an FRD vector block always has 3 components; in 2D the
  third (`D3`, `V3`, ...) is now written as 0 instead of reading
  `iuknwn + 2*nder_` (the next dof). Result names and component names
  are unchanged (`VELO`, `DISP`, `STRESS`/`SXX..`, `TOSTRAIN`/`EXX..`,
  `NDTEMP`; other dofs: dof label truncated to 8 characters).
- `print_gmsh()`: the field names are now built from the dof TYPE
  (`DOF_TYPE`, newly read): `node_materi_velocity_0`,
  `node_materi_stress_xx`, `node_condif_temperature`, ... (same base
  names as the Professional). Before, they used `dof_label[ipuknwn]`
  with the wrong index (`node_velx_1` for the y velocity, `node_sigxx_yy`
  ...).
- Verification: `validation-suite/test-2026/print_frd_gmsh_dof_loop1.dat`
  (single quad4 with `derivatives`, velocity, displacement, stress,
  total strain, temperature). Old binary: 68 FRD `-4` blocks and 324
  gmsh `node_` fields; new: 5 FRD blocks (`VELO DISP STRESS TOSTRAIN
  NDTEMP`) and 17 `node_` fields, no `e-3xx` garbage. Checked by
  `scripts/build_safe.sh` (the `.frd` is appended to when it exists, so
  the check removes it before the run). The values equal the FIRST block
  of each dof of the Professional 02-08-2026 (`velx` 0.01, `disx`
  0.01, `sigxx` 10, `epsxx` 0.01, `temp` 0/1).
- Remaining layout differences with the Professional (BACKLOG C27, not
  redesigned here): the Professional writes one `1PSTEP`/`100CL` header
  per result block, vector blocks with an extra `ALL` component, the
  names `VELOCITY`, `S11..S13`, `T11..T13`, extra non-dof blocks
  (`BOUNDA_*`, `BOUNDARY`, `REACTION`) and, with `derivatives`, one
  block per value AND derivative entry (whose values are themselves
  shifted, the same defect class as C24). Its gmsh file writes vector
  and tensor fields as one 3/9-component field (`node_materi_velocity`,
  `node_materi_stress`) plus bounda/boundary/rhside/group/flux fields;
  the GNU writes one scalar field per component.

## Pending

- BACKLOG C27: FRD/gmsh layout differences listed above.

## Differences with the Professional version

- **`control_print_dof -separate_index` file layout** (known difference, not
  decided; BACKLOG C26). GNU: One `dof.N` file, components one after
  another, short `%g`. Professional 02-08-2026: One file per dof (`disx.21`,
  `sigxx.21`, ...), `%e` with 12 digits. Status: Open, LOW; follow the Pro
  layout or keep it documented. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **FRD / gmsh output layout** (known difference, not decided; BACKLOG C27).
  GNU: FRD: no PSTEP per block, 3-component vectors, names `VELO`, `SXX`..;
  gmsh: one scalar field per component, nodal averages. Professional
  02-08-2026: FRD: `1PSTEP` per block, `ALL` component,
  `VELOCITY`/`S11`/`T11`, extra boundary/reaction blocks, and with
  `derivatives` value AND derivative blocks shifted into other dofs (Pro
  defect, not to be copied); gmsh: multi-component fields, integration-point
  element fields. Status: Open, LOW; decide per format which Pro conventions
  to follow. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
