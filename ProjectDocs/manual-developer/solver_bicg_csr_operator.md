# solver_bicg_csr_operator — developer notes (F4 T9, 2026-10-10)

## Where

`so_bicg.cc`:

- `bicg_csr_build()` — builds the operator at the start of
  `solve_iterative_bicg()` (after the solve setup, before the first
  product); `bicg_csr_free()` at its end.
- `bicg_csr_visit_element()` — template visitor of the matrix entries of
  one element mapped to solver dofs, with the same
  `ELEMENT_MATRIX_VALUES` / `GROUP_MATRIX_VALUES` selection as
  `solve_iterative_bicg_element()`.
- `bicg_csr_sys()`, `bicg_csr_rows()`, `parallel_solve_iterative_bicg_csr()`
  — the products; `solve_iterative_bicg_sys()` calls `bicg_csr_sys()`
  when `bicg_csr_on`.

Owner decision 2026-10-10 (BACKLOG C43 item 4): assembled CSR operator
by default, element-by-element fallback, threads, differences within the
Bi-CG tolerance.

## Pattern without a triplet buffer

1. Node set of every element: the nodes of its entries whose row and
   column are solver dofs (marker array, one pass over the entries).
2. Node -> element incidence (counting sort of the node sets).
3. Rows: every solver dof of node `a` gets the solver dofs of all nodes
   sharing an element with `a` (union of the node sets of its elements,
   marker array), sorted; the list is computed once per node and copied
   to each of its dof rows (`rowptr` in local dof order).
4. Values: every element entry is added at its column, found by binary
   search in its row (`std::lower_bound`). An entry outside the pattern
   (impossible by construction) aborts the build and the solve uses the
   element-by-element products.
5. `NODE_LHSIDE` is kept apart as a diagonal vector (`bicg_csr_lhside`),
   with `bicg_csr_bfac` = 1 for the dofs of nodes with `NODE_LHSIDE`: the
   element-by-element residue adds the right-hand side only in that node
   loop, and the CSR residue reproduces it.

Memory: `rowptr` (8 bytes per row), `col` (int) + `val` (double) = 12
bytes per nonzero, plus temporary node sets and incidence (a few longs
per element). Professional `force11`: 31104 rows, 2.16 M nonzeros
(26 MB).

## Products

Per row `i`: `Ad1 = p_i sum_j A_ij p_j d1_j + lhs_i p_i^2 d1_i`,
`residue = bfac_i p_i b_i - p_i sum_j A_ij p_j x_j - lhs_i p_i^2 x_i`,
`p_tmp = A_ii + lhs_i` (1 when below `EPS_P`); for the Bi-CG the
transpose product `Ad2_j += A_ij p_i p_j d2_i` is scattered into the
per-thread buffer `Ad2_thread` and reduced with `lhs_j p_j^2 d2_j`. The
same quantities as the element-by-element products; only the summation
order changes (round-off).

Threads: `options_processors` row blocks through `parallel_sys_routine`
/ `parallel_sys_next_of_loop` (`bicg_csr_nblock` = number of threads;
the serial `sysother` build runs all blocks in one call). Ad1, residue
and p_tmp are written per row (no reduction), so they do not depend on
the thread count.

## Fallbacks (element-by-element products)

- `mpc_tie_elimination_active()` (the tie expansion of DIAG-SOLVE-MIXTO
  §16 maps slave entries to master pairs; only small systems,
  `MPC_TIE_TIGHT_NMAX`);
- an entry not found in the pattern;
- environment variable `TN_BICG_EBE` (diagnostics, A/B of the operator).

## Hardcoded parameters

- none new (`EPS_P` diagonal guard as before).

## Verification

- `bicg_csr_threads1` (test-2026): CSR line in stdout (old GNU: none),
  two threads, analytic targets; same iteration count as the
  element-by-element build (75 residual lines).
- build_safe 177 OK / 0 FALLO. Suite `.dbs` vs the element-by-element
  build (F4 T9 b3): 99 of 300 files differ at round-off; node_dof
  difference relative to the field maximum (components of one quantity)
  <= 6e-6 except fields of round-off size (|field| < 1e-6 of the largest
  field of the file, or stress fields ~1e-9 of the support tests
  `tsup_damp`/`tsup_auto`): `dpoint_time` 2.1e-5 and `msf_cant3d_hex27`
  1.4e-5 of tiny-load models where the absolute floor
  `bicg_error_minimum` binds.
- Corpus A/B on omen (element-by-element build `3349515` -> CSR, 4 GB,
  timeout 120, 6 jobs): Pro other 274 -> 276 OK (`ground11_nonsaturated`
  and `ground11_phreatic_level` now finish: 102 s / 82 s, were over the
  120 s slot), sfnet 192 -> 192, t26 111 -> 111, vs 232 -> 232. Wall time
  of the tests finished on both: Pro other 696 -> 545 s (-22 %), sfnet
  288 -> 220 s (-24 %), t26 37 -> 33 s, vs 28 -> 24 s. `.dbs` differ at
  round-off (Pro other 111, sfnet 78, t26 30, vs 82 files); node_dof
  relative to the field maximum <= 1e-5 except: `mesh_activate1` 3.3e-3
  (eptyy), `distri3` 1.8e-4, `update2` 9e-6, `truss13` 5e-6, sfnet
  `examp16` 1.4e-2, `examp10` 2.8e-3, `taylor3` (stalled test, rc 1 on
  both, kap 0). Explained as round-off amplification: the same inputs
  with the element-by-element build and `options_processors 2` (only the
  summation order changes) differ from the one-thread run by a similar or
  larger amount (`mesh_activate1` 2.8e-2, `distri3` 1.4e-4, `update2`
  5.7e-4, `truss13` 5.2e-6, `examp16` 5.0e-3); on `examp10` the first
  divergent solve follows a solve stopped by the absolute floor
  (`bicg_error_minimum` 1e-16 on the squared residual, i.e. a relative
  residual of ~2e-4 there), and a tighter `solver_bicg_error` does not
  change the difference (the floor binds). These are path-dependent
  plasticity / activation models; no class changes besides the two
  timeouts.

## Pending refactorings

- `solve_iterative_bicg_symmetric()` still walks the element matrices
  once per solve; it could test the CSR matrix.
- The pattern is rebuilt on every solve (the mesh rarely changes); it
  could be cached per mesh version.
- The SuperLU and band-retry fills could reuse the CSR pattern.
