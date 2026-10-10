# Bi-CG direct-LU retry band width — developer notes

> **F4 T7 (2026-10-09):** the `NODE_NODE` band loop of the Bi-CG branch
> described below no longer exists. The retry builds its matrix lazily in
> `solve_band_retry()` with the EXACT band measured on the element matrix
> entries (node-number or Cuthill-McKee ordering, whichever is narrower);
> see [solver_memory](solver_memory.md). The history below is kept for the
> examp2 diagnosis.

## Where it is implemented

- `so.cc`, function `solve()`, the `bicg_solver` branch that computes the
  band width (`band = 0;` loop over `NODE_NODE`, ~line 517).
- Consumer: the direct-LU retry inside `solve()` after
  `solve_iterative_bicg()` returns unsuccessful (`dgbsv_` on `mat`, ~line
  958), and the band fill in the same function (~line 830).

## Implementation details

The band workspace (`mat`, `kl = ku = band`, `ldmat = 2*kl+ku+1`) is
allocated for every solver route because the Bi-CG failure path retries with
`dgbsv`. The band width for the Bi-CG route is computed from the node-to-node
connectivity (`NODE_NODE`) with the identity node ordering:

```c
long int diff = ( inod>jnod ) ? inod-jnod : jnod-inod;
if ( (1+diff)*nprinc > band ) band = (1+diff)*nprinc;   // fixed
```

The previous code used `diff*nprinc`. Two nodes separated by `diff` node
numbers occupy `diff*nprinc + (nprinc-1)` consecutive local indices: the last
dof of one node to the first dof of the other. The old bound was short by
`nprinc-1`, so for the first-to-last coupling of adjacent nodes the band
index `ii = kl+ku+1+ilocal-jlocal` exceeded `ldmat` and the write
`mat[indx-1] += val` landed in the next column (silent heap corruption).

The Cuthill-McKee ordering branch (`else` of the same `if (bicg_solver)`)
already used `itmp = (1 + labs(global_ordered[inod]-global_ordered[jnod])) * nprinc`,
i.e. the same `(1+diff)*nprinc` bound; the Bi-CG branch (added later for the
direct-LU retry) omitted the `+1`. The fix makes both branches agree.

### How it was diagnosed

1. `examp2` with `inertia_apply -yes` + `groundflow_consolidation_apply -yes`
   gave `pres = -0.7206` (default Bi-CG) but `-0.634519741323` with
   `options_solver -matrix_lapack` / `-matrix_superlu`, exactly the
   Professional `-0.6345197413352`.
2. A `git bisect` (build + run at each step, default solver) identified
   `6cafafe` (tighter Bi-CG defaults) as the first commit that changes the
   result, because it makes the Bi-CG break down and take the faulty retry.
3. Temporary instrumentation in `so.cc` printed
   `band=2 solve_nlocal=128 nprinc=2 npuknwn=14` and
   `OUT_OF_BAND ilocal=3 jlocal=0 ii=8 ldmat=7 diff=3`: the consolidation
   coupling (pressure at node i+1 to velocity at node i) is a dof distance 3
   while the band only allowed 2.

## External features / libraries

- LAPACK `dgbsv_` (band LU) and `dgesv_` (dense LU); declared in `so.cc`.
- SuperLU for `options_solver -matrix_superlu`.

## Hardcoded parameters

None.

## Pending refactorings

- Done in F4 T7: the retry band is derived exactly from the element matrix
  unknown pairs (`solve_band_retry_entries`).
- The out-of-bounds write was silent. A bounds assertion in the band fill
  (debug builds) would have caught it earlier.
