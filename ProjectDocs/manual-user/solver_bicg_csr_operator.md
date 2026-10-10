# solver_bicg_csr_operator (assembled matrix of the iterative solver)

## What it does

The default iterative solver of the GNU (`-matrix_iterative_bicg`: CG for
symmetric systems, Bi-CG otherwise) multiplies the system matrix with a
vector on every iteration. Until F4 T9 it did that element by element,
re-reading every element matrix on every iteration. Since 2026-10-10 the
element matrices are assembled once per solve into one sparse matrix
(compressed sparse rows) and the iterations multiply with it.

## Why it is useful

The iterations are about five times faster on 3D models (Professional
`force11`: ~5000 iterations in the time of ~1000 before) and every model
that uses the iterative solver gains. The results are the same up to
round-off (the same sums in another order); the number of iterations
is unchanged.

## Input

Nothing to set: it is automatic for `-matrix_iterative_bicg` (the
default), and for the Bi-CG fallback of `-matrix_pardiso`.

- `options_processors <n>` (or `processors <n>`) splits the products in
  `n` row blocks, one per thread.
- `print_solver -yes` prints, per solve, `iterative solver: assembled CSR
  operator, rows <n> nonzeros <nnz>` (or `element-by-element operator`
  when the fallback below is used).
- The element-by-element product is still used for systems with the
  `mpc` tie elimination (small tied systems, DIAG-SOLVE-MIXTO §16), and
  when the environment variable `TN_BICG_EBE` is set (diagnostics).

## Example

```
options_processors 2
print_solver -yes
```

Test: `validation-suite/test-2026/bicg_csr_threads1.dat` (3D hex8 block,
two threads, analytic targets; build_safe.sh checks the CSR line).

## Differences with the Professional version

None in the results: the Professional default solver is PARDISO (see
[matrix_pardiso](matrix_pardiso.md)); this is an internal speed-up of
the GNU iterative solver.
