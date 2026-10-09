# Solver memory (default Bi-CG solver)

## What it does

The default iterative solver (Bi-CG) no longer allocates direct-solver
workspaces that it does not use. Before F4 T7 every solve allocated and
zero-filled a dense `n x n` copy of the system matrix (8 n² bytes: 20 GB
for 50 000 unknowns), although nothing ever read it.

## Why it is useful

Large models (tens of thousands of unknowns) ran out of memory or were
killed by the 4 GB process cap long before the solver started, and every
solve of every model paid the O(n²) zero-fill. The results of the solver
are unchanged: the removed matrix was never an input of any solve.

## Input syntax

No keyword. The behaviour is automatic for the default solver
(`-matrix_iterative_bicg`, and `-matrix_pardiso`, which falls back to it).

## Parameters

None.

## Minimal example

Any input that uses the default solver, for example the Professional
`test/other/ground19_water_under_dam.dat`: before the fix it stopped with
`std::bad_alloc` under a 4 GB cap (5.2 GB RSS without the cap); it now
runs.

## Differences with the Professional version

None in the results. The Professional uses its own direct solver
(PARDISO).
