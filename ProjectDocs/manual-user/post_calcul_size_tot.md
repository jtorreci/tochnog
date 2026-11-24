# post_calcul_size_tot

## Description

`post_calcul` operators `-size_tot` and `-size_dev` (GNU spellings
`-sizetot`, `-sizedev`) give the size of a calculated matrix or vector
(manual Professional, `post_calcul`):

- `-size_tot`: sqrt(Aij Aij) for a matrix, sqrt(Ai Ai) for a vector.
- `-size_dev`: sqrt(Bij Bij) with Bij = Aij - delta_ij (A11+A22+A33)/3
  the deviatoric part (matrices only).

The result is stored as the scalar item `s<name>` / `sd<name>` of
`node_dof_calcul` (e.g. `smateri_stress`, `sdmateri_stress`).

Since F4 T6a (2026-10-11) these are the manual definitions. Before, the
GNU used sqrt(0.5 Aij Aij) and sqrt(0.5 Bij Bij), 1/sqrt(2) of the
manual value (e.g. 0.707107 instead of 1 for a unit strain).

## Usage

```
post_calcul -materi_strain_total -size_tot -materi_strain_total -size_dev
```

## Example

Bar with a strain of 1 in x (Professional `calcul1.dat`,
`validation-suite/test-2026/post_calcul_size_tot1.dat`): `-size_tot`
1, `-size_dev` 0.816497.

## Differences with the Professional version

None for the sizes. The `-prival` order still differs (DIVERGENCES K12).
Inputs written for the old GNU values (the sfnet `calcul1`/`calcul2`
targets 0.707 / 0.577) are stale.
