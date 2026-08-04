# check_solver

## Descripción

`check_solver` makes the solver inspect the diagonal terms of the matrix before
solving. If the absolute value of a diagonal term is smaller than the given
tolerance `eps`, a warning is printed with the equation number. A small diagonal
term usually indicates a problem in the input file, such as an undefined
stiffness or incorrect boundary conditions, and would otherwise lead to a
singular or nearly singular matrix.

The check only applies to the direct LAPACK band solver
(`options_solver -matrix_lapack`). It is disabled by default.

## Uso

Place it in the data part, as a keyword line:

```
check_solver 1.e-6
```

## Parámetros

| Parameter | Meaning                                                             |
|-----------|---------------------------------------------------------------------|
| `eps`     | Tolerance for the matrix diagonal. Diagonal terms with `|diag| < eps` trigger a warning with the equation number. |

## Ejemplo

Minimal input that enables the LAPACK band solver and checks its diagonal:

```
control_geometry
   cartesian
control_time
   end 1.0
   dt 0.1
control_print
   history 0
   step 10
options_solver
   -matrix_lapack
materi_elasti_young
   0 2.e7
materi_elasti_poisson
   0 0.3
materi_density
   0 2500.
check_solver
   1.e-6
```

For every equation whose diagonal term is smaller than `eps`, the output
contains a warning similar to:

```
Warning: small diagonal term in equation <n> (<value> < eps=1e-06).
This normally indicates a problem in the input file.
```
