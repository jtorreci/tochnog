# post_calcul_multiply

## Description

`post_calcul_multiply factor_0 factor_1 ...` (manual Professional 6.916). A
multiplication factor for each calculated item, for example to use another
definition of a result.

An item is one dof/operat pair of `post_calcul`; all values of a
multi-valued item (for example the three `-prival` values) share its
factor. Give a factor for every item; a shorter record is an error
("post_calcul_multiply needs a factor for each item of post_calcul.");
extra values are ignored.

The factor is applied first, before `post_calcul_limit` and
`post_calcul_absolute`.

## Usage

```
post_calcul -materi_strain_total -average -materi_strain_total -size_tot
post_calcul_multiply 2. 5.
```

## Example

`validation-suite/test-2026/post_calcul_modify1.dat` (factors 2, -3, 5) and
`post_calcul_modify_error1.dat` (two factors for three items).

## Differences with the Professional version

None measured. The Professional prints only "Error detected for data item :
post_calcul_multiply" for a short record; the GNU adds the reason.
