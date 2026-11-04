# post_calcul_limit

## Description

`post_calcul_limit lower_0 upper_0 lower_1 upper_1 ...` (manual Professional
6.905). The lower and upper allowed value of each calculated item: a result
below `lower_i` becomes `lower_i`, a result above `upper_i` becomes
`upper_i`.

An item is one dof/operat pair of `post_calcul`; all values of a
multi-valued item (for example the three `-prival` values) share its pair.
Give a pair for every item; a shorter record is an error ("post_calcul_limit
needs a lower and an upper value for each item of post_calcul."); extra
values are ignored.

The limit is applied after `post_calcul_multiply` and before
`post_calcul_absolute`.

## Usage

```
post_calcul -materi_strain_total -average -materi_strain_total -size_tot
post_calcul_limit -0.5 0.5 0. 10.
```

## Example

`validation-suite/test-2026/post_calcul_modify1.dat` (-0.666667 limited to
-0.5), `post_calcul_modify2.dat` (limit before absolute) and
`post_calcul_modify_error2.dat` (one pair for three items).

## Differences with the Professional version

None measured. The Professional prints only "Error detected for data item :
post_calcul_limit" for a short record; the GNU adds the reason.
