# check_element_shape

## Description

`check_element_shape factor` checks the Jacobian distortion of isoparametric
elements during the integration loop. For each element it computes the average
of `|det_j - det_avg| / det_avg` over the integration points. If that average
distortion exceeds `factor`, a warning is printed to the standard output.
Heavily distorted elements usually point to a bad mesh: inverted or nearly
degenerate elements that degrade accuracy and may stall convergence.

## Usage

Place it in the data part, as a keyword line:

```
check_element_shape 0.01
```

## Parameters

| Parameter | Meaning                                                    |
|-----------|------------------------------------------------------------|
| `factor`  | Maximum allowed average Jacobian distortion. A value `>0` enables the check; `0` (default) disables it. |

## Example

Minimal input that runs the check on an elastic block:

```
control_geometry
   cartesian
control_time
   end 0.1
   dt 0.01
control_print
   history 0
   step 10
materi_elasti_young
   0 2.e7
materi_elasti_poisson
   0 0.3
materi_density
   0 2500.
check_element_shape
   0.01
```

Run it; for every element whose average Jacobian distortion exceeds the
tolerance the log will show a line like
`Warning: element with distorted shape, distortion=... > factor=...`.
