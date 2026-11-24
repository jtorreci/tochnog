# group_materi_plasti_hypo_wolffersdorff_niemunis

## Description

Niemunis extension of the von Wolffersdorff hypoplastic law (manual
Professional 6.719; Niemunis 2003, "Extended hypoplastic models for soils",
section 4.3). Two changes on top of `group_materi_plasti_hypo_wolffersdorff`:

- **Shear stiffness.** The reference law corresponds to a Poisson ratio of
  about 0.38 at the hydrostatic axis, which is too soft in shear. A term
  `b2 (I - 1/3 1 1)` is added to the linear stiffness, with `b2` chosen so
  that the hydrostatic response has the Poisson ratio `nu`:
  `(1 + a^2/3 + a/sqrt(3)) / (1 + b2) = (1 + nu) / (1 - 2 nu)`. The flow
  rule is kept, so the nonlinear part follows the new stiffness.
- **Consistency at e = e_d.** The density factor becomes
  `f_d = sign(e - e_d) |r_e|^alpha`, `r_e = (e - e_d)/(e_c - e_d)`, and if
  `f_d < 1`: `f_d := f_d + (1 - f_d)^z f_d_bar`, with `f_d_bar` the value
  that keeps the state above the lower bound `e_d` (eq. 4.221). This keeps
  dense states from passing below `e_d` in extension.

## Syntax

```
group_materi_plasti_hypo_wolffersdorff_niemunis index nu z
```

`nu`: Poisson ratio at the hydrostatic axis, typically 0.2 or 0.3; it must
be larger than -1 and smaller than 0.5. `z`: smoothing exponent of the f_d
interpolation, typically 5.

## Example

```
group_materi_plasti_hypo_wolffersdorff 0 33. 1000. 0.25 0.95 0.55 1.05 0.25 1.
group_materi_plasti_hypo_wolffersdorff_niemunis 0 0.2 5.
```

Tests: `validation-suite/test-2026/group_materi_plasti_hypo_wolffersdorff_niemunis1`-`3`
(`1` is the Professional test `other/hypo11`), `_error1`.

## Differences with the Professional version

- `nu` of 0.5 or more (b2 of -1 or less) is an error; the Professional
  02-08-2026 runs and diverges (sigyy -5.5e19, PARDISO error -4).
- Only the von Wolffersdorff route reads the record (not lowangles, not the
  Masin routes).
