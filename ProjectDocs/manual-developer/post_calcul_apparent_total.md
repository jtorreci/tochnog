# post_calcul_apparent_total

## Implementation (F5 T2 V8)

- **Record:** `database.cc`, INTEGER switch, length 1, no index, class POST.
- **Consumer:** `calcul.cc`, `calculate_operat()`, branch
  `YOUNG_APPARENT` / `POISSON_APPARENT`. With `-yes` the previous state
  `prev_dof` is set to 0 (instead of the `NODE_DOF_PREVIOUS_STEP`
  snapshot), so the "increments" are the total stresses and strains. The
  rest of the branch is unchanged.

## Verification

`post_calcul_apparent_total1` (Professional `ground17` plus the switch):
Professional 02-08-2026 1.333333333334e+07 / 3.333333333338e-01, GNU
13333333.3333 / 0.333333333334 (targets in the input).
