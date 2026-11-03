# element_dof_initial_specific_number / element_dof_initial_specific_value

## Implementation

- **Records** (`database.cc`, enums appended in `tochnog.h` /
  `tochnog-mod.h`, F5 T2 V6):
  - `ELEMENT_DOF_INITIAL_SPECIFIC_NUMBER`: INTEGER, length 1,
    fixed_length 1, class ELEMENT. It holds a dof label.
  - `ELEMENT_DOF_INITIAL_SPECIFIC_VALUE`: DOUBLE_PRECISION, up to 1+MDIM
    values, class ELEMENT.
- **Consumption** in `general()` (`general.cc`), right after the
  `element_dof_initial` block. It runs only when that record is not
  active and `ELEMENT_DOF_INITIAL_APPLIED` is not yet set:
  - The label is mapped to `edis_iuknwn` with `array_member()` over
    `DOF_LABEL`.
  - The field `value_0 + sum grad_j * x_j` is evaluated at each node
    (`db_dbl(NODE, ...)`) into `edis_node_values[inol]`.
  - In the lumped inertia loop, `past_dof` of that dof takes the node
    value instead of `old_dof`.
- **Marker:** `step_close` (`top.cc`, task==YES) writes
  `ELEMENT_DOF_INITIAL_APPLIED` for elements with the specific number
  record, in the same way as for `element_dof_initial`. The field
  therefore applies only during the birth step.
- **Validation:**
  - value count other than 1 and less than 1+ndim: `db_error` on the value record;
  - missing value record or unknown dof label: `db_error` on the number record.

## Verification

Professional 02-08-2026 probes, run on the minipc (f5/v6, same inputs as the tests):

| Input | Result | Professional |
|---|---|---|
| 0.1 + 0.8 x | 0.95 | 0.95 |
| 0.8 x, 2 steps | 0.95 (birth step only) | 0.95 |
| `element_dof_initial` 0.4 + specific | 0.7 (precedence) | 0.7 |
| quad4 0.1 + 0.2 x + 0.8 y, nodes 3 / 4 | 0.9625 / 1.0375 | 0.9625 / 1.0375 |
| constant 0.4 | 0.7 (manual) | NaN |
| 2D with one gradient | input error | NaN |

The Professional evaluates the field at the nodes. Evaluating it at the
integration point would give 0.7 for the 0.8 x bar, not 0.9.

Tests: `element_dof_initial_specific1..5` and `_error1`
(`scripts/build_safe.sh`). RED on fe664a3: the keyword did not parse.

## Pending

- `group_dof_initial_specific_number` / `_value` (6.610 / 6.611) are
  not registered.
- Only one dof per element (the record has one number), as in the manual.
