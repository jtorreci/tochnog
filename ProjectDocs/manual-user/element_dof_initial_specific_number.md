# element_dof_initial_specific_number

## Description

`element_dof_initial_specific_number` together with
`element_dof_initial_specific_value` gives an element an initial "past"
value for ONE dof (manual Professional 6.423 / 6.424). It works like
`element_dof_initial`: when the element comes the first time to live it
assumes it had this value in the past, so the transient inertia of its
birth step (heat capacity, mass, water storage) integrates from it.
After the birth step the element continues from its regular dofs.

Unlike `element_dof_initial`, the value can vary linearly in space:

```
value(x) = value_0 + value_grad_x * x + value_grad_y * y + value_grad_z * z
```

It is evaluated at the coordinates of each element node.

## Usage

```
element_dof_initial_specific_number <element> <dof>
element_dof_initial_specific_value  <element> value_0 [value_grad_x value_grad_y value_grad_z]
```

- `dof`: the dof label, for example `-temp`, `-velx` or `-sigxx`.
- `value_0`: the value at x = y = z = 0.
- `value_grad_*`: the gradients. Give one per space dimension (only
  `value_grad_x` in 1D, `value_grad_x value_grad_y` in 2D), or none at
  all for a constant value `value_0`. Any other count is an input error.
- When the element also has `element_dof_initial`, that record is used
  and the specific records are ignored (as in the Professional).
- A specific number without a specific value, or a dof that is not
  active, is an input error.

## Example

Bar of `phase1.dat` (left node prescribed 1, lumped inertia 1 per node):

```
element_dof_initial_specific_number 1 -temp
element_dof_initial_specific_value 1 0.1 0.8
```

The right node starts from the past value 0.9, and the dynamic
equilibrium `(1 - T) = (T - 0.9)` gives `node_dof 2 = 0.95`
(`validation-suite/test-2026/element_dof_initial_specific1.dat`).

## Differences with the Professional version

- `value_0` alone (a constant field, as the manual allows) makes the
  Professional 02-08-2026 stop with a NaN; GNU uses the constant value
  (DIVERGENCES DOF-INITIAL-SPECIFIC-CONSTANT).
- A gradient count other than 0 or ndim also gives a NaN in the
  Professional; GNU stops with an input error.
