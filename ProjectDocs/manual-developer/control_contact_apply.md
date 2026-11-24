# control_contact_apply (developer)

## Files and functions

- `tochnog.h` / `tochnog-mod.h` — enum `CONTROL_CONTACT_APPLY` (after
  `CONTROL_CHANGE_DATAITEM_APPLY`; headers in sync).
- `database.cc` — registration: `INTEGER`, `data_length = 1`,
  `data_class = CONTACT` (indexed by the control_timestep index).
  No `data_required`: the family convention (contact_apply et al.) has
  no hard combination checks, and requiring the literal `contact`
  keyword rejects valid inputs that only use `contact_geometry` etc.
- `contact.cc` — `parallel_contact()`: the existing `contact_apply`
  gate also reads `ICONTROL` and `CONTROL_CONTACT_APPLY[icontrol]`;
  any `-NO` returns before any contact work. Same precedence shape as
  the `groundflow_consolidation_apply` family (any -no wins).

## Verification

Tests `contact_block` vs `contact_ctrl_apply` (suite 62/62), a clean
A/B on a falling block against a contact face:

- contact active (default): the velocity penalty stops the block,
  disy(node 1) ~ -0.003 after t=0.1.
- `control_contact_apply 0 -no`: free fall, disy = -0.055 exactly
  (implicit integration of g=10, dt=0.01), i.e. numerically identical
  to a model without contact data.
