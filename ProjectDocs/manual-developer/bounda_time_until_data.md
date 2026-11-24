# bounda_time_until_data and bounda_time_until_value_minimum — developer notes

Manual Professional 6.40/6.41.

## Where implemented

- `bounda.cc`, routine `bounda()`, in the per-`iboun` block that reads
  the `bounda_time` record (after the `BOUNDA_TIME`/`BOUNDA_TIME_FILE`
  dispatch, before `BOUNDA_TIME_ON_OFF`).

### Reading the monitor (per iboun, once per step)

The record `bounda_time_until_data iboun name index number` is read with
`db( BOUNDA_TIME_UNTIL_DATA, iboun, ... )`:

- `until_data[0]` — data item name (negative enum, e.g. `-post_node_result`)
- `until_data[1]` — data item index
- `until_data[2]` — data item number (dof label, e.g. `-velx`)

The component slot is resolved exactly like `target_item` in
`miscel.cc` `exit_tn()`: `array_member(dof_label, number, nuknwn, indx)`
and, when the data item length equals `npuknwn`, `indx /= nder`.

The monitor value is the data item value of the PREVIOUS time step:
`post_node_result` is written by `post()` at the end of the step, and
`bounda()` runs at the start of the next step, so the value available
here is the one of the previous step (verified: matches the
Professional series exactly).

### The reduction factor

Computed on the first step in which the `bounda_time` record is active
(`found==1` in the time loop):

- `first` (written once as `bounda_time_until_first`): the current
  monitor value.
- `factor = clamp(((monitor/first - wanted)/(start - wanted))^2, 0, 1)`.
- `load *= factor`.
- `bounda_time_until_used` is written with the applied factor (output
  parity with the Professional .dbs).

The factor is applied inside the `for (inc...)` loop right after the
load is computed and `found` is set, and before the node loop.

## Verified against the Professional

Two runs of `until1.dat` with the user-supplied Professional binary
25-10-2023 (E=1 and E=2, i.e. young modulus 1.0 and 2.0) were compared
series by series (time_current, node_dof 2, post_node_result):

- E=1: first=1.0; factor at t=1.992 = 0.81 = (9e-3/1e-2)^2 — exact.
- E=2: first=2.0; factor at t=1.992 = 0.81 = (1.8e-2/(1e-2*2))^2 — exact.
- Full series match to 1e-6.

## Records added in database.cc

- `BOUNDA_TIME_UNTIL_DATA` — INTEGER, length 3 (name, index, number),
  fixed length, `data_required = BOUNDA_TIME`.
- `BOUNDA_TIME_UNTIL_VALUE_MINIMUM` — DOUBLE_PRECISION, length 2
  (wanted, start), fixed length, `data_required = BOUNDA_TIME_UNTIL_DATA`.
- `BOUNDA_TIME_UNTIL_FIRST` / `BOUNDA_TIME_UNTIL_USED` — DOUBLE_PRECISION,
  length 1; output records written by bounda() (Professional .dbs parity).
- `BOUNDA_TIME_UNTIL_VALUE` — DOUBLE_PRECISION, length 3; REGISTERED as a
  known keyword for the 02-08-2026 corpus (validation_14_mesh) but the
  consumption is PENDING (the 25-10-2023 Professional binary rejects the
  record, so its semantics could not be verified).

## Keywords registered for the corpus parser only

The GNU parser ends a variable-length record (`fixed_length=0`) when it
meets a KNOWN keyword. The corpus tests put the following families right
after `bounda_time`/`bounda_dof` records; without their registration the
parser attributed the error to bounda_time ("Problem reading :
bounda_time / I don't know what to do with : <keyword>"). They are
registered in database.cc with the correct type/length but WITHOUT
consumption (see SEGUIMIENTO-CONVERGENCIA.md):

- `mpc_node_number`, `mpc_node_factor`, `mpc_geometry`,
  `mpc_geometry_dof`, `mpc_geometry_switch`, `mpc_linear_quadratic`,
  `mpc_element_group` (family mpc_*, manual 6.856-6.875).
- `control_mesh_truss_distribute_mpc`, `..._exact` (manual 6.245/6.250).
- `post_calcul_length` (written by the Professional in the .dbs).
- `strain_volume_absolute_time`, `strain_volume_element`,
  `post_strain_volume_absolute`, `post_strain_volume_relative`
  (family strain_volume, manual 6.96x).
- `bounda_used` (output record of the Professional .dbs).

## Pending

- `bounda_time_until_value` consumption (3-value variant).
- The `mpc_*`, `control_mesh_truss_distribute_mpc`, `strain_volume_*`
  families: registered, consumption pending.
- force16/force17/mpc8/mpc9/post7 segfault in `post_element_force`
  (pre-existing, not related to this lot; the corpus counts them as
  RUNFAIL as before).
