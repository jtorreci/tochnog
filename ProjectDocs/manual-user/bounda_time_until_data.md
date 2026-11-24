# bounda_time_until_data and bounda_time_until_value_minimum

Manual Professional sections 6.40 and 6.41.

## What it does

`bounda_time_until_data` monitors a data item (for example the reaction
force of a node) and reduces the load prescribed by the `bounda_time`
record with the same index when the monitored value falls. This is used
to slowly release a prescribed velocity/displacement so that the
structure reaches a minimum force response.

A typical application (from the Professional manual 6.40): prescribe a
velocity on a structure, monitor its force response, and reduce the
prescribed velocity when the force response drops below a threshold:

```
bounda_dof 210 ... -vely                 (apply velocity on the structure)
bounda_time 210 ...
bounda_time_until_data 210 -post_node_result 10 -vely   (monitor force response)
bounda_time_until_value_minimum 210 0.1 0.3             (reduce when below 30%,
                                                         zero when 10%)
```

## Input syntax

```
bounda_time_until_data index data_item_name data_item_index data_item_number
bounda_time_until_value_minimum index wanted start
```

- `index` — the index of the `bounda_time` record that is reduced.
- `data_item_name` — the data item to monitor, e.g. `-post_node_result`
  (a `post_node` record that sums a nodal quantity), `-node_rhside`, etc.
- `data_item_index` — index of the monitored data item.
- `data_item_number` — the component to monitor, e.g. `-velx`.
- `wanted` — monitored value at which the applied load becomes 0.
- `start` — monitored value below which the load starts to be reduced.

## Semantics (verified against the Professional binary 25-10-2023)

The load of the `bounda_time` record with the same index is multiplied
by a quadratic reduction factor:

```
factor = clamp( ((monitor/first - wanted) / (start - wanted))^2 , 0, 1 )
```

where:

- `monitor` is the value of the data item of the PREVIOUS time step;
- `first` is the value of the monitor when the mechanism starts (written
  to the database as `bounda_time_until_first`).

So while `monitor >= start*first` the full load is applied, and the load
decreases quadratically to zero as `monitor` approaches `wanted*first`.
The applied factor of the last step is written to the database as
`bounda_time_until_used` (same output records as the Professional .dbs).

## Example (corpus test until1.dat, passes rc=0)

```
bounda_dof   10  1 -velx
bounda_dof   20  2 -velx
bounda_time  20  0.0 1.0 1.0 1.0
bounda_dof   30  2 -velx
bounda_time  30  1.001 -1.0 5.0 -1.0
bounda_time_until_data 30  -post_node_result 10 -velx
bounda_time_until_value_minimum 30  0.0 1.e-2
...
target_item  0   -node_rhside 1 -velx
target_value 0   0.0 1.e-3
```

The node 2 is first pulled at velocity +1.0 (bounda_time 20) and then
pushed at -1.0 (bounda_time 30). The `until` mechanism monitors the
reaction of node 1 and reduces the -1.0 load as the reaction falls from
`start*first` to `wanted*first`, so that the final reaction is ~0.

## Related records

- `bounda_time_until_value` (3 values, corpus 02-08-2026 variant) is
  registered as a known keyword but its consumption is PENDING (the
  25-10-2023 binary rejects the record, so it could not be verified).
- `bounda_time_until_force` (GNU, 2026-08-05) covers the reaction-force
  case with a simpler linear reduction.
