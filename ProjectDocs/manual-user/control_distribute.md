# control_distribute

## Description

Apply a random number, based on a `-normal` or `-lognormal` distribution,
to the selected data item records (manual Professional 6.127):

```
control_distribute <index> <distribution_type> <data_item_name> <data_item_index> <data_item_number>
```

- `data_item_index` selects the record index, or `-all` for all existing
  indices.
- `data_item_number` selects the value inside the record (0 for the first
  value, 1 for the second, ...); for dof-based records a dof label like
  `-temp` can be used.
- `control_distribute_parameters` (same index) gives the mean value and
  standard deviation of the distributed value itself.
- `data_item_name` can be `group_*` or `node_*`. For a `group_*` item the
  record itself is NOT changed: each element using the record gets its own
  drawn value (per-element data). For node items the drawn value REPLACES
  the selected number of the record (e.g. nodal initial temperatures, or
  node coordinates).

Optional companions (same index):

- `control_distribute_seed` — makes the random sequence reproducible:
  equal seeds give identical fields.
- `control_distribute_correlation_length` — spatially correlated field
  (1..ndim values, one per direction). A length larger than 1.e12 gives a
  CONSTANT field (all components the same value). Data further apart than
  the correlation distance are not correlated; the distance defaults to
  4 times the correlation length (`control_distribute_correlation_distance`
  overrides it).
- `control_distribute_minimum_maximum` — cut-off bounds for the drawn
  values (e.g. keep a void ratio inside the range required by a
  hypoplasticity law).

The distribution fires at every timestep of the `control_timestep` block
with the same index (the field is redrawn per step); it should be placed
BEFORE any `control_reset_dof` records (lower index).

### Differences with the GNU legacy layout

The GNU keyword of the same name uses a different record layout
(triplets with `-uniform`/`-normal` plus a `control_distribute_values`
record of deltas, added to the current values). Both layouts are
supported: a 4-value record is the Professional layout, triplets are the
GNU one. New models should use the Professional layout.

## Example

```
control_distribute             10  -normal -group_materi_elasti_young 0 0
control_distribute_parameters  10  1000. 50.
control_distribute_seed        10  7.
```

Every element of group 0 gets Young = a draw from N(1000, 50), identical
across reruns (seed 7).

See also the tests `cdist_normal`, `cdist_corr`, `cdist_clamp`
(validation-suite).
