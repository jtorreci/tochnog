# change_dataitem_geometry

## Description

For element group data (`group_*`) you can restrict the application of
the `change_dataitem` (with the same index) to only those elements which
are part of the geometry specified by `geometry_entity_name
geometry_entity_index` (manual Professional 6.48).

The restriction is materialized by splitting the element group: the
first time the change fires, a clone group is created with identical
copies of all `group_*` records, and the elements of the original group
that are FULLY inside the geometry move to the clone. The value change
is applied to the clone only: elements inside the geometry see the new
value, elements outside keep the original value (and the original group
record is left untouched). Element internal state (stresses, history
variables) is not affected by the move.

Typical use: reduce strength parameters (phi-c reduction) only in a
zone, e.g. under a footing, keeping the rest of the model at full
strength.

## Usage

```
change_dataitem_geometry <index> -<geometry_entity_name> <geometry_entity_index>
```

## Parameters

| Parameter | Meaning                                                            |
|-----------|--------------------------------------------------------------------|
| `index`   | Index of the `change_dataitem` record the restriction applies to.  |
| geometry  | Geometry entity (e.g. `-geometry_brick`, `-geometry_circle`) and its index. `geometry_brick` (a box) is the natural region selector. |

Notes:
- Only valid for `group_*` targets (using it with another data item is an
  input error).
- Elements must be FULLY inside the geometry (all nodes) to move to the
  clone; partially covered elements keep the original parameters.

## Example

```
geometry_brick  1  0.5 0.5 0.  1.5 1.5 3.0  0.01

change_dataitem           10  -group_materi_plasti_mohr_coul_direct 0 1 -use
change_dataitem_time      10  0.0 1.0  0.05 0.0  100. 0.0
change_dataitem_geometry  10  -geometry_brick 1
```

The cohesion (value 1) of `group_materi_plasti_mohr_coul_direct` drops
to 0 from t=0.05, but ONLY for the elements inside the brick: the clone
group gets c=0, the original group keeps c=1.0.
