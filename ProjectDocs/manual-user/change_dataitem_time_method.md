# change_dataitem_time_method

## Description

Requires that the cosinus, sinus or tangent of a data value will be
changed, instead of the data value directly itself. The time table of
the `change_dataitem_time` record (with the same index) contains
cosinus, sinus or tangent VALUES; the value actually stored in the
changed data item is the inverse angle: `acos(val)`, `asin(val)` or
`atan(val)` respectively (manual Professional 6.51).

This is typically convenient for geotechnical safety factor
calculations where you want a Mohr-Coulomb law's cohesion and the
tangent of the friction angle to be decreased at the same ratio in
time (phi-c reduction): the tables give the tangent directly.

## Usage

```
change_dataitem_time_method <index> <-cosinus|-sinus|-tangent>
```

## Parameters

| Parameter | Meaning                                                          |
|-----------|------------------------------------------------------------------|
| `index`   | Index of the `change_dataitem`/`change_dataitem_time` records.   |
| method    | `-cosinus` store `acos(val)`, `-sinus` store `asin(val)`, `-tangent` store `atan(val)`. Default: store `val` unchanged. |

## Example

```
group_materi_plasti_mohr_coul_direct 10  0.5236 50. 0.

(tangent of friction angle reduction)
change_dataitem              10  -group_materi_plasti_mohr_coul_direct 10 0 -use
change_dataitem_time         10  0.0 0.577   100. 0.40
change_dataitem_time_method  10  -tangent

(cohesion reduction)
change_dataitem              20  -group_materi_plasti_mohr_coul_direct 10 1 -use
change_dataitem_time         20  0.0 50.  100. 35.
```

With the `-tangent` method the first table's values (0.577 = tan 30 deg,
0.40 = tan ~21.8 deg) are stored as the ANGLES phi = atan(val); the
cohesion table is stored directly. Both are reduced at consistent ratios.
