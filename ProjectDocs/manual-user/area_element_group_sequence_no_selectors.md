# area_element_group_sequence (without selectors)

## Description

`area_element_group_sequence` (manual Professional 6.9) can be used as a
completely separate option WITHOUT the `_geometry`/`_element` selectors:
then the elements of the **previous group number in the sequence**
(`group_(i-1)`) get the new group number `group_i` at `time_i`. The
previous group selects the elements (the element must currently have the
group of the previous time window).

The GNU previously REQUIRED one of the selectors and errored out
("AREA_ELEMENT_GROUP_SEQUENCE_GEOMETRY or AREA_ELEMENT_GROUP_SEQUENCE
should be specified"); dam_building uses the sequence with only
`area_element_group_sequence_element_group` + `_time`.

## Usage

```
area_element_group_sequence_element_group 0  7   107
area_element_group_sequence_time          0  start_time start_layer1
area_element_group_sequence_element_group 1  8   108
area_element_group_sequence_time          1  start_time start_layer1
...
```

## Notes

- When neither `_geometry` nor `_element` is used, the element must have
  the group of the previous time window (`group_(i-1)` at `time_(i-1)`),
  otherwise it is left unchanged.
- dam_building (layered dam construction, ~7000 element groups) parses
  and runs with this mode; the full model exceeds the corpus time
  budget (~10 minutes vs 45 s).
