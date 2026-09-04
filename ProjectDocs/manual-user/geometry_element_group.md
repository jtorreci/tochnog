# geometry_element_group

Restrict a geometry to nodes of specific element groups (manual Professional
6.524): nodes on the geometry record with the SAME index that are also a node
of elements of one of the listed element groups belong to the geometry; other
nodes do not. Used together with any geometry selector (merge geometries,
force_edge geometries, bounda/delete/excavation geometries, ...).

```
geometry_point 10 0. 0. 1.e20
geometry_element_group 10  1 2
geometry_element_group_method 10  -all
```

Group-attachment rules (see geometry_element_group_method): `-all` requires the
node to be attached to ALL the listed groups (and no other); the default and
`-any`/`-only` accept a node whose elements belong ONLY to the listed groups.

EMPIRICAL NOTE (Professional binary 25-10-2023, measured on merge2 of the
corpus): under every method a node attached to elements OUTSIDE the listed
groups is excluded, even where the manual-literal `-any` would include it. The
GNU reproduces the measured behavior (node of a group-0 element AND of a
group-1 element is not merged when only group 0 is listed).

Consumed inside the node tests of `geometry()` (geometry.cc); the filter is
opt-in (a model without the record pays no scan).
