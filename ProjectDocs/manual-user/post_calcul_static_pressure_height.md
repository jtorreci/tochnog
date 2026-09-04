# post_calcul_static_pressure_height

Reference-height static pressure for the groundflow post-processing split
(manual Professional 6.921): when no groundwater level applies to a node, the
static pressure of `post_calcul -groundflow_pressure -static_pressure` is
determined relative to a reference height instead:
`p_static = rho * g * (height_ref - coord_vertical)`.

Record values: region triples `coord_min coord_max height_ref` along the
VERTICAL coordinate; several regions can be given. A node belongs to the first
region that contains its vertical coordinate; nodes outside every region keep
the phreatic-level behavior (static 0 when no level is defined).

```
post_calcul -groundflow_pressure -total_pressure -groundflow_pressure -static_pressure
post_calcul_static_pressure_height 0. 1. 123.
post_calcul_static_pressure_height_element_group -all
```
Unlocks ground13 (st_pres = -1220 EXACT; further corpus blockers of the
family: bounda_dof -topres conventions, see SEGUIMIENTO).
