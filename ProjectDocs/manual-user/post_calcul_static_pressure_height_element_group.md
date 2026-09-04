# post_calcul_static_pressure_height_element_group

Restrict the regions of post_calcul_static_pressure_height to element groups
(manual Professional 6.922): the i-th value is the element group for which the
i-th region is valid; the special value `-all` makes the region valid for all
element groups. A node is served by the region whose group matches one of the
node's element groups.
