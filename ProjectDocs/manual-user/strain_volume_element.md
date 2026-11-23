# strain_volume_element (and the strain_volume_* family)

## Description

Prescribed volume strain of a set of elements (manual 6.1062-6.1066), for
example grouting:

```
strain_volume_relative_time index time_0 percentage_0 time_1 percentage_1 ...
strain_volume_absolute_time index time_0 volume_increase_0 ...
strain_volume_element       index element_0 element_1 ...
strain_volume_element_group index element_group_0 element_group_1 ...
strain_volume_geometry      index geometry_item_name geometry_item_index
```

## Effect on the calculation

- A third of the volume strain eps_v is an extra normal strain in x, y and
  z (also z in plane strain) that the stress does not see:
  stress = D (eps - eps_v/3 I). The actual deformation depends on the
  boundary conditions: a fully fixed hex8 with 3 % gets
  sxx = syy = szz = -E/(1-2nu)*0.01; a free 3D element grows by eps_v.
- eps_v is the value of the time table at the time point (not summed over
  the steps); outside the table it is 0.
- Relative table: percentage of the initial volume (100 = double volume).
- Absolute table: volume increase divided by the initial volume of the
  whole set (`post_strain_volume_initial`), the same strain for every
  element of the set.
- The set: the elements that satisfy every selector you give, as in the
  Professional 02-08-2026:
  - listed in `strain_volume_element`;
  - in one of the groups of `strain_volume_element_group`;
  - with all nodes inside the geometry of `strain_volume_geometry`.

  For example, element 1 with groups 0 and 1 selects element 1 only. When no
  element satisfies all the selectors, the set is empty and every result is
  0. Without any of the three selectors, every isoparametric element of a
  `-materi` group is in the set.
- Exactly one of the two tables per index; both is an input error.
- An element in several indices takes the last index.
- Small deformation, isoparametric elements only.

Output (same index): `post_strain_volume_initial` (initial volume of the
set), `post_strain_volume_absolute` (summed actual volume increase, the
integral of the strain trace) and `post_strain_volume_relative`
(100 * absolute / initial).

## Example

```
strain_volume_absolute_time 10  0. 0. 3. 0.1234 1.e10 0.1234
strain_volume_element       10  1
```

Self-checking inputs in `validation-suite/test-2026`: `strain_volume1`
(Professional example), `strain_volume_fixed1`, `strain_volume_plane1`,
`strain_volume_all1`, `strain_volume_element_list1`, `strain_volume_group1`,
`strain_volume_geometry1`, `strain_volume_index1`, `strain_volume_index2`,
`strain_volume_error1`.

## Differences with the Professional version

Measured on the Professional 02-08-2026; every test above matches it to
1e-8 except two cases where it contradicts the manual:

- `strain_volume_element` with several elements: the Professional uses only
  the first one (DIVERGENCES STRAIN-VOLUME-ELEMENT-LIST).
- Both tables in one index: the Professional silently uses the absolute
  one (DIVERGENCES STRAIN-VOLUME-BOTH-TABLES).
