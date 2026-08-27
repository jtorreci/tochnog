# control_print_dof_line_element_group

## Implementación

- Stored as a CONTROL INTEGER record with `data_length` DATA_ITEM_SIZE,
  `fixed_length 0` (variable-length group list). Read in
  `print_dof_line_point()` of `print_dl.cc` into `group_list[ngroups]`.
- Applied in `dof_line_interpolate()` BEFORE the point-in-element test:
  ```
  if ( ngroups>0 ) {
    db( ELEMENT_GROUP, element, &element_group, ddum, ldum,
      VERSION_NORMAL, GET_IF_EXISTS );
    in_group = 0;
    for ( igroup=0; igroup<ngroups; igroup++ )
      if ( group_list[igroup]==element_group ) in_group = 1;
    if ( !in_group ) continue;
  }
  ```

## Diseño / decisiones

- An element without an `ELEMENT_GROUP` record defaults to group 0
  (`GET_IF_EXISTS` leaves `element_group` at its initialised 0).
- A point whose containing elements are all outside the listed groups is
  NOT found: it is omitted from the files (documented decision, same
  behaviour as an out-of-mesh point).
- No list -> all groups are searched (the filter is only applied when the
  record exists).

## Detalles

- The filter runs before `point_el()`, so excluded elements never
  participate in the "first element wins" race.

## Pendiente

- None.
