# control_change_dataitem_apply

## Implementación

- **Logic**: in `data()` in `data.cc`, inside the `change_dataitem` block.
  The switch is read once per control record:
  ```
  db( CONTROL_CHANGE_DATAITEM_APPLY, icontrol, &change_dataitem_apply,
    ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  if ( change_dataitem_apply!=-NO ) { ... loop over change_dataitem ... }
  ```
  With `-no` the whole `change_dataitem` loop is skipped.
- **Keyword**: `control_change_dataitem_apply` (data_class CONTROL, type
  INTEGER, data_length 1) registered in `database.cc`. `ival[0]` holds the
  switch.
- **New enum**: `CONTROL_CHANGE_DATAITEM_APPLY` in `tochnog.h` /
  `tochnog-mod.h` (kept in sync), placed between `CONTROL` and
  `CONTROL_CRACK`.

## Diseño / decisiones

- The default is `-YES` (change_dataitem applied), matching the
  Professional behavior ("Default switch is set to -yes").
- The check wraps the existing `change_dataitem` loop, so no other part of
  `data()` is affected.

## Detalles

- `change_dataitem_apply` is initialized to `-YES` and overwritten by
  `GET_IF_EXISTS`.

## Pendiente

- The manual also lists `change_dataitem_geometry` (restrict application to
  elements in a geometry); not implemented here.
