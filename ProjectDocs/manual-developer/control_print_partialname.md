# control_print_partialname

## Implementación

- **Keyword**: `control_print_partialname` (INTEGER, data_length
  DATA_ITEM_SIZE, `fixed_length = 0`, data_class CONTROL) registrado en
  `database.cc`. Los valores son data item names negativos (e.g.
  `-element`); un valor >= 0 -> `db_error`.
- **Dispatch**: bloque nuevo en `step_close()` (`top.cc`):
  ```
  if ( db_active_index( CONTROL_PRINT_PARTIALNAME, icontrol, VERSION_NORMAL ) ) {
    db( CONTROL_PRINT_PARTIALNAME, icontrol, ival, ddum, nval, VERSION_NORMAL, GET );
    print_partialname( icontrol, VERSION_NORMAL, ival, nval );
  }
  ```
- **print_partialname()** (`print_db.cc`): bucle sobre todos los data
  items (`idat` 0..MDAT-1); para cada PREFIJO (`ival[iprefix]` negativo
  -> `db_name(labs(...))` = la cadena del prefijo, e.g. "element") se
  testea `!strncmp( db_name(idat), prefix, strlen(prefix) )` — match de
  PREFIJO, como el manual ("all records starting with data_item_name").
  Para cada item que casa se imprimen TODOS los records activos a
  stdout (mismo formato que `control_print`: `nombre indice valores`).
  Se aplica `control_print_filter` (mismo icontrol).
- **Prefijo vs helper**: el helper GNU `db_partialname()` usa `strstr`
  (substring en cualquier posición); el manual exige PREFIX ("starting
  with") — el nuevo camino usa `strncmp` (strstr habría impreso también
  `group_element*` con el prefijo `element`).
- **Salida a stdout**: igual que `control_print` (no a .dbs).

## Enums nuevos

- `CONTROL_PRINT_PARTIALNAME` (bloque CONTROL_PRINT, junto a
  `CONTROL_PRINT_NUMBER_ITERATIONS`) en `tochnog.h` / `tochnog-mod.h`
  (sync).

## Pendiente

- El `control_print_filter` se aplica por data item (data_number =
  idat); el filtro no distingue entre los distintos prefijos del record.
- Nada más: el formato de impresión es el del print genérico.
