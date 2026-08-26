# control_print_database_method

## Implementación

- **Keyword**: `control_print_database_method` (INTEGER, data_length 1,
  data_class CONTROL) registered in `database.cc`. `ival[0]` holds the
  method switch: `-all`, `-size_tot` o `-size_tot_large`.
- **Dispatch**: bloque nuevo en `step_close()` (`top.cc`), junto al de
  `CONTROL_PRINT_DATABASE`:
  ```
  if ( db_active_index( CONTROL_PRINT_DATABASE_METHOD, icontrol, ... ) ) {
    db( CONTROL_PRINT_DATABASE_METHOD, icontrol, ival, ddum, ldum, ... );
    if      ( ival[0]==-ALL ) print_database( icontrol, VERSION_NORMAL, -EVERYTHING );
    else if ( ival[0]==-SIZETOT ) print_database( icontrol, VERSION_NORMAL, -SIZETOT );
    else if ( ival[0]==-SIZE_TOT_LARGE ) print_database( icontrol, VERSION_NORMAL, -SIZE_TOT_LARGE );
    else db_error( CONTROL_PRINT_DATABASE_METHOD, icontrol );
  }
  ```
  El `-all` del manual se mapea a `-EVERYTHING` (imprime TODOS los
  records base, no solo los `db_external` como hace `-YES`).
- **print_database()** (`print_db.cc`): rama `-SIZE_TOT_LARGE` nueva —
  mismo cálculo de `size = max * db_data_length * sizeof(type)` que
  `-SIZETOT`, pero solo imprime si `size > 1024*1024` (1 Mb). El
  filename (`.dbs`) y la línea "Total size is" se comparten.
- **Tamaño de la matriz del sistema**: `solve_nlocal` (global definido
  en `initia.cc`, `extern` ya declarado en `tochnog.h`) — el número de
  ecuaciones del último `solve()`. Es la misma métrica que el reporte
  del band solver en `so.cc` ("number of equations/size of matrix").
  DECISIÓN: se imprime el número de ecuaciones (no bytes) — es la única
  métrica de tamaño disponible de forma fiable fuera del solver; el
  almacenamiento real depende del solver (band/sparse).

## Enums nuevos

- `CONTROL_PRINT_DATABASE_METHOD` (bloque CONTROL_PRINT) y
  `SIZE_TOT_LARGE` (junto a `SIZETOT`) en `tochnog.h` / `tochnog-mod.h`
  (sync), con `name[SIZE_TOT_LARGE] = "size_tot_large"` y
  `name[SIZETOT] = "sizetot"`.
- **GOTCHA de nombre**: el manual Professional escribe `-size_tot`
  (con guion bajo) pero el GNU registra `sizetot` (sin guion). Se añadió
  la traducción `db_number("size_tot") -> SIZETOT` (patrón Sprint 9:
  la traducción en `db_number`, no solo en el punto del keyword — el
  detector de fin-de-valores de records variable-length también llama
  `db_number`).

## Pendiente

- La línea del sistema matrix usa `solve_nlocal` del último solve; en
  un run sin solve (nuknwn == 0) imprime 0.
