# control_print_history_relative_time

## Implementación

- **Keyword**: `control_print_history_relative_time` (DOUBLE_PRECISION,
  data_length 1, data_class CONTROL,
  `data_required = CONTROL_PRINT_HISTORY`) registrado en `database.cc`.
- **Lógica**: en `print_history()` (`print_hi.cc`), tras leer
  `TIME_CURRENT` e `ICONTROL`:
  ```
  double time_relative = 0.;
  db( CONTROL_PRINT_HISTORY_RELATIVE_TIME, icontrol, idum, &time_relative,
    ldum, VERSION_NORMAL, GET_IF_EXISTS );
  time_current -= time_relative;
  ```
  `time_relative` INICIALIZADO a 0 (regla GET_IF_EXISTS). El `time_current`
  local se desplaza y se usa en las dos ramas de salida (INTEGER y
  DOUBLE) de `print_history`.
- **Aplicación también en `print_history_smooth()`**: el archivo smooth
  (`smooth<dof><index>.his`) comparte el mismo eje de tiempo relativo
  (mismo record, mismo icontrol) — el manual 6.322 habla de "the history
  files" en general.
- **Afecta solo al tiempo impreso**: `time_current` del db NO se
  modifica (la variable local se resta).

## Enums nuevos

- `CONTROL_PRINT_HISTORY_RELATIVE_TIME` (bloque CONTROL_PRINT, tras
  `CONTROL_PRINT_HISTORY_SMOOTH`) en `tochnog.h` / `tochnog-mod.h`
  (sync).

## Pendiente

- Si el record aparece con `icontrol` distinto del de
  `control_print_history`, no aplica (por diseño: mismo índice, patrón
  `data_required`).
