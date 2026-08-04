# control_mesh_switch

## Archivos y funciones

- `mesh.cc` — `void mesh_switch( long int control_mesh_switch[], long int length )` (línea 131). Implementación completa de la permutación de ejes.
- `top.cc` — llamado en `step_start` (líneas 486–491), justo después del bloque de `adjust_geom` (líneas 480–484).
- `tochnog.h` — enum `CONTROL_MESH_SWITCH` (línea 251).
- `database.cc` — registro del keyword (líneas 718–722).

## Detalles de implementación

- Parseo de la permutación: para cada dimensión `idim < ndim && idim < length`,
  `-X/-Y/-Z` se traduce a un índice 0/1/2 guardado en `order[idim]`
  (`order[idim]` = eje OLD que alimenta al eje NEW `idim`). `ndim` viene del modelo.
- Validación (errores fatales con `exit(TN_EXIT_STATUS)`):
  - los tokens deben ser `-x`, `-y` o `-z`;
  - deben cubrir TODOS los ejes de `ndim` (`order[idim] < 0` → error);
  - deben formar una permutación (repetir eje → error, array `has[]`).
- Aplicación: recorre todos los nodos activos (`db_max_index(NODE,...)` +
  `db_active_index`), lee `NODE`, permuta con `new_coords[idim] = coords[order[idim]]`
  y hace `PUT`. Repite la operación sobre `NODE_START_REFINED` si el nodo está activo allí.
- Cierra con `mesh_has_changed( VERSION_NORMAL )` para forzar el refresco de la malla.
- `database.cc`: `type = INTEGER`, `data_length = DATA_ITEM_SIZE`,
  `fixed_length = 0` (longitud variable = número de ejes dados), `data_class = CONTROL`.

## Dependencias externas

Ninguna. Usa solo la API interna del database (`db`, `db_active_index`,
`db_max_index`) y `mesh_has_changed`. `-X/-Y/-Z` y `MDIM` vienen de las constantes
del programa.

## Parámetros hardcodeados / refactorizaciones pendientes

- Solo usa los primeros `ndim` valores; los sobrantes del `DATA_ITEM_SIZE` se ignoran
  silenciosamente. Se podría validar que `length == ndim`.
- No soporta signos (reflexión/rotación con cambio de signo, p. ej. `-y -x` negado);
  solo permutaciones puras. Una rotación de 90° requiere extensión con signo.
- El código de permutación de `NODE` y `NODE_START_REFINED` está duplicado; se
  podría extraer un helper `permute_node_coords(inod, order, item)`.
