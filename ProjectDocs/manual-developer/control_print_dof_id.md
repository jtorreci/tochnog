# control_print_dof_id

## Implementación

- **Keyword**: `control_print_dof_id` (INTEGER, data_length 1,
  data_class CONTROL, `data_required = CONTROL_PRINT_DOF`) registrado en
  `database.cc` — el manual lo liga a `control_print_dof` ("works in
  combination with the control_print_dof record").
- **Lógica**: dentro de `print_dof()` (`print_hi.cc`), antes del
  `db_version_copy(VERSION_NORMAL, VERSION_PRINT)`:
  1. `db( CONTROL_PRINT_DOF_ID, icontrol, &dof_id, ... GET_IF_EXISTS )`
     con `dof_id` INICIALIZADO a `-YES` (default del manual 6.270 y
     regla GET_IF_EXISTS: inicializar el buffer — gotcha del NaN).
  2. Si `-yes`, se captura el MAPEO posición-compacta -> número de nodo
     ORIGINAL: `node_number_of_position[pos++] = inod` recorriendo
     `NODE` en VERSION_NORMAL (los nodos activos en orden de índice).
  3. En el bucle de salida, tras `x y z <dof>` se añade
     `node_number_of_position[inod]` (inod = posición compacta).
- **Por qué el mapeo**: `renumbering(VERSION_PRINT, ...)` compacta los
  índices de nodo a 0..n-1; el número de nodo ORIGINAL (el 'identity'
  del manual) se pierde. El mapeo es exacto porque la renumeración
  asigna los números nuevos en el MISMO orden de índice creciente.
  `print_mesh_dof` (data.cc) NO renumera y por eso imprime el índice
  crudo directamente (que ES el número de nodo).
- **Convención de numeración (descubrimiento)**: en tochnog el índice
  interno de `node N` es N (1-based; el índice 0 queda sin usar). Tras
  renumbering, la posición compacta p corresponde al p-ésimo nodo
  activo.

## Enums nuevos

- `CONTROL_PRINT_DOF_ID` (bloque CONTROL_PRINT, tras
  `CONTROL_PRINT_DOF`) en `tochnog.h` / `tochnog-mod.h` (sync).

## Detalles

- `dof_id` inválido (!= -YES/-NO) -> `db_error`.
- El archivo `coord.<index>` no cambia con el switch.

## Pendiente

- Con `control_print_dof` vía `-separate_sequential`, el contador
  estático `dof_seq` y el mapeo se calculan por llamada — correcto pero
  no validado en secuencia larga.
