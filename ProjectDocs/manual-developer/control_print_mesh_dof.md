# control_print_mesh_dof

## Implementación

- **Alias puro**: `control_print_mesh_dof` es el nombre Professional del
  GNU `print_mesh_dof` (ya implementado en Sprint 9, data.cc). Se
  registra como traducción de prefijo en `db_number()` (`database.cc`):
  ```
  else if ( !strncmp( str, "control_print_mesh_dof", 22 ) ) {
    strcpy( translated, "print_mesh_dof" );
    strncat( translated, &str[22], MCHAR-40 );
    return db_number( translated );
  }
  ```
  La traducción de PREFIJO cubre también las variantes
  `control_print_mesh_dof_geometry` / `_values` (mismo patrón que el
  alias `bounda_print_mesh_dof` de Sprint 9). Sin enum nuevo, sin
  dispatch nuevo.
- **Por qué en db_number**: mismo gotcha que dof_rhside/f178849 — el
  detector de fin-de-valores de records variable-length llama
  `db_number`; traducir solo en el punto del keyword rompería el parseo.
- **Semántica**: el record GNU `PRINT_MESH_DOF` es `no_index = 1` y
  toma la lista de dofs directamente; el manual Professional dice
  `index switch` ("See print_mesh_dof"). El alias mantiene el layout
  GNU (sin índice) — diferencia documentada en el manual-user.

## Detalles

- El dump es one-shot (static `print_mesh_dof_done`), escribe
  `print_mesh_dof.dat` con `inod` (número de nodo == índice interno,
  1-based) + coords + dofs listados.
- LIMITACIÓN descubierta: en un modelo recién creado los dofs listados
  imprimen el valor en la PRIMERA evaluación (todo 0 antes de cargas);
  además, el mapeo `array_member(dof_label, -disy, nuknwn, ...)` depende
  de que el label exista en DOF_LABEL del modelo (comportamiento
  preexistente de Sprint 9, no tocado).

## Pendiente

- La variante `_values` sigue registrada sin uso (Sprint 9).
