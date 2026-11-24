# control_print_element_method

## Implementación

- **Keyword**: `control_print_element_method` (INTEGER, data_length 1,
  data_class CONTROL, `data_required = CONTROL_PRINT_ELEMENT`)
  registrado en `database.cc`. `ival[0]` es el método (`-middle` /
  `-node`).
- **Lógica**: dentro de `print_element()` (`print_el.cc`):
  `db( CONTROL_PRINT_ELEMENT_METHOD, icontrol, &method, ... GET_IF_EXISTS )`
  con `method` INICIALIZADO a `-MIDDLE` (default del manual 6.286).
  `method` inválido -> `db_error`.
  - `-middle`: una línea por elemento — coordenada media
    (`suma(coords de los nodos del elemento)/nnol` por dimensión) +
    `dval[ival]` (el valor del data item; para los records
    element-data es ya el valor promediado del elemento).
  - `-node`: el comportamiento GNU previo — por cada nodo del elemento
    (con la lógica `ok` de truss/beam: el truss imprime ambos nodos; el
    beam imprime el nodo 0 para `ival < nval/2` (q) y el nodo 1 para
    `ival >= nval/2` (m)), coordenadas nodales + `dval[ival]`.
- **CAMBIO de comportamiento del default**: antes el GNU imprimía
  siempre el formato por-nodo; ahora el default es `-middle` (manual).
  Ningún test legacy usa `control_print_element` (verificado con grep en
  validation-suite y sfnet) — el cambio es seguro.
- **Naming GNU vs manual**: los archivos siguen la convención GNU
  `element_truss_force_<ival>.<index>` /
  `element_beam_force_moment_<ival>.<index>`; el manual Professional
  usa `element_truss_force_n.<index>` y
  `element_beam_force_moment_q/m.<index>`. Diferencia documentada, sin
  renombrar (fuera de alcance).

## Enums nuevos

- `CONTROL_PRINT_ELEMENT_METHOD` (bloque CONTROL_PRINT, tras
  `CONTROL_PRINT_ELEMENT`) y `MIDDLE` (switch, junto a METHOD2/MINIMAL)
  en `tochnog.h` / `tochnog-mod.h` (sync), con
  `name[MIDDLE] = "middle"` — el parser convierte `-middle` via
  `db_number("middle")`. `-node` ya existía (`NODE`).

## Detalles

- GOTCHA de validación: la fuerza de un `element_truss_force` inicial es
  un self-stress que SE RELAJA con la deformación (new_force =
  old_force + k*dL; en equilibrio final queda 0). El test elmethod
  discrimina por ESTRUCTURA (conteo de líneas 2 vs 4 y la relación
  exacta media = promedio de coords nodales), no por el valor impreso.
