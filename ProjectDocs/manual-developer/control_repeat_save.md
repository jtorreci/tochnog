# control_repeat_save

## Implementación

- **Record**: `control_repeat_save` (INTEGER, variable length — one
  triplet `(data_item_name, data_item_index, data_item_number)` per
  item to save, `fixed_length = 0`, data_class CONTROL,
  `data_required = CONTROL_REPEAT`) registered in `database.cc`.
- **Lógica**: en `repeat.cc`, dentro de la rama `CONTROL_REPEAT` de
  `repeat()` (llamada desde `top.cc` tras cada `step_close`): cuando el
  contador del repeat es `> 0` (se va a saltar de vuelta a
  `start_control`), se llama a la función estática
  `control_repeat_save_data( icontrol )`:
  1. Lee el record `CONTROL_REPEAT_SAVE` en `icontrol` (el mismo índice
     del `control_repeat`, como en el ejemplo de Monte Carlo del manual
     Professional 2.6.4) y calcula `ndata = length/3`.
  2. Para cada triplete lee el valor actual del data item (mismo patrón
     de `control_repeat_until_item`: `data_item_number < 0` → dof label
     vía `array_member(dof_label, ...)` con la corrección
     `number /= nder` si el record tiene longitud `npuknwn`; si no,
     número positivo directo). Bounds estrictos `0 <= number < length`
     (el código GNU del until_item usa `<=`, permitiendo un valor fuera
     de rango; aquí se usa `<` por seguridad).
  3. Escribe `REPEAT_SAVE_RESULT[isave] = [v_0, v_1, ..., v_{ndata-1}]`
     con `isave` = número de índices ACTIVOS ya escritos.
- **GOTCHA (db_max_index)**: `db_max_index` devuelve el máximo
  ALOCADO, que incluye el margen heurístico de `db_allocate`
  (`increase = index/10`, mínimo 1). Usarlo como "siguiente índice
  libre" salta huecos (saves en 0, 2, 3, ...) y el cálculo posterior
  aborta al leer un índice inactivo (`db_error`). El índice libre se
  obtiene CONTANDO los índices activos (`db_active_index` de 0 a
  `db_max_index`).
- **Timing**: el save captura el estado AL FINAL de la pasada que
  acaba de completarse (el `control_repeat_save` está en el mismo
  índice que `control_repeat`, que se procesa después del bloque
  repetido). Con contador inicial N hay N saltos → N saves (pasadas
  0..N-1); la pasada N (contador 0) no guarda.

## Layout de REPEAT_SAVE_RESULT

- `repeat_save_result` (DOUBLE, `fixed_length = 0`,
  `data_length = MCALCUL` — cap de 20 data items, ver PENDIENTE).
- Un índice por repeat: `repeat_save_result[irepeat]` = valores de
  todos los items en orden de tripletas en el momento del salto
  (repeat 0 → índice 0, repeat 1 → índice 1, ...). Los repeats
  subsiguientes escriben en índices subsiguientes (manual 6.348).
- La escritura va a VERSION_NORMAL directamente (fuera de loops
  paralelos; `repeat()` corre en el flujo principal de `top.cc` tras
  el `db_version_copy(NEW → NORMAL)` del cierre de paso, así que no la
  pisa el copy de versión).

## PENDIENTE

- `data_length = MCALCUL (20)`: si el usuario pide más de 20 data
  items el PUT aborta con "Length too small of repeat_save_result".
- No se soporta con `control_repeat_until_*` (el manual la liga solo a
  `control_repeat`).
- Repeats anidados: cada paso por el repeat interior añade saves
  (contador global); el `control_repeat_save_calculate` se re-ejecuta
  al completar cada repeat interior y sobreescribe
  `REPEAT_CALCULATE_RESULT` (último write gana) — comportamiento
  documentado, no validado.
