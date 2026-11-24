# bounda_time_on_off

## Archivos y funciones

- `bounda.cc` — bloque dentro de `bounda()` (líneas 125–134), tras el encadenado
  if/else-if que fija el tipo de `bounda_time` (líneas 88–124). Variable de estado
  `bounda_on_off` (líneas 45 y 84). Aplicación periódica en el bucle de interpolación
  (líneas 214–217).
- `tochnog.h` — enum `BOUNDA_TIME_ON_OFF` (línea 141).
- `database.cc` — registro del keyword (líneas 183–187).

## Detalles de implementación

- Lectura: `db( BOUNDA_TIME_ON_OFF, iboun, ... )` devuelve `on_off_tmp[2]`
  (`on_time`, `off_time`). Valida `on_time > 0` y `off_time >= 0`, si no
  `db_error( BOUNDA_TIME_ON_OFF, iboun )`. Marca `bounda_on_off = 1`.
- Aplicación: dentro de la interpolación de `bounda_time` (rama `time`), cuando
  el segmento ya se encontró (`found`), calcula
  `phase = fmod( time_total, on_time + off_time )`; si `phase >= on_time` pone
  `load = 0.` y `found = 0` (carga cero en la ventana "off").
- IMPORTANTE: el bloque es INDEPENDIENTE de la cadena if/else-if de tipos de
  `bounda_time`. Insertarlo DENTRO de esa cadena sobreescribiría
  `ninc`/`time`/`length_bounda_time` y rompería la interpolación — bug conocido,
  no reintroducir. El bloque actual solo lee los valores y activa una bandera.
- `database.cc`: `type = DOUBLE_PRECISION`, `data_length = 2`,
  `data_class = BOUNDA`, `data_required = BOUNDA_TIME` (exige el keyword base).

## Dependencias externas

Ninguna. Solo la API interna del database y `fmod` de la librería estándar de C.

## Parámetros hardcodeados / refactorizaciones pendientes

- Los valores de `on_off_tmp[2]` y `until_tmp[2]` se leen con `ldum` (longitud no
  verificada explícitamente); confiar en `data_length = 2` de `database.cc`.
- Al marcar `found = 0` la carga "off" sigue buscando segmentos posteriores en el
  mismo `time`; el comportamiento en ciclos con varios segmentos debe revisarse.
- Los tres bloques nuevos de `bounda()` (on_off, until_force) duplican el patrón
  "leer bandera + validar + guardar en variable local"; se podría unificar la lectura.
