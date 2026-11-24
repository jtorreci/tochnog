# bounda_time_until_force

## Archivos y funciones

- `bounda.cc` — bloque dentro de `bounda()` (líneas 136–145) que lee la
  configuración; variable de estado `bounda_until_force` (línea 45).
  Aplicación (líneas 417–428) justo después de fijar
  `new_node_dof[iuknwn] = factor * load` (línea 413).
- `top.cc` — copia de las reacciones a `NODE_RHSIDE_PREVIOUS` al final de
  `step_close` (líneas 790–803).
- `tochnog.h` — enum `BOUNDA_TIME_UNTIL_FORCE` (línea 142) y
  `NODE_RHSIDE_PREVIOUS` (línea 751).
- `tochnog-mod.h` — enum espejo `NODE_RHSIDE_PREVIOUS` (línea 744), debe mantenerse
  en sincronía con `tochnog.h`.
- `database.cc` — registros `BOUNDA_TIME_UNTIL_FORCE` (líneas 189–193) y
  `NODE_RHSIDE_PREVIOUS` (líneas 3363–3366).

## Detalles de implementación

- Lectura: `db( BOUNDA_TIME_UNTIL_FORCE, iboun, ... )` devuelve `until_tmp[2]`
  (`until_force`, `until_factor`). Valida `0 < until_factor <= 1`, si no
  `db_error( BOUNDA_TIME_UNTIL_FORCE, iboun )`. Marca `bounda_until_force = 1`.
- Aplicación: solo para velocidades/desplazamientos prescritos
  (`bounda_until_force && force==0`, es decir rama `unknown`, no `bounda_force`).
  Lee la reacción del paso anterior con `db_dbl( NODE_RHSIDE_PREVIOUS, inod, ... )`
  y toma su valor absoluto en el índice `ireac = ( iuknwn - vel_indx ) / nder`.
  Si `reaction > until_force`:
  `new_node_dof[iuknwn] *= until_factor * ( until_force / reaction )`.
- `NODE_RHSIDE_PREVIOUS` es un item NUEVO de tipo `DOUBLE_PRECISION` con
  `data_length = npuknwn`. Se copia al final de `step_close` porque
  `iteration_start` (top.cc, línea 632) hace `db_set_dbl( NODE_RHSIDE, ... )`
  (resetea a 0) en cada iteración, y el bloque de `bounda()` corre DENTRO de la
  iteración. Solo se copia si `NODE_RHSIDE` está activo (`db_active_index` /
  `db_max_index`), por lo que el item solo existe cuando hay reacciones.
- `database.cc`: `BOUNDA_TIME_UNTIL_FORCE` es `DOUBLE_PRECISION`, `data_length = 2`,
  `data_class = BOUNDA`, `data_required = BOUNDA_TIME`.

## Dependencias externas

Ninguna. Solo la API interna del database (`db`, `db_dbl`, `db_active_index`) y
`scalar_dabs`/`fabs` estándar.

## Parámetros hardcodeados / refactorizaciones pendientes

- La reacción se lee con `scalar_dabs` (valor absoluto); un cambio de signo de la
  fuerza no se distingue (solo se limita la magnitud).
- `NODE_RHSIDE_PREVIOUS` se copia en `step_close` para TODOS los nodos con
  `NODE_RHSIDE` activo, incluso si no hay `bounda_time_until_force`; coste de
  memoria/copia en modelos grandes — podría condicionarse a que exista el keyword.
- El índice de reacción `ireac = ( iuknwn - vel_indx ) / nder` asume que la
  reacción está alineada con `vel_indx`; verificar para `-ROTATION_*_AXIS` y otros
  dof acoplados (los `rotate` del bloque de `BOUNDA_UNKNOWN`).
