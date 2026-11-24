# control_print_frequency_timestep

## Implementación

- **Keyword**: `control_print_frequency_timestep` (INTEGER,
  data_length 1, data_class CONTROL, data_required CONTROL_TIMESTEP)
  registrado en `database.cc`. `ival[0]` es el número de pasos N. Solo
  el valor 0 o negativo da `db_error`.
- **Estado**: record dedicado `control_print_frequency_timestep_count`
  (INTEGER, data_length 1, external 0, data_class CONTROL) en
  `database.cc`: pasos desde el último print, por icontrol. Vive en
  VERSION_NORMAL (sobrevive a restart). Igual que el record de tiempo
  del timeinterval, NO lleva `no_index=1` (db() rechaza PUT con
  index>0 en records no-indexables; el estado es por icontrol).
- **Gate**: MISMO helper `control_print_frequency_allowed( icontrol )`
  de top.cc que el timeinterval (ver manual-developer de
  control_print_frequency_timeinterval). Rama timestep:
  - Se lee el contador con GET_IF_EXISTS (ausente → 0), se incrementa
    UNA vez por step_close (step_close corre una vez por paso: el bucle
    `while ( time_current<time_new )` de top() llama a step_close al
    final de cada paso, línea 472).
  - Permitido si `count >= N` o fin del incremento (misma detección que
    el timeinterval: `time_current >= time_new` con CONTROL_TIMESTEP
    activo). Al imprimir se resetea el contador a 0.
  - El fin de incremento SIEMPRE imprime y resetea el contador, así que
    cada incremento arranca su cadencia desde 0 (en el ejemplo del
    manual con 2 incrementos, la cadencia 5-5-1 se reproduce exacta).
- **Dispatch**: en step_close (top.cc) cada item `control_print_*`
  (EXCEPTO control_print, control_print_history + _history_smooth y
  control_print_data_versus_data) se envuelve con `frequency_allowed &&`.

## Enums nuevos

- `CONTROL_PRINT_FREQUENCY_TIMESTEP` y
  `CONTROL_PRINT_FREQUENCY_TIMESTEP_COUNT` (bloque CONTROL_PRINT, tras
  CONTROL_PRINT_DOF_ID) en `tochnog.h` / `tochnog-mod.h` (sync).

## GOTCHA del test (secuencia de pasos del GNU)

El bucle de timestep del GNU CLAMPEA el último paso parcial de un
incremento: con `control_timestep 10 0.04 0.41` la secuencia real es
0.04..0.36 y 0.41 (0.36+0.04 = 0.40, tmp = 0.01 < 0.5*dtime → el paso
se funde en 0.41). El ejemplo del manual escribe 0.40, que exige un
paso exacto en 0.40: el test `freq_timestep` usa DOS incrementos
(`0.04 0.40 0.01 0.01`) para reproducir la secuencia 0.04..0.40, 0.41
del manual. El timeinterval NO tiene este problema (cadencia por
tiempo, no por paso).

## Pendiente

- El contador no se reinicia explícitamente en restart (sigue desde el
  valor guardado); con un control_repeat que rebobina a un control
  anterior el conteo puede desfasarse. Caso no cubierto por el manual.
