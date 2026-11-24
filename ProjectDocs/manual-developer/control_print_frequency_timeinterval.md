# control_print_frequency_timeinterval

## Implementación

- **Keyword**: `control_print_frequency_timeinterval` (DOUBLE_PRECISION,
  data_length 1, data_class CONTROL, data_required CONTROL_TIMESTEP)
  registrado en `database.cc`. `dval[0]` es el intervalo de tiempo.
  Solo el valor 0 o negativo da `db_error` (intervalo inválido).
- **Estado**: record dedicado `control_print_frequency_timeinterval_time`
  (DOUBLE_PRECISION, data_length 1, external 0, data_class CONTROL) en
  `database.cc`: el último tiempo en que se imprimió, por icontrol.
  VIVE en VERSION_NORMAL (patrón de `control_print_gid_time`), así que
  sobrevive a restart; si `last_print_time > time_current` (restart que
  rebobina el tiempo) se re-ancla al inicio del incremento actual
  (TIME_OLD). **OJO con no_index**: el record de estado NO puede llevar
  `no_index=1` porque db() rechaza `PUT` con `index>0` para records
  no-indexables (control_print_gid_time usa index 0 por eso); el estado
  es por icontrol, así que es indexable (no_index 0) y external 0.
- **Gate**: helper `control_print_frequency_allowed( icontrol )` en
  `top.cc` (antes de step_close). Se llama UNA vez por step_close al
  inicio del bloque de dispatch de prints y devuelve 1/0:
  - Sin record de frecuencia para el icontrol → 1 (comportamiento
    previo, sin coste adicional).
  - FIN del incremento: se reutiliza la MISMA condición del bucle de
    timestep en top(): TIME_NEW se escribe al inicio de cada incremento
    de control_timestep y el último paso se clamp EAXACTAMENTE a él
    (`time_current = time_new`, top.cc:340-347), así que
    `time_current >= time_new - 1e-9*|time_new|` (con CONTROL_TIMESTEP
    activo para el icontrol y TIME_NEW existente) ⇒ fin del incremento
    ⇒ SIEMPRE se imprime.
  - timeinterval: permitido si `time_current >= last_print_time +
    interval`. El ANCLAJE inicial es TIME_OLD (inicio del incremento),
    NO el primer paso: el ejemplo del manual (intervalo 0.15, dt 0.04)
    imprime en 0.16, 0.32, 0.41 — anclar al primer paso (0.04) daría
    0.20. Al imprimir (por intervalo o por fin de incremento) se
    re-ancla `last_print_time = time_current`.
  - timestep: contador `control_print_frequency_timestep_count`
    (INTEGER, external 0) por icontrol, incrementado una vez por
    step_close (step_close corre una vez por paso); permitido si
    `count >= N` o fin de incremento; al imprimir se resetea a 0.
  - Si ambos records existen para el mismo icontrol, gana el
    timeinterval (decisión documentada; el manual no define el caso).
- **Dispatch**: en step_close (top.cc) cada item `control_print_*` del
  bloque (EXCEPTO `control_print`, `control_print_history` +
  `_history_smooth` y `control_print_data_versus_data`) se envuelve con
  `frequency_allowed &&`. NO se gatean cosas que no son control_print_*:
  `print_lastdatabase`, `exit_tn`, `post()`.

## Enums nuevos

- `CONTROL_PRINT_FREQUENCY_TIMEINTERVAL` y
  `CONTROL_PRINT_FREQUENCY_TIMEINTERVAL_TIME` (bloque CONTROL_PRINT,
  tras CONTROL_PRINT_DOF_ID) en `tochnog.h` / `tochnog-mod.h` (sync).

## Pendiente

- El gate evalúa por icontrol en step_close; los prints lanzados desde
  otras rutinas (p.ej. print_gid al final del cálculo, print_dx desde
  check_end) no pasan por el gate.
- Con restart, el contador de timestep NO se reinicia (el conteo sigue
  desde el valor guardado); el timeinterval sí se re-ancla (clamp). El
  comportamiento con restart se considera aceptable y queda documentado.
