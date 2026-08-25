# control_repeat_save_calculate

## Implementación

- **Record**: `control_repeat_save_calculate` (INTEGER, length 1,
  data_class CONTROL, `data_required = CONTROL_REPEAT`) registrado en
  `database.cc`.
- **Lógica**: en `repeat.cc`, en la rama `CONTROL_REPEAT` de
  `repeat()`: cuando el contador del repeat llega a 0 (la pasada final,
  `repeat_return == 0`, no hay salto), se llama a la función estática
  `control_repeat_save_calculate( icontrol )`:
  1. `ndata = db_len(CONTROL_REPEAT_SAVE)/3` (si `ndata < 1` no hace
     nada).
  2. Cuenta los índices ACTIVOS de `REPEAT_SAVE_RESULT` (mismo GOTCHA
     de `db_max_index` que en `control_repeat_save` — ver su manual) →
     `nsave`. Si `nsave < 1` no hace nada (un `control_repeat 0` no
     habría guardado nada).
  3. Lee todas las filas `[nsave × ndata]` y por cada data item
     calcula media y varianza en DOS pasadas (media primero, luego
     Σ(x−mean)² — estable numéricamente, sin cancelación catastrófica
     de `E[x²]−mean²`).
  4. Escribe `REPEAT_CALCULATE_RESULT[idata] = [mean, variance]`.

## Layout de REPEAT_CALCULATE_RESULT

- `repeat_calculate_result` (DOUBLE, length 2 fija — un índice por
  data item, cada uno con `[average, variance]`; data_class CALCUL).
- Varianza POBLACIONAL: Σ(xᵢ−mean)² / n (la varianza descriptiva del
  conjunto de datos; el manual Professional 6.349 solo dice "average
  value and variance", sin fórmula — decisión documentada).
- Se escribe SOLO cuando hay al menos un save y el repeat ha
  completado; si el record no existe, no se crea (A/B verificado).

## PENDIENTE

- No se soporta con `control_repeat_until_*`.
- Repeats anidados: se re-ejecuta al completar cada repeat interior
  (último write gana sobre `REPEAT_CALCULATE_RESULT`), no validado.
- El resultado se lee vía `target_item` / `.dbs`; no hay
  `control_print` específico (decisión: records de resultado, no de
  salida).
