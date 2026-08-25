# group_materi_elasti_young_power

## Implementación

- **Registro** (database.cc): `GROUP_MATERI_ELASTI_YOUNG_POWER`,
  DOUBLE_PRECISION, `data_length` **3 -> 6** (Sprint 10 lote 6),
  data_class MATERI, data_required GROUP_TYPE: `E0 E1 E2 E3 p1 alpha`.
- **Bloque en `set_stress()`** (stress.cc), reescrito sobre el patrón
  GNU:
  ```c
  if ( get_group_data( GROUP_MATERI_ELASTI_YOUNG_POWER, gr, element,
      new_unknowns, young_power, ldum, GET_IF_EXISTS ) ) {
    array_set( &C[0][0][0][0], 0., MDIM*MDIM*MDIM*MDIM );
    array_set( &Cmem[0][0][0][0], 0., MDIM*MDIM*MDIM*MDIM );  // <- vease Gotchas
    ...k0 hook...
    young0=young_power[0]; young1=young_power[1]; young2=young_power[2];
    young3=young_power[3]; p1=young_power[4]; alpha=young_power[5];
    if ( p1<=0. ) db_error( GROUP_MATERI_ELASTI_YOUNG_POWER, gr );
    p = - ( new_sig[0] + new_sig[4] + new_sig[8] ) / 3.;
    if ( control_materi_gate_off( CONTROL_MATERI_ELASTI_YOUNG_POWER_APPLY ) )
      young = young0;                       // 6.801: ley apagada
    else {
      young = young0 + young1 * scalar_power( p/p1, alpha );
      if ( young<young2 ) young = young2;   // E >= E2
      if ( young>young3 ) young = young3;   // E <= E3
    }
    C_matrix( young, poisson, ..., C, task );
    C_matrix( young, poisson, ..., Cmem, task );
  }
  ```
  **Decisión de signo de p**: `p = -sig_mean`, positiva en compresión
  (misma convención que el bloque GNU original y que el poisson_power
  del lote 5). La compresión (sig_mean < 0 -> p > 0) SUBE E. NO se usa
  `scalar_dabs` en `(p/p1)` (a diferencia de poisson_power): el manual
  Professional da `E = E0 + E1*(p/p1)^alpha` sin valor absoluto; la
  rama `E >= E2` protege el lado bajo para alpha=1. Para alpha NO
  entero con p < 0 (tracción), `scalar_power(p/p1, alpha)` devuelve NaN
  (dominio de pow) — comportamiento del modelo crudo; documentado.
- **check.cc**: requiere `materi_stress` + `materi_velocity` +
  `materi_strain_total` (heredado del registro GNU).

## Física

Manual 6.662 / teoría 2.2.2: `E = E0 + E1*(p/p1)^alpha`, caps
`E >= E2` y `E <= E3`; `p = -sig_mean`. La ley ES el módulo de Young:
sustituye a `group_materi_elasti_young` (no se suma).

## Validación

- `myoung6.dat`: oedometro E0=1000 E1=500 E2=800 E3=3000 p1=1 alpha=1,
  nu=0.3, eps_zz=-0.001, 100 pasos dt=0.005 + iteraciones 8.
  Punto fijo autoconsistente (`E = 1000+500p`, `p = E*eps/(3(1-2nu))`):
  E=1714.29, p=1.4286, sigma_zz=-2.3077, sigma_xx=-0.9890.
  Medido: sigma_xx=-0.713886, sigma_yy=-1.66573 -> E_eff ~ 1237 =
  **72% del fijo secante** (ver Gotchas).
- `myoung6_e2.dat`: E2=3000 -> E forzado a 3000: sigma_xx -1.7308 /
  sigma_yy -4.0385 EXACTO (cap E>=E2 activo, E constante -> elástico
  lineal exacto).
- `myoung6_e3.dat`: E3=800 -> E forzado a 800: sigma_xx -0.4615 /
  sigma_yy -1.0769 EXACTO (cap E<=E3 activo).
- `myoung6_apply.dat`: gate -no -> E=E0=1000: sigma_xx -0.5769 /
  sigma_yy -1.3462 EXACTO (la base lineal; A/B contra myoung6).

## Gotchas

- **CAMBIO SEMÁNTICO 3 -> 6 params**: el registro ya NO acepta la forma
  GNU `young0 * |p/p0|^alpha` (3 valores). Ningún test legacy ni de
  sfnet usa la forma vieja (verificado con grep en validation-suite y
  external-downloads/sfnet/extracted/test) — el upgrade no rompe nada.
- **GOTCHA MAYOR — C_matrix ACUMULA**: `C_matrix` hace
  `array_add(C_local, Ctot, Ctot)` (elasti.cc:166): NO sobreescribe.
  El GNU 2014 con `group_materi_elasti_young` + `young_power` presentes
  sumaba C(young) + C(E_power) -> rigidez DOBLADA (medido: 2x exacto en
  el oedometro). La semántica Professional es que young_power ES la ley
  (y el apply 6.801 usa E0 DEL PROPIO RECORD), así que el bloque
  young_power LIMPIA C y Cmem antes de construir. Verificado con
  probes: solo young_power -> -0.6980; young+young_power con la
  acumulación -> -1.3961 (= 2x); con la limpieza -> igual a solo
  young_power.
- **El estado final NO es el punto fijo secante** (72% en vez de 100%):
  la C se evalúa con la tensión del PASO PREVIO (acoplamiento
  explícito; new_sig se reinicializa desde old_unknowns en CADA
  iteración, así que las iteraciones del paso no actualizan la ley), y
  la tensión acumulada usa el PROMEDIO de los E pasados — la trayectoria
  queda por debajo del fijo (modelo E_{n+1} = 1.00417·E_n verificado a
  1.4%). Mismo patrón que mpower (84%, lote 5): GOTCHA documentado, no
  corregido.
- `control_timestep_iterations` > 1 NO cambia el resultado (la ley se
  evalúa con el stress del paso previo en todas las iteraciones).
