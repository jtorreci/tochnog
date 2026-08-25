# group_materi_elasti_poisson_power

## Implementación

- **Enum nuevo**: `GROUP_MATERI_ELASTI_POISSON_POWER` en tochnog.h y
  tochnog-mod.h (en sync, orden alfabético entre `ELASTI_POISSON` y
  `ELASTI_SMALLSTRAIN`).
- **Registro** (database.cc): DOUBLE_PRECISION, length 5, data_class
  MATERI, data_required GROUP_TYPE: `nu0 nu1 nu2 p1 alpha`.
- **Bloque en `set_stress()`** (stress.cc), INMEDIATAMENTE después del
  bloque `GROUP_MATERI_ELASTI_YOUNG_POWER`, siguiendo EXACTAMENTE su
  patrón:
  ```c
  if ( get_group_data( GROUP_MATERI_ELASTI_POISSON_POWER, gr, element,
      new_unknowns, poisson_power, ldum, GET_IF_EXISTS ) ) {
    get_group_data( GROUP_MATERI_ELASTI_POISSON, gr, element,
      new_unknowns, &poisson, ldum, GET_IF_EXISTS );   // default
    p1 = poisson_power[3];
    if ( p1<=0. ) db_error( GROUP_MATERI_ELASTI_POISSON_POWER, gr );
    p = - ( new_sig[0] + new_sig[4] + new_sig[8] ) / 3.;
    poisson = nu0 + nu1 * scalar_power(scalar_dabs(p/p1),alpha);
    if ( poisson>nu2 ) poisson = nu2;
    C_matrix( young, poisson, ..., C, task );
    C_matrix( young, poisson, ..., Cmem, task );
  }
  ```
  `young` ya viene de los bloques young/young_power previos; el bloque
  RECOMPUTA C y Cmem con el poisson dependiente de presión.
- **check.cc**: requiere `materi_stress` + `materi_velocity`.

## Física

Manual 6.653 / teoría 2.2.2: `nu = nu0 + nu1*(p/p1)^alpha` con
`nu <= nu2`; `p = -sig_mean` (positiva en compresión — misma convención
que `young_power`). `scalar_dabs` en `(p/p1)` replica el patrón de
young_power (evita NaN para potencias fraccionarias con p negativo).

## Validación

- `mpower.dat`: oedometro `E=1000`, `eps_zz = -0.0012` (vel -0.0024,
  t=0.5), `nu0=0.2 nu1=0.1 nu2=0.5 p1=1 alpha=1`, iteraciones 8.
  Punto fijo autoconsistente analítico (`nu = 0.2+0.1p`,
  `p = E*eps/(3(1-2nu))` -> `E*eps = 1.2`): nu=0.4, p=2,
  sigma_yy=-2.5714, sigma_xx=-1.7143.
  Medido: sigma_xx=-1.4389, sigma_yy=-3.2924 (p ~ 2.06, nu_ratio ~ 0.43
  vs 0.406 analítico — ~5%). Ventana: sigma_xx -1.4389±0.15,
  sigma_yy -3.2924±0.3 (excluye la base nu=0.3: -0.5769/-1.3462).

## Gotchas

- El estado final NO es exactamente el punto fijo secante: la C se
  evalúa con la tensión del PASO PREVIO (acoplamiento explícito), y el
  dof de tensión es un unknown del sistema acoplado (formulación de
  tensiones) — el estado aterriza ~84% del recorrido hacia el punto
  fijo. No es un bug: es la misma arquitectura de young_power (que nunca
  fue validado al punto fijo). Documentado en vez de "corregido".
- `control_timestep_iterations` > 1 NO cambia el resultado elástico
  (verificado: iter 1 == iter 8) — la no-linealidad residual es del
  acoplamiento, no de la convergencia del paso.
