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

## K06 fix (2026-10-06): the block replaces the C, it does not add to it

- **Where**: `set_stress()` in `stress.cc`, the
  `GROUP_MATERI_ELASTI_POISSON_POWER` block.
- **Defect**: `C_matrix()` ACCUMULATES into its target (`array_add`). The
  block called it on the `C`/`Cmem` already built by the
  `GROUP_MATERI_ELASTI_YOUNG` block, so the stiffness was
  `C(young, poisson) + C(young, nu(p))`. The "84 %" and "nu ~ 0.43" of the
  validation notes above were artefacts of that double stiffness; a
  forward-Euler replay of the oedometer with the summed matrix reproduces
  the old GNU numbers to 1e-9 (-1.43887 / -3.29236).
- **Fix**: `array_set()` clears `C` and `Cmem` before the rebuild, the
  same pattern as the `_young_power` and `_hardsoil` blocks. `young` keeps
  the modulus of the previous block (constant, polynomial or power law).
- **Reference**: Professional 02-08-2026 on `mpower` (converted) gives
  -0.44173 / -1.41230; the forward-Euler replay with `C(1000, nu(p))`
  gives -0.43887 / -1.41012, the GNU now -0.438872 / -1.410125. The
  residual 0.2-0.6 % is the Professional's pressure evaluation inside the
  step (not identified; below the target window).
- **Tests**: `validation-suite/test-2026/elasti_poisson_power_replace1.dat`
  (RED -1.4389/-3.2924, GREEN -0.4389/-1.4101), `mpower` retargeted to
  -0.4389 / -1.4101 +/- 0.02.
- **Pending**: the `_young_polynomial` block still accumulates on top of
  `group_materi_elasti_young` when both are given (not in any corpus); the
  `_young_power` + `_young` combination is an open decision (BACKLOG C36).
