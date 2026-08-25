# group_materi_elasti_stress_pressure_history_factor

## Implementación

- **Enums nuevos** (tochnog.h y tochnog-mod.h, en sync): el dof
  `MATERI_STRESS_PRESSURE_HISTORY` (enum de initia, entre
  `MATERI_STRESS` y `MATERI_VELOCITY`) y el group
  `GROUP_MATERI_ELASTI_STRESS_PRESSURE_HISTORY_FACTOR` (entre
  `ELASTI_SMALLSTRAIN` y `ELASTI_TRANSVERSE_ISOTROPY`).
- **Globals** (initia.cc): flag `materi_stress_pressure_history` e
  índice `sph_indx` (patrón kap_indx).
- **Parser** (input.cc): initia `materi_stress_pressure_history`
  (manual 4.50) — n=1, `dof_type -MATERI_STRESS_PRESSURE_HISTORY`,
  `dof_scal_vec_mat -SCALAR`.
- **Registros** (database.cc): initia (name) + basename `sph` en la
  tabla DOF_LABEL (para `target_item -post_point_dof ... -sph`) +
  group DOUBLE_PRECISION length 1 MATERI GROUP_TYPE.
- **Historia del dof** (dof.cc, `parallel_new_dof_diagonal`, patrón del
  clamp de kap):
  ```c
  if ( materi_stress_pressure_history && materi_stress ) {
    iuknwn = sph_indx;
    tmp = scalar_dabs( -( node_dof_new[stres_indx] +
      node_dof_new[stres_indx+3*nder] + node_dof_new[stres_indx+5*nder] ) / 3. );
    if ( node_dof_new[iuknwn] < tmp ) node_dof_new[iuknwn] = tmp;
  }
  ```
  node_dof_new fue copiado de VERSION_NORMAL al inicio del paso, así que
  es `max(viejo, |p_nuevo|)` — el máximo histórico running. El dof NO
  tiene ecuación (no está en general.cc: mismo patrón que
  materi_velocity_integrated) — el solver no lo toca; dof.cc es su único
  escritor. `stres_indx+{0,3,5}*nder` = sigxx/sigyy/sigzz (con nder
  correcto si hay derivatives).
- **Bloque en `set_stress()`** (stress.cc), DESPUÉS del bloque
  shear_factor (punto de aplicación: la C/Cmem elástica FINAL, así
  escala la rigidez completa de young/young_power/poisson_power y
  combina con shear_factor):
  ```c
  if ( materi_stress_pressure_history &&
       get_group_data( GROUP_MATERI_ELASTI_STRESS_PRESSURE_HISTORY_FACTOR,
         gr, element, new_unknowns, &sph_factor, ldum, GET_IF_EXISTS ) ) {
    p = - ( new_sig[0] + new_sig[4] + new_sig[8] ) / 3.;
    matrix_a4b( C, inc_ept, work );          // work = C:inc_ept
    p += -( work[0] + work[4] + work[8] ) / 3.;   // p del paso ACTUAL
    if ( p==0. ) p = 0.;                     // normaliza IEEE -0.0
    if ( scalar_dabs(p) < old_unknowns[sph_indx] ) {
      C *= sph_factor; Cmem *= sph_factor;   // descarga/recarga
    }
  }
  ```
- **check.cc**: el group requiere `materi_stress` + `materi_velocity` +
  `materi_stress_pressure_history`; la initia requiere `materi_stress`.

## Física

Manual 6.655: mientras la presión actual es MENOR que el máximo
histórico el material está descargando/recargando -> rigidez × factor;
si es el nuevo máximo, se actualiza la historia y NO se multiplica.
`p = -sig_mean` (compresión positiva).

**Decisión de diseño — presión del paso actual y máximo a inicio de
paso**: los bloques elásticos evalúan la rigidez con la tensión del
PASO PREVIO (new_sig parte de old_unknowns), así que la decisión usa la
presión que el paso actual VA a alcanzar: `p_est = p_old + dp` con
`dp = -mean(C:inc_ept)` (el incremento elástico de presión del paso).
Esto captura el PRIMER paso de descarga (con p_old solo, en el pico
`p_old == sph` y el factor arrancaría un paso tarde). El máximo
histórico se lee de **old_unknowns[sph_indx]** (VERSION_NORMAL, valor a
inicio de paso): durante el paso dof.cc sube sph con el |p| corriente, y
la decisión debe comparar contra la historia EXCLUYENDO el paso actual —
si se compara contra new_unknowns[sph] (subido medio paso), en el pico
`p_est == sph` EXACTAMENTE y cualquier redondeo decide carga vs
descarga: factor espurio -> la tensión se dispara -> sph se infla ->
runaway (observado: sph 0.3333 -> 0.5 y stress 6x en el paso del pico).

## Validación

- `msph.dat`: oedometro E=1000 nu=0.3, factor 3. Fase carga: 4 pasos
  vely=-0.002 (pico sigyy -0.5385, |p|=0.3333). Fase descarga: 2 pasos
  vely=+0.002. Con factor 3 la tangente de descarga es 3E: la tensión
  SOBREPASA a +0.2692 (recupera 3x0.2692 desde el pico); sph se queda
  en 0.3333 (target -sph 0.3333±0.03). Paso a paso verificado:
  pico -0.5385 (1E por paso), un paso de descarga -> -0.1346 (3E),
  dos -> +0.2692.
- `msph_flat.dat`: mismo input con factor 1 -> -0.2692 (A/B: la
  descarga recupera elásticamente 1E por paso); sph 0.3333.
- El factor NO se aplica durante la carga (verificado paso a paso:
  1E, 2E, 3E, 4E hasta el pico — sin factor espurio).

## Gotchas

- **scalar_dabs(-0.0) devuelve -0.0**: `if (a<0.) result=-a; else
  result=a;` — para -0.0 la rama `a<0.` es falsa (IEEE) y devuelve -0.0;
  `-0.0 < sph` es TRUE para cualquier sph -> el factor se aplicaría
  desde el primer paso de carga (el stress viejo es 0 y p = -0.0).
  Normalizado con `if (p==0.) p = 0.;` ANTES del scalar_dabs.
- **matrix_a4b NO es in-place seguro**: `matrix_a4b(C, work, work)`
  corrompe el buffer de salida (escribe C[0][0] y luego relee work[0][0]
  ya sobrescrito). Usar `matrix_a4b(C, inc_ept, work)` (fuente intacta).
- **El dof sph no tiene ecuación**: no está en general.cc
  (unknown_belongs_to_type) — patrón de materi_velocity_integrated. El
  solver no lo toca; solo dof.cc lo escribe. Con `derivatives`, el valor
  está en el slot primario (sph_indx), la derivada en sph_indx+1.
- La decisión usa `inc_ept` (incremento TOTAL, incluye plástico si lo
  hay) como estimación del incremento elástico de presión; el dof sph
  se actualiza con el |p| REAL resuelto, así que cualquier desviación de
  la estimación se autocorrige en el paso siguiente.
- El target del dof usa el basename `-sph` (patrón del GOTCHA -kap del
  lote 4).
