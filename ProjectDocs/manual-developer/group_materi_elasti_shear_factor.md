# group_materi_elasti_shear_factor

## Implementación

- **Enum nuevo**: `GROUP_MATERI_ELASTI_SHEAR_FACTOR` en tochnog.h y
  tochnog-mod.h (en sync, entre `ELASTI_POISSON_POWER` y
  `ELASTI_SMALLSTRAIN`).
- **Registro** (database.cc): DOUBLE_PRECISION, length 1, data_class
  MATERI, data_required GROUP_TYPE: `factor`.
- **Bloque en `set_stress()`** (stress.cc), después del bloque
  poisson_power (cubre los young / young_polynomial / young_power /
  poisson_power):
  ```c
  if ( get_group_data( GROUP_MATERI_ELASTI_SHEAR_FACTOR, gr, element,
      new_unknowns, &shear_factor, ldum, GET_IF_EXISTS ) ) {
    for ( idim=0; idim<MDIM; idim++ )
      for ( jdim=0; jdim<MDIM; jdim++ ) {
        if ( idim==jdim ) continue;
        for ( kdim=0; kdim<MDIM; kdim++ )
          for ( ldim=0; ldim<MDIM; ldim++ ) {
            if ( kdim==ldim ) continue;
            C[idim][jdim][kdim][ldim] *= shear_factor;
            Cmem[idim][jdim][kdim][ldim] *= shear_factor;
          }
      }
  }
  ```
  Escala SOLO las entradas de corte de C y Cmem (el tensor 4º orden de
  C_matrix: los bloques (0,1),(0,2),(1,2) — las entradas diagonales
  (3,3),(4,4),(5,5) de Voigt en 3D y la (2,2) en 2D plane strain/plane
  stress — y sus gemelos simétricos). Las entradas normales (i==j o
  k==l) NO se tocan. Multiplicar ceros (p.ej. C[0][1][1][0]) es inofensivo.
- **check.cc**: requiere `materi_stress` + `materi_velocity`.

## Física

Manual 6.654: la rigidez cortante de young+poisson se multiplica por
factor (G_eff = G*factor). La tensión cortante se actualiza con C
(matrix_a4b en el predictor elástico) y la tangente ddsdde se ensambla
desde Cmem (stress.cc:1046) — ambas escaladas, así que tensión Y matriz
de rigidez del elemento ven el factor.

## Validación

- `mshf.dat` (factor 2) / `mshf_nof.dat` (sin record): cizalla pura
  (bottom -velx, top +velx), 1 paso dt=0.05, `E=1000 nu=0`, iteraciones 8.
  - sin record: sigma_xy = 33.33, eptxy = 0.03216;
  - factor 2: sigma_xy = 50, eptxy = 0.025;
  - factor 0: sigma_xy = 0 EXACTO (toda la tensión cortante fluye por
    las entradas escaladas — prueba dura de que el bloque actúa).
- A/B verificado también sin record == factor 1.0 (sin efectos laterales).

## Gotchas

- **La respuesta del dof post_point NO es lineal en el factor**: el dof
  de tensión es un unknown del sistema acoplado (formulación de
  tensiones: momentum + constitutiva + cinemática ensambladas juntas),
  y su valor resuelto responde a un cambio de tangente con
  sigma_xy(f) ~ 100*f/(f+2): ratio 1.5 con factor 2 (no 2.0). El dof de
  deformación también se desplaza (eptxy 0.0322 -> 0.025). El rig de
  fuerza (bounda_force + options_inertia -no) NO sirve como sonda de
  rigidez: el desplazamiento no responde a E ni al factor (el dof veli
  queda desacoplado de la rigidez del material en esa formulación).
  Documentado: la validación es la direccionalidad + el cero exacto +
  el A/B, no la proporcionalidad exacta.
- El orden de los bloques importa: shear_factor debe ir DESPUÉS de
  young/young_power/poisson_power (que recomputan C) y ANTES de las
  leyes volumétricas/camclay/lade (modelos separados que sobrescriben C).
