# materi_stress_pressure_history

## Implementación

- **Enum**: `MATERI_STRESS_PRESSURE_HISTORY` en tochnog.h y
  tochnog-mod.h (en sync; entre `MATERI_STRESS` y `MATERI_VELOCITY`,
  orden alfabético del enum de initia).
- **Globals** (initia.cc): flag `materi_stress_pressure_history` e
  índice `sph_indx` (patrón kap_indx/f_indx).
- **Parser** (input.cc, entre `materi_stress` y `materi_velocity`):
  ```c
  materi_stress_pressure_history = 1;
  sph_indx = unknown_indx;
  n = 1;
  array_set( &dof_type[sph_indx], -MATERI_STRESS_PRESSURE_HISTORY, n*nder );
  array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
  ```
- **Registro** (database.cc): name[] de la initia + basename `sph` en la
  tabla DOF_LABEL (tras el basename del `MATERI_STRESS`), para que el
  dof sea direccionable en `target_item -post_point_dof ... -sph` y en
  post-proceso.
- **Actualización del dof** (dof.cc, `parallel_new_dof_diagonal`,
  patrón del clamp de kap, ~línea 185): cada paso (y cada iteración del
  paso) después del solve,
  ```c
  node_dof_new[sph_indx] = max(node_dof_new[sph_indx], |p_new|)
  ```
  con `p_new = -(sigxx+sigyy+sigzz)/3` leído de
  `node_dof_new[stres_indx + {0,3,5}*nder]` (componentes Voigt:
  sigxx/sigyy/sigzz — ver la tabla de basenames del MATERI_STRESS).
  node_dof_new se copió de VERSION_NORMAL al inicio del paso, así que el
  máximo es running (`max(viejo, |p|)`).
- **check.cc**: la initia requiere `materi_stress` (la presión se
  computa de los dofs de tensión; sin materi_stress, stres_indx = -1 ->
  crash).

## Física

Manual 4.50: "The maximum of the absolute value of the pressure which
occurs over time is added to the node_dof records." `p = -sig_mean`,
positiva en compresión (misma convención que
`group_materi_elasti_young_power`). El dof es un máximo running: nunca
decrece, así que marca la presión pico histórica en cada punto.

## Validación

- `msph.dat` / `msph_flat.dat` (ver
  group_materi_elasti_stress_pressure_history_factor.md): el dof llega a
  0.3333 (el pico |p| de la fase de carga) y se MANTIENE en 0.3333
  durante la descarga — target `-sph 0.3333 ± 0.03` en ambos.

## Gotchas

- El dof NO tiene ecuación en general.cc (no está en
  `unknown_belongs_to_type`): es el patrón de
  materi_velocity_integrated. Sin lhside, el solver no lo toca; el único
  escritor es dof.cc. Si algún día necesita interpolación/derivadas
  espaciales, habrá que añadirlo a la tabla de general.cc (conv_part=1).
- Con `derivatives` (nder=2), el VALOR está en el slot primario
  `sph_indx`; el slot `sph_indx+1` es la derivada temporal (computada
  por el bloque `derivatives` de dof.cc).
- El consumo (factor de descarga) vive en set_stress() y lee
  `old_unknowns[sph_indx]` (máximo a INICIO de paso) — ver
  group_materi_elasti_stress_pressure_history_factor.md para la
  justificación.
