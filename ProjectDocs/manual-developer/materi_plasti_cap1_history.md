# materi_plasti_cap1_history

## Implementación

- **Enum**: `MATERI_PLASTI_CAP1_HISTORY` en tochnog.h y tochnog-mod.h
  (en sync; antes de `MATERI_PLASTI_F`, orden alfabético del enum de
  initia).
- **Globals** (initia.cc): flag `materi_plasti_cap1_history` e índice
  `cap1_indx` (patrón kap_indx/f_indx, inicializados a 0 y -1).
- **Parser** (input.cc, tras `materi_plasti_kappa`):
  ```c
  materi_plasti_cap1_history = 1;
  cap1_indx = unknown_indx;
  n = 1;
  array_set( &dof_type[cap1_indx], -MATERI_PLASTI_CAP1_HISTORY, n*nder );
  array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
  ```
- **Registro** (database.cc): name[] de la initia + basename `pc` en la
  tabla DOF_LABEL (tras el basename de `MATERI_PLASTI_KAPPA`), para que
  el dof sea direccionable en `target_item -post_point_dof ... -pc`.
- **Ecuación** (general.cc): casos `MATERI_PLASTI_CAP1_HISTORY` en
  `unknown_belongs_to_type` (`inertia = 1.`) y en `conv_part = 1.`
  (patrón exacto de `MATERI_PLASTI_KAPPA`).
- **RHS** (materi.cc): `element_rhside += volume*h*(new_cap1pc -
  old_cap1pc)/dtime` (patrón kappa); `old_cap1pc`/`new_cap1pc` se leen
  de old/new_unknowns y se pasan a `set_stress()` por referencia.
- **Actualización** (stress.cc, `set_stress()`): `new_cap1pc` se calcula
  tras las iteraciones plásticas con el `inc_epp` convergido (ver
  group_materi_plasti_cap1.md para la ley y la decisión de diseño).
- **Clamp** (dof.cc, `parallel_new_dof_diagonal`): `pc >= 0` tras cada
  paso (patrón del clamp de kap).

## Física

`pc` es la variable de historia del cap1 (posición del cap a lo largo
del eje p: `p*c = pc + c*cot(phi)`). Endurece con la deformación
volumétrica plástica del cap (forma inversa de la ley del manual):

```
pc_dot = deps_p_cv_dot * K_ref/(lambda*/kappa* - 1) * ((pc + c*cot(phi))/p_ref)^m
```

En el paso discreto: `new_cap1pc = old_cap1pc + deps_p_cv*K_ref/
(lambda*/kappa* - 1)*((old_cap1pc + c*cot(phi))/p_ref)^m` con
`deps_p_cv = -trace(inc_epp)` clampado a >= 0 (compresión endurece,
descarga/dilatación no). El valor inicial lo da el usuario en `node_dof`
(manual 4.17: "You need to give an initial value for it in the node_dof
records").

## Validación

- `mcap1.dat`: pc = 100 inicial -> 109.9996 tras 20 pasos plásticos
  (punto fijo discreto analítico 110.0) y constante en la descarga
  (target `-pc 110.4 ± 2.5`). Ver
  group_materi_plasti_cap1.md para la validación completa.

## Gotchas

- `node_dof` espera **nuknwn valores** (con `derivatives`, nder =
  1+ndim+1 slots por dof físico): en 3D con 25 dofs físicos son 125
  valores; el valor primario de `pc` va en el slot 90 (dof físico 18 =
  cap1_history, tras vel(3) veli(3) ept(6) epp(6)). Sin `derivatives`
  (nder=1) son 25 valores y `pc` va en la posición 19.
- La initia es obligatoria para `group_materi_plasti_cap1`: el bloque
  de plasti.cc sale con error si falta (patrón CAMCLAY) y check.cc
  falla el chequeo de combinación.
- El basename es `pc` (NO `materi_plasti_cap1_history`, que es el nombre
  de la initia y `array_member()` sobre dof_label no lo encuentra —
  mismo GOTCHA que kap vs materi_plasti_kappa).
