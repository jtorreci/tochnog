# group_materi_failure_crunching (+ alias del typo GNU)

## Implementación

- **Renombre del name[]**: `database.cc` —
  `name[GROUP_MATERI_FAILURE_CRUCHING]` cambia de
  `"group_materi_failure_cruching"` (typo GNU, sin la "n") a la
  ortografía Professional `"group_materi_failure_crunching"`.
- **Alias legacy en `db_number()`** (database.cc):

  ```c
  else if ( !strcmp( str, "group_materi_failure_cruching" ) )
    return GROUP_MATERI_FAILURE_CRUCHING;
  ```

  Ningún test (test-2014 ni sfnet) usa el typo, así que el renombre es
  seguro y el alias es solo por compatibilidad con inputs externos.
- **Física** (`failure()` en failure.cc, llamada desde top.cc por
  timestep): el criterio se activa con
  `db_partialname_any("group_materi_failure")`. Por elemento:
  `tmp = workval[2]` tras `matrix_jacobi` + `sort()` — **el autovalor
  principal MÁS COMPRESIVO** (sort() ordena descendente: val[0]=máx,
  val[2]=mín); si `tmp > 0` -> `tmp = 0` (tracción anula el criterio).
  Falla cuando `tmp > threshold`; el elemento se registra en
  `ELEMENT_DELETE_TIMES` [time_current, time_current+delete_time] y se
  borra gradualmente (`ELEMENT_DELETE_FACTOR` rampa 1 -> EPS_DELETE_FACTOR).
- **check.cc**: `GROUP_MATERI_FAILURE_CRUCHING` requiere
  `materi_strain_total`.

## Validación

- `mcrunch.dat` (threshold 0.5): oedometro `E=1000 nu=0.3 eps_zz=-0.001`
  -> tmp = -0.001 > 0.5 falso -> intacto, `sigma_xx = -0.5769` EXACTO.
- `mcrunch_low.dat` (threshold -0.5): el elemento se marca para borrado
  en el PRIMER timestep — `element_delete_times[0] = 0.05` (target sobre
  el record genérico `-element_delete_times 1 0`).

## Gotchas

- **Semántica de signo del threshold**: tmp <= 0 (compresión), luego
  `tmp > threshold` con threshold NEGATIVO dispara **incluso a
  deformación cero** (tmp = 0 > -0.5 en el primer failure()). Es el
  comportamiento fiel del GNU (idéntico al jan-2014: failure.cc:102-103,
  120) — documentado, NO corregido: el usuario debe dar el "crushing
  strain" como número negativo.
- `failure()` usa `db_dbl( NODE_DOF, inod, VERSION_NORMAL )` y hace la
  media por nodos (`array_add` + `1/nnol`); con `node_dof -ra N M -ra`
  (rango plano) los valores se asignan solo a algunos nodos — usar la
  sintaxis `-ra -from N -to M -ra` (lección aprendida en mvoid).
