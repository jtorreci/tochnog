# group_materi_failure_void_fraction (+ alias del nombre GNU)

## Implementación

- **Renombre del name[]**: `database.cc` —
  `name[GROUP_MATERI_FAILURE_VOIDFRACTION]` cambia de
  `"group_materi_failure_voidfraction"` (sin underscore) a la ortografía
  Professional `"group_materi_failure_void_fraction"`.
- **Alias legacy en `db_number()`** (database.cc):

  ```c
  else if ( !strcmp( str, "group_materi_failure_voidfraction" ) )
    return GROUP_MATERI_FAILURE_VOIDFRACTION;
  ```

- **Física** (`failure()` en failure.cc): `tmp = |void_fraction|`
  (dof `materi_void_fraction`, `void_indx`); falla cuando
  `tmp > threshold` -> `ELEMENT_DELETE_TIMES` + borrado gradual con
  `delete_time`. La fracción de vacíos evoluciona con la deformación
  plástica volumétrica en materi.cc
  (`tmp = volume*h*(1-void)*void*(inc_epp volumétrico)`); sin plasticidad
  se queda en su valor inicial (node_dof).
- **check.cc**: `GROUP_MATERI_FAILURE_VOIDFRACTION` requiere
  `materi_void_fraction`.

## Validación

- `mvoid.dat` (threshold 0.9): void inicial 0.3 via
  `node_dof -ra -from 1 -to 4 -ra` (11 valores por nodo: velx vely velix
  veliy + 6 de stress + void) -> `0.3 > 0.9` falso -> intacto,
  `sigma_xx = -0.5769` EXACTO.
- `mvoid_low.dat` (threshold 0.2): `0.3 > 0.2` en el primer timestep ->
  `element_delete_times[0] = 0.05` (target genérico
  `-element_delete_times 1 0`).

## Gotchas

- **node_dof con rango**: `node_dof -ra 1 4 -ra` asigna los valores solo
  a los nodos 1 y 4 (parser); la forma correcta es
  `node_dof -ra -from 1 -to 4 -ra`. Diagnóstico con print del
  average_node_dof dentro de failure() (void 0.15 = media de nodos 0.3/0).
- El post_point interpola el dof escalar void dando 0.15 (mitad del valor
  nodal 0.3) — el criterio de fallo usa los valores NODALES
  (average_node_dof), que son los correctos.
- El void dof participa en el solve (unknown): su residual es 0 sin
  plasticidad, por lo que se conserva el valor inicial.
