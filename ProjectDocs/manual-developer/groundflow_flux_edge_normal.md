# groundflow_flux_edge_normal (familia)

## Files and functions

- `area.cc` — `area()` (integrales de área sobre bordes de elementos). La
  familia se añadió como `type[5] = GROUNDFLOW_FLUX_EDGE_NORMAL` con
  `type_area[5] = GROUNDFLOW_FLUX_EDGE_NORMAL_GEOMETRY` (`MTYPES` 5→6). Carga
  temporal por `_SINE`/`_TIME` (igual que `FORCE_ELEMENT_EDGE_NORMAL`); la
  aplicación al rhs del dof de presión:
  ```c
  force_factor( GROUNDFLOW_FLUX_EDGE_NORMAL_FACTOR, ind, &new_coord[inol*ndim], factor );
  db( GROUNDFLOW_FLUX_EDGE_NORMAL, ind, idum, values, ldum, VERSION_NORMAL, GET );
  ipuknwn = pres_indx/nder;
  tmp = factor * node_factor * load * weight[inol_side] * area_size * values[0];
  element_rhside[inol*npuknwn+ipuknwn] += tmp;
  ```
- Restricciones evaluadas ANTES del bucle de lados (`_element`,
  `_element_group`, `_element_side`) y por nodo dentro del bucle (`_node`,
  `_element_node`, `_element_node_factor`).
- `database.cc` — 11 keywords registradas (bloque GROUNDFLOW_FLUX_EDGE_NORMAL_*,
  entre `GROUNDFLOW_DENSITY` y `GROUNDFLOW_NONSATURATED_APPLY`).
- `check.cc` — requiere `groundflow_pressure` (y `ndim` 2-3 para la principal).
- Enums `GROUNDFLOW_FLUX_EDGE_NORMAL*` en `tochnog.h`.

## Implementation details

- El flux se añade al rhs del `pres_indx` (el dof de presión del groundflow);
  no a la velocidad. El signo: flux positivo inyecta agua.
- `_geometry` selecciona el área; la detección de aristas usa el mismo
  `border_nodes_*` y el cálculo de normal/área de las fuerzas de borde
  existentes.
- `_element_node_factor` es `DOUBLE_PRECISION`: `values_fac[0]=element`,
  `values_fac[1..]=factores` por nodo local (leído con `db(..., idum,
  values_fac, ...)`).
- `_element_side` guarda pares `(elemento, lado)`; se comprueba por elemento.
- El `continue` de las restricciones salta el registro completo (no se evalúa
  ninguna arista).

## External dependencies

- `force_factor()`, `force_time()` (force.cc), `array_member()`, db accessors.
- Globals `pres_indx`, `ndim`, `npuknwn`, `time_total`.

## Hardcoded parameters / pending refactorings

- `_element_side` comprueba solo por elemento (no por número de lado); el lado
  concreto se deja al `_geometry`. Refinamiento posible.
- El signo del flux (entrada = positivo) se tomó de la convención de fuerzas de
  borde; verificar contra un ejemplo de Professional si aparece.
