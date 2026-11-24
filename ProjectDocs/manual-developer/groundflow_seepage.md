# groundflow_seepage_geometry / groundflow_seepage_node / groundflow_seepage_eps

## Files and functions

- `bounda.cc` — in `bounda()` (procesamiento de `bounda_dof`/`bounda_unknown`),
  just after `node_bounded[ipuknwn] = 1`. When the dof is `pres_indx` and
  `GROUNDFLOW_SEEPAGE_EPS` exists, the seepage check runs:
  ```c
  if ( iuknwn==pres_indx && groundflow_pressure &&
       db_active_index( GROUNDFLOW_SEEPAGE_EPS, 0, VERSION_NORMAL ) ) {
    db( GROUNDFLOW_SEEPAGE_EPS, 0, idum, &seep_eps, ldum, ... );
    // node list (GROUNDFLOW_SEEPAGE_NODE) or geometry (GROUNDFLOW_SEEPAGE_GEOMETRY)
    // determine if the node is on a seepage edge (geometry() gives the outward
    // normal of the seepage geometry).
    // flow_dot = gvel . seep_normal (Darcy velocity of the previous step)
    if ( in_seep && flow_dot < -seep_eps )
      node_bounded[ipuknwn] = 0;   // flow entering: close the edge
  }
  ```
- `database.cc` — 3 keywords registradas (bloque GROUNDFLOW_SEEPAGE_*, entre
  `GROUNDFLOW_SATURATION` y `GROUNDFLOW_VELOCITY`). `_eps` es `no_index`;
  `_geometry` guarda `[entity, index]`; `_node` es variable-length.
- `check.cc` — `_geometry`/`_node` requieren `groundflow_pressure`.
- Enums `GROUNDFLOW_SEEPAGE_*` en `tochnog.h`.

## Implementation details

- El seepage se aplica SOLO cuando existe `GROUNDFLOW_SEEPAGE_EPS` (es el
  gatillo de activación); el manual da default 0.1.
- La normal exterior sale de `geometry()` sobre el geometry del
  `_geometry`; para `_node` no hay normal propia (el flujo se comprueba con la
  normal del geometry del bounda asociado, vía `geometry()` con el mismo
  `geometry_entity`).
- `flow_dot` usa la velocidad de Darcy del paso anterior (`NODE_DOF` VERSION_NORMAL,
  componente `gvel_indx`), que apunta en la dirección del flujo. `flow_dot>0`
  → flujo saliente (normal apunta afuera); `flow_dot < -eps` → flujo entrante
  → el borde se cierra (no se impone la presión).
- El check se hace ANTES de imponer `new_node_dof[iuknwn]`, desmarcando
  `node_bounded` para que el dof quede libre (borde cerrado).

## External dependencies

- `geometry()`, `array_normalize()`, `array_member()`, db accessors.
- Globals `groundflow_pressure`, `groundflow_velocity`, `pres_indx`,
  `gvel_indx`, `nder`, `ndim`.

## Hardcoded parameters / pending refactorings

- El criterio usa la velocidad de Darcy del paso anterior; en el primer paso
  es 0, por lo que el borde no se cierra hasta que el flujo se establece.
- `_node` no lleva normal propia; depende del geometry del bounda asociado.
