# groundflow_phreatic_level_multiple (familia)

## Files and functions

- `groundfl.cc`:
  - `groundflow_phreatic_level_multiple_find( inod )` — devuelve el índice del
    registro multiple que es dueño del nodo (o -1). El dominio se selecciona
    por exactamente uno de `_element`, `_element_group`, `_element_geometry` o
    `_node` (no combinables, según el manual). El nodo pertenece al dominio si
    está en `_node`, o si uno de sus elementos (`NODE_ELEMENT`) está en
    `_element`, en `_element_group` (via `ELEMENT_GROUP`) o cumple la geometría
    `_element_geometry` (todos o cualquiera de los nodos del elemento en la
    geometría).
  - `groundflow_phreatic_coord()` — si existe `GROUNDFLOW_PHREATICLEVEL_MULTIPLE`,
    usa el nivel del nodo (vía `_find`) en lugar del único
    `GROUNDFLOW_PHREATICLEVEL`; el resto de la rutina (static/total pressure,
    clamps) es idéntico.
  - `groundflow_phreatic_apply()` — bloque final: para cada índice con
    `GROUNDFLOW_PHREATICLEVEL_MULTIPLE_STATIC -yes`, impone en los nodos del
    dominio `node_dof[pres_indx] = static_pressure` y marca el dof como bounded.
- `database.cc` — 7 keywords registradas como `groundflow_phreatic_level_multiple*`
  (con guion, coincidiendo con el manual de Professional 2024; la única sigue
  siendo `groundflow_phreaticlevel` sin guion, convención del GNU).
- `check.cc` — requiere `groundflow_pressure`; `_n` además `ndim==3`.
- Enums `GROUNDFLOW_PHREATICLEVEL_MULTIPLE*` en `tochnog.h`.

## Implementation details

- `_find` lee `NODE_ELEMENT` (lista de elementos por nodo) y `ELEMENT_GROUP`
  por elemento, y `geometry()` para `_element_geometry`.
- `_static` se aplica en `groundflow_phreatic_apply` (barrido de nodos), no
  dentro de `groundflow_phreatic_coord`, porque requiere poner el dof como
  bounded.
- El `_static` usa el `static_pressure` ya computado por `groundflow_phreatic_coord`
  (que resuelve el nivel del nodo), consistente con el único.

## External dependencies

- `geometry()`, `table_xy()`, `table_xyz()`, db accessors.
- Globals `pres_indx`, `ndim`, `nder`, `groundflow_pressure`.

## Hardcoded parameters / pending refactorings

- `_find` hace un barrido O(nodos × niveles) dentro de `groundflow_phreatic_apply`;
  para mallas muy grandes podría cachearse la pertenencia por elemento.
- El `_element_geometry` acepta "all OR any" nodo del elemento en la geometría
  (patrón `geometry_method`); el manual no especifica el método.
