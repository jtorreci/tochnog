# groundflow_phreatic_level_multiple (familia)

## Files and functions

- `groundfl.cc`:
  - `groundflow_phreatic_level_multiple_active()` — TRUE cuando existe AL MENOS
    un registro `GROUNDFLOW_PHREATICLEVEL_MULTIPLE` en cualquier índice
    (`db_max_index >= 0`). GOTCHA corregido 2026-09-07: los registros del
    corpus (ground8/ground19) viven en los índices 10/20/30, y la rama multiple
    de `groundflow_phreatic_coord` y `groundflow_phreatic_apply` se gateaba con
    `db_active_index(...,0,...)` (solo ve un registro en el índice 0, el layout
    del test interno `groundflow_phreatic_multiple`) → con niveles en índices
    != 0 la familia entera estaba MUERTA (ground8 evaluaba to_pres = dof − ρgz
    en vez de dof + static del nivel: −60 vs −10 del Pro).
  - `groundflow_phreatic_level_multiple_find_element( elnum )` — devuelve el
    índice del nivel dueño de UN elemento (o -1). Núcleo de pertenencia
    compartido: `_element` (miembro de la lista), `_element_group` (grupo del
    elemento), `_element_geometry` (todos o cualquiera de los nodos del
    elemento dentro de la geometría, patrón `geometry_method`); un dominio
    `_node` no puede ser dueño de un elemento. Barrido ascendente: gana el
    primer (menor) índice.
  - `groundflow_phreatic_level_multiple_find( inod )` — resuelto SOBRE
    `find_element`: primero los dominios `_node` (barrido ascendente), luego el
    mínimo índice entre los dueños de los elementos del nodo (equivalente al
    orden del barrido original por niveles). GOTCHA: los nodos sin registro
    `NODE_ELEMENT` (mallas control_mesh_macro donde el registro solo se guarda
    para nodos con elementos) NO pertenecen a ningún dominio → return -1 (sin
    el guard, `db_len` de un índice inactivo → "Error detected for data item:
    node_element").
  - `groundflow_phreatic_level_multiple_find_coord( coord[] )` — para las
    evaluaciones con inod<0 (post points): localiza el elemento que contiene la
    coordenada (point_el sobre todos los elementos activos, frame
    NODE_START_REFINED, igual que `parallel_post_point`) y devuelve
    `find_element` del elemento hallado. Solo se usa con inod<0 (coste O(malla)
    por llamada; los paths por nodo usan `_find`, sin localización).
  - `groundflow_phreatic_coord()` — la rama multiple se activa con
    `_multiple_active()` (cualquier índice) y resuelve el nivel del nodo
    (inod>=0, vía `_find`), del punto (inod<0, vía `_find_coord`); el resto
    (static/total pressure, clamps a la atmosférica) es idéntico.
  - `groundflow_phreatic_apply()`:
    - la condición de exclusión del bloque de nivel SIMPLE usa
      `!_multiple_active()` (antes `!db_active_index(...,0,...)`: con niveles
      multiple en índices != 0 el bound de superficie libre del nivel simple se
      habría aplicado encima de los dominios multiple);
    - BLOQUE NUEVO: superficie libre POR NIVEL multiple SIN `_static -yes`,
      análoga al nivel simple: para cada nodo del dominio del nivel, si
      `coord >= level − EPS` se boundea `pres_dof = 0` (zona seca; el static de
      `phreatic_coord` se clampea a 0 sobre el nivel) → la zona saturada queda
      confinada bajo el nivel con su Dirichlet de superficie libre;
    - bloque `_static -yes`: gateado con `_multiple_active()` (antes índice 0).
- `materi.cc` (~línea 560) — sustitución de la presión de poro en la mecánica:
  `new_pres = total_pressure` del nivel SOLO cuando NO hay registros multiple
  (`!groundflow_phreatic_level_multiple_active()`). La rama multiple queda
  PENDIENTE de calibración (ver "pending"): con dominios multiple + cargas de
  head `bounda_dof -pres` (estáticas de embalse de ground19), sumar el static
  del nivel al head resuelto duplica ρg·z_L y rompe el sistema acoplado
  (Bi-CG breakdown → retry LU que no termina). Medido: con la sustitución
  activa ground19 entra en breakdown; sin ella rc=0 con el MISMO flujo
  (−0.00557923). El nivel simple (ground15/16) mantiene la sustitución.
- `database.cc` — 7 keywords registradas como `groundflow_phreatic_level_multiple*`
  (con guion, coincidiendo con el manual de Professional 2024; la única sigue
  siendo `groundflow_phreaticlevel` sin guion, convención del GNU).
- `check.cc` — requiere `groundflow_pressure`; `_n` además `ndim==3`.
- Enums `GROUNDFLOW_PHREATICLEVEL_MULTIPLE*` en `tochnog.h`/`tochnog-mod.h`
  (sincronizados) + prototipo de `groundflow_phreatic_level_multiple_active`
  en ambos headers (lo consume materi.cc).

## Implementation details

- La pertenencia por elemento vive en `find_element` (reutilizada por `_find`,
  `_find_coord` y el futuro uso por elemento en materi()); `_find` ya no
  duplica la lógica de selectores.
- `_static` y la superficie libre por nivel se aplican en
  `groundflow_phreatic_apply` (barrido de nodos), no dentro de
  `groundflow_phreatic_coord`, porque requieren poner el dof como bounded.
- La superficie libre por nivel se salta los niveles con `_static -yes` (el
  modo estático prescribe la presión sobre TODO el dominio, no hay zona seca
  que boundear).

## External dependencies

- `geometry()`, `table_xy()`, `table_xyz()`, `point_el()`, db accessors.
- Globals `pres_indx`, `ndim`, `nder`, `groundflow_pressure`.

## Hardcoded parameters / pending refactorings

- `_find_coord` localiza por fuerza bruta (point_el sobre todos los elementos);
  solo corre para evaluaciones con inod<0 (post points, pocas por paso). Si un
  día se evalúan muchos puntos contra mallas grandes, cachear el último
  elemento acertado o indexar por cajas.
- El acoplamiento MECÁNICO de los niveles múltiples (sustitución
  total_pressure en materi.cc para dominios multiple) está PENDIENTE de
  calibración: requiere alinear primero la convención del dof pres entre los
  loads `-pres` (head completo) y el split dinámico del nivel (ver registro de
  ground8 en SEGUIMIENTO-CONVERGENCIA.md). ground8 (targets de to_pres, rc=0)
  NO necesita la sustitución: sus pres dofs se resuelven a 0 (agua en reposo
  bajo el nivel) y el perfil hidrostático sale del static del nivel en el
  post-cálculo.
- `_element_geometry` acepta "all OR any" nodo del elemento en la geometría
  (patrón `geometry_method`); el manual no especifica el método.
