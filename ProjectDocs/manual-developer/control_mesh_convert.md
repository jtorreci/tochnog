# control_mesh_convert

## Implementación

- **Output**: `interface_convert()` in `interface.cc`. Invoked from
  `step_start()` in `top.cc` (after `generate_spring`), so it runs before
  the element assembly of each step.
- **Keywords** (data_class CONTROL) registered in `database.cc`:
  - `control_mesh_convert` (INTEGER, length 1) — switch.
  - `control_mesh_convert_element_group` (INTEGER, DATA_ITEM_SIZE,
    required control_mesh_convert) — groups on one side of the interface.
- **New enums**: `CONTROL_MESH_CONVERT`, `CONTROL_MESH_CONVERT_ELEMENT_GROUP`
  in `tochnog.h` / `tochnog-mod.h` (kept in sync), between
  `CONTROL_MESH_GENERATE_TRUSS_BEAM_MACRO` and `CONTROL_MESH_KEEP_ELEMENT`.

## Algoritmo (bar2 -> quad4, bar3 -> quad6)

Para cada elemento `-bar2`/`-bar3` cuyo grupo tiene `group_interface -yes`:

1. Leer los nodos {a,b} (lado 1 de la interfaz) — 2 nodos para el bar2,
   3 para el bar3 (2D cuadrático).
2. Calcular la normal perpendicular a la línea a-b (2D).
3. Crear ns1 nodos nuevos {a',b',...} = copias de los del lado 1
   desplazados `shift` en la normal (`shift = 0.01 * length(a-b)`;
   el tangent se normaliza antes, así que shift queda 0.01 constante —
   la posición exacta de los duplicados es físicamente irrelevante para
   la penalización). Copiar NODE, NODE_START_REFINED, NODE_DOF,
   NODE_DOF_START_REFINED (patrón de generate.cc).
4. Reescribir el elemento como `-quad4` {a,b,a',b'} o `-quad6`
   {a,b,c,a',b',c'} (3+3 nodos).
5. Reconectar vecinos: para cada elemento que comparte TODOS los nodos
   del lado 1 y NO está en `control_mesh_convert_element_group`,
   reemplazar esos nodos por sus duplicados. El test del centroide (dot
   con la normal > 0) selecciona solo el bloque del lado +normal: el
   bloque del lado -normal conserva los nodos originales.
6. `mesh_has_changed(VERSION_NORMAL)` si se convirtió al menos un elemento.

NOTA (verificada contra el binario del Professional 25-10-2023, .dbs de
interface_bar2_quad4/interface_bar3_quad8): el Professional crea los
duplicados CON LAS MISMAS coordenadas (interfaz de espesor cero) y los
asigna al bloque del lado -normal (el bloque +normal conserva los nodos
originales); el GNU los desplaza 0.01 y los asigna al bloque del lado
+normal. La física de penalización por pares es la misma (los tests del
corpus pasan con rc=0 y valores idénticos al Pro en ~1e-7); solo cambia
el layout de numeración del .dbs.

## Detalles

- `interface_element()` admite `-bar2` no convertido (degenerado: nodo 0 =
  lado 1, nodo 1 = lado 2) para no abortar si un bar2 de interfaz llega al
  ensamblaje sin conversión. El `-quad6` (3+3, pesos Lobatto 1/6,4/6,1/6,
  historias por IP) ya existía en `interface_element()`.
- La conversión es idempotente: en pasos posteriores el elemento ya es
  `-quad4`/`-quad6`/`-prism6`/`-hex8`/`-hex18`, no se re-convierte.
- El bar3 2D necesita que los sólidos adyacentes sean cuadráticos: los
  `-quad8` del input se convierten a `-quad9` ANTES (mesh_convert_quad8,
  ver manual quad8.md) — el orden en step_start es convert quad8 ->
  convert hex20 -> extrude -> interface_convert.
- En 3D un `-bar3` de interfaz se trataría como línea (rama bar2 3D de la
  normal), pero no hay caso de corpus; el elemento cuadrático facial 3D es
  la familia `-quad8` interface -> `-hex18` (implementado, ver
  hex18.md).

## Pendiente

- `control_mesh_convert_quad9_quad6` y `_tria6_tria3` no están
  implementados.
