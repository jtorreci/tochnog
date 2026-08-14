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

## Algoritmo (bar2 -> quad4)

Para cada elemento `-bar2` cuyo grupo tiene `group_interface -yes`:

1. Leer los nodos {a,b} (lado 1 de la interfaz).
2. Calcular la normal perpendicular a la línea a-b (2D).
3. Crear 2 nodos nuevos {a',b'} = copias de {a,b} desplazados `shift` en
   la normal (`shift = 0.01 * length(a-b)`). Copiar NODE,
   NODE_START_REFINED, NODE_DOF, NODE_DOF_START_REFINED (patrón de
   generate.cc).
4. Reescribir el elemento como `-quad4` {a,b,a',b'}.
5. Reconectar vecinos: para cada elemento que comparte {a,b} y NO está en
   `control_mesh_convert_element_group`, reemplazar {a,b} por {a',b'} en
   su conectividad.
6. `mesh_has_changed(VERSION_NORMAL)` si se convirtió al menos un elemento.

## Detalles

- `interface_element()` admite `-bar2` no convertido (degenerado: nodo 0 =
  lado 1, nodo 1 = lado 2) para no abortar si un bar2 de interfaz llega al
  ensamblaje sin conversión.
- La conversión es idempotente: en pasos posteriores el elemento ya es
  `-quad4`, no se re-convierte.

## Pendiente

- Solo `-bar2` -> `-quad4` (2D). Los casos 3D (`-tria3` -> `-prism6`,
  `-quad4` -> `-hex8`) y `-bar3`/`-quad8`/`-quad9` -> `-quad6` no están
  implementados.
- `control_mesh_convert_quad9_quad6` y `_tria6_tria3` no están
  implementados.
