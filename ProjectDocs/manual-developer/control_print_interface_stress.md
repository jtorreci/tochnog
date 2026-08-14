# control_print_interface_stress

## Implementación

- **Output**: `print_interface_stress()` in
  `print_interface_stress.cc` (new file, added to the makefile as
  `PRINT_IFACE_STRESS_OBJ`). Invoked from the control loop in `top.cc`
  (same pattern as `print_frd`).
- **Keywords** (data_class CONTROL) registered in `database.cc`:
  - `control_print_interface_stress` (INTEGER, length 1).
  - `control_print_interface_stress_2d_coordinates` (DOUBLE, length 4,
    required control_print_interface_stress) — cut line xstart ystart
    xend yend.
  - `control_print_interface_stress_3d_geometry` (INTEGER, length 2),
    `_3d_order` (INTEGER, length 1) — registered but not used.
- **New enums**: `CONTROL_PRINT_INTERFACE_STRESS`,
  `_2D_COORDINATES`, `_3D_GEOMETRY`, `_3D_ORDER` in `tochnog.h` /
  `tochnog-mod.h`, between `CONTROL_PRINT_HISTORY_SMOOTH` and
  `CONTROL_PRINT_MATLAB`.

## Algoritmo (2D)

- Recorre los elementos de interfaz (grupo con `GROUP_INTERFACE`).
- Lee el strain normal acumulado `ELEMENT_INTERFACE_STRAIN_NORMAL` y la
  rigidez `kn` del grupo; `sign = kn * strain_normal`.
- Para cada nodo del elemento, proyecta la coordenada sobre el vector del
  corte: `proj = ((x-xstart)*dx + (y-ystart)*dy)/len`.
- Escribe `proj sign sigt` por nodo en `interface_stress.<index>`
  (append por paso, serie temporal).

## Detalles

- `sigt` (tensión tangencial) se reporta 0: la fuerza tangencial
  acumulada no se almacena por elemento. Pendiente de implementar.
- El archivo crece por paso (append), como los otros prints.

## Pendiente

- 3D (`_3d_geometry`, `_3d_order`) registrados pero sin implementar.
- `sigt` tangencial no se calcula.
- El strain normal se lee de VERSION_NORMAL (el del paso anterior); el
  valor mostrado es el del último paso completado.
