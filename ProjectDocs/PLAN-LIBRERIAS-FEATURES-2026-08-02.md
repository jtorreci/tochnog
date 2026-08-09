# Plan: librerías de solver optimizadas + features de Tochnog Professional

Fecha: 2026-08-02
Rama: documentation-improvement

## 1. Librerías de solución de ecuaciones

### Estado actual (verificado)

El binario enlaza contra:
- `liblapack.so.3` / `libblas.so.3` del sistema → apuntan a **OpenBLAS 0.3.26** (via alternatives)
- **SuperLU 4.3** local: `external-downloads/superlu-4.3/lib/libsuperlu_4.3.a`

Solvers implementados en el código (`so.cc`, `so_bicg.cc`, `so_suplu.c`):

| Solver | Código | Estado |
|---|---|---|
| `-matrix_iterative_bicg` | `so_bicg.cc` | OK (por defecto) |
| `-matrix_lapack` (band) | `so.cc` | OK, usa LAPACK |
| `-matrix_superlu` (directo) | `so_suplu.c` | OK, SuperLU 4.3 |
| `-matrix_superlu_dist` | `so_suplu.c` | Compila, requiere MPI+SuperLU_DIST |
| `-matrix_superlu_mt` | `so_suplu.c` | Compila, requiere SuperLU_MT |
| `-matrix_iterative_petsc` | `so.cc` | NO compila (PETSC_USE=0, no instalado) |

### Hallazgos del sistema

- **OpenBLAS 0.3.26** ya activo (LAPACK/BLAS del sistema).
- **SuperLU 6.0.1** disponible en apt (`libsuperlu-dev`) y por GitHub — API **drop-in compatible** con la 4.3 que usa `so_suplu.c`.
- **SuperLU_DIST 8.2** instalado pero con API distinta (requiere reescribir `so_suplu.c`).
- **MUMPS 5.6** instalado con headers dev — no hay integración en el código.
- **ScaLAPACK 2.2, CombBLAS 2.0** presentes, sin uso.
- **PETSc**: NO instalado; el código lo soporta pero desactivado.

### Acción realizada

- Actualizado `scripts/build_safe.sh` para auto-detectar SuperLU (sistema → 6.0.1 local → 4.3).
- Añadida regla `so_suplu.o` al makefile (antes se compilaba a mano).
- Añadido `scripts/install_solvers.sh` (requiere sudo) para instalar `libsuperlu-dev` (SuperLU 6) y OpenBLAS.
- Compilado SuperLU 6.0.1 localmente y **verificado**: todos los tests de referencia pasan
  (hypo1-4, examp7/15, refine4, ho_othr2, examp22, truss1, elasti1, spring1, wave1).
- `examp9` falla igual (bad_alloc) con SuperLU 4.3 y 6.0.1 → problema del test (límite de memoria), no de la librería.

### Siguiente paso recomendado

Correr `sudo ./scripts/install_solvers.sh` para instalar SuperLU 6 del sistema (más limpio que el local).
Opcional a futuro: integrar MUMPS 5.6 o SuperLU_DIST 8 (ambos más modernos, requieren adaptar `so_suplu.c`).

---

## 2. Features de Tochnog Professional que faltan

Fuente: changelog oficial de tochnogprofessional.nl (captura web.archive.org 20240114), verificado contra `database.cc`.

### Leyenda
- **FALTA** = no existe la keyword en la implementación.
- Difícil ~ área del código donde iría (estimado).

### 2.1 Modelos de material (mayor dificultad)

| Feature | Keyword | Área |
|---|---|---|
| Hipoplasticidad Masin con OCR | `group_materi_plasti_hypo_masin` | `plasti.cc`, similar a hypo_wolfersdorff |
| Permeabilidad dependiente de tensión vertical | `group_groundflow_permeability_vertical_stress` | `groundda.cc` / flujo |
| Tension directa en plano normal | `group_materi_plasti_tension_direct_normal` (+ `_automatic`) | `plasti.cc` |
| Mohr-Coul en plano normal | `group_materi_plasti_mohr_coul_direct_normal` (+ `_automatic`) | `plasti.cc` |
| Gap de interfaz | `group_interface_gap` (cambio de default) | `interface` / `contact.cc` |

### 2.2 Comandos de control (dificultad media)

| Feature | Keyword | Área |
|---|---|---|
| Reset de valor de dof | `control_reset_value_dof` | `data.cc` / `top.cc` |
| Movimiento de malla | `control_mesh_move` | `mesh.cc` / `adjust.cc` |
| Cambiar coordenadas | `control_mesh_switch` | `mesh.cc` |
| Generar interfaz por geometría (2d/3d) | `control_mesh_generate_interface_geometry` | `generate.cc` |
| Límite de fuerza por reducción de velocidad | `bounda_time_until_*` | `bounda.cc` |
| Activación periódica de bounda_time | `bounda_time_on_off` | `bounda.cc` |
| Factor multiplicador en prints | `control_print_history_factor`, `control_print_data_versus_data_factor` | `print*.cc` |
| Print suavizado de historia | `control_print_history_smooth` | `print_hi.cc` |
| Desplazamiento relativo de material | `materi_displacement_relative` (i.c.w. `materi_velocity_integrated`) | `materi.cc` / `elem.cc` |
| Apply/no-apply de change_dataitem | `change_dataitem_apply` | `change.cc` |
| Check de uso de datos | `check_used` | `check.cc` |
| Iteraciones plásticas máximas | `materi_plasti_maximum_iterations` | `plasti.cc` |
| Factor de presión de poro | `groundflow_pressure_factor` | `groundda.cc` |

### 2.3 Post-proceso y utilidades (dificultad baja-media)

| Feature | Keyword | Área |
|---|---|---|
| Seguir partícula material | `post_point_move` | `print*.cc` |
| Fuerza puntual en el espacio | `force_point` | `force.cc` |
| Diagrama de asiento con tensiones | `strain_settlement_diagram*` | `print*.cc` |
| Factor de borde parabólico | `bounda_factor_parabolic_x` | `bounda.cc` |
| Slide axisimétrico | `slide_axisymmetric` | `slide.cc` |
| Fuerza de borde multi-lineal | `force_edge_multi_linear_factor_x` | `force.cc` |

### 2.4 Correcciones de bugs que también faltan en nuestra copia

| Bug (changelog) | ¿Existe la feature? |
|---|---|
| `control_reset_dof*` i.c.w. `geometry_method -any` | FALTA (tenemos `control_unknown_reset`, distinto) |
| `print_group_data` valores erróneos | FALTA (no existe `print_group_data`) |
| `control_print_gid_smooth_dof` con malla cambiante | FALTA |

---

## 5. Estado de implementación (actualizado 2026-08-02 tarde)

### Implementadas y verificadas (9)
Fase 1: `check_used`, `materi_plasti_maximum_iterations`,
`control_print_history_factor`, `control_print_data_versus_data_factor`,
`force_point`, `post_point_move`.
Fase 2: `control_mesh_switch`, `bounda_time_on_off`, `bounda_time_until_force`.

### Falta implementar (20, del changelog profesional)

**Materiales / mecánica de suelos:**
- `group_materi_plasti_hypo_masin` (+ OCR) — Fase 3
- `group_materi_plasti_tension_direct_normal` (+ `_automatic`)
- `group_materi_plasti_mohr_coul_direct_normal` (+ `_automatic`)
- `group_interface_materi_plasti_tension_direct` (re-introducida)
- `group_groundflow_permeability_vertical_stress`
- `groundflow_pressure_factor`
- `materi_displacement_relative`

**Controles / malla:**
- `control_mesh_move`
- `control_mesh_generate_interface_geometry` (2d/3d)
- `control_reset_value_dof`
- `control_print_history_smooth`
- `change_dataitem_apply`

**Post-proceso / utilidades:**
- `strain_settlement_diagram*`
- `force_edge_multi_linear_factor_x`
- `bounda_factor_parabolic_x`
- `slide_axisymmetric`
- `print_group_data`
- `bounda_dof ... -geometry_list ... -veln`

### Nota sobre "input"
No existe una keyword `input` como data record en el changelog profesional ni en el
manual 2011. "input" se refiere al **archivo de entrada** (`tn.dat`) y al mecanismo
`input_runtime` (leer registros en cada paso), que **ya está implementado** en
`input.cc:1221`. No hay feature "input" pendiente.

---

## 3. Priorización propuesta

### Fase 1 (rápido, alto valor, ~cada una 0.5-2 días)
1. `check_used` — verificación de que todos los datos del input se usan (auditoría).
2. `materi_plasti_maximum_iterations` — control de convergencia plástica.
3. `control_print_history_factor` / `control_print_data_versus_data_factor` — multiplicadores en salida.
4. `force_point` — fuerza puntual fuera de nodo.
5. `post_point_move` — seguimiento de partícula.

### Fase 2 (medio, ~2-5 días c/u)
6. `control_mesh_switch` — rotación de ejes para mallar. [x]
7. `control_mesh_move` — mueve la malla con desplazamiento lineal. [x]
8. `bounda_time_on_off` / `bounda_time_until_*` — control temporal de condiciones de contorno. [x]
9. `control_mesh_mirror` — refleja la malla. [x]
10. `control_mesh_copy` — copia la malla desplazada. [x]
11. `control_mesh_rotate` — rota 2D→3D. [x] DESBLOQUEADO (2026-08-04): implementados los elementos 3D. quad4→hex8, tria3→prism6. Con `number_of_integration_points` en initia.
12. `control_mesh_delete_element`, `control_mesh_keep_element`, `control_mesh_keep_element_group`, `control_mesh_change_element_group`. [x]
13. `control_mesh_keep_node`, `control_mesh_rotate_angle` (rotación 2D plana). [x]
14. **Elementos 3D**: TET4/TET10 (ya estaban), HEX8 (fórmula general polynom.cc), PRISM6 (implementado 2026-08-04). `control_mesh_rotate` usa estos.
15. `control_mesh_extrude` — extruye 2D→3D por capas (tria3→prism6, quad4→hex8). [x] Implementado 2026-08-04 en `mesh.cc` (mesh_extrude) + `extrude.cc`.
16. `control_mesh_remove` — borra elementos por método (method1: dentro de otros grupos). [x] Implementado 2026-08-04 en delete.cc (mesh_remove).
17. `control_mesh_convert` — convierte interfaces (requiere concepto de "interface", pendiente).
18. `groundflow_pressure_factor`.

**P2 COMPLETA (2026-08-04)**: switch, move, mirror, copy, delete_element, keep_element, keep_element_group, change_element_group, keep_node, rotate (2D→3D), rotate_angle, extrude, remove. Solo `control_mesh_convert` pendiente (depende de interfaces).

### Fase 3 (alto esfuerzo, ~1-2 semanas c/u)
11. `group_materi_plasti_hypo_masin` — modelo hipoplástico de Masin con OCR (reusar estructura de hypo_wolfersdorff).
12. `group_groundflow_permeability_vertical_stress`.
13. Integrar MUMPS 5.6 o SuperLU_DIST 8 como solver directo paralelo.

---

## 4. Nota de fiabilidad

La comparación "manual vs implementación" por extracción de texto del docx produce mucho ruido
(variables matemáticas, tokens partidos). La fuente fiable para features reales es el **changelog
de tochnogprofessional.nl** (archivado en `external-downloads/professional/changes-site-archive.txt`),
que lista features concretas introducidas tras la versión 2014. El plan de la sección 3 se basa en esa
lista, no en el docx.

---

## 6. PLAN DE TRABAJO COMPLETO (2026-08-03)

### Fuente definitiva
El **índice (outline) del manual profesional 2024** (419 págs, descargado de archive.org:
`tochnogprofessional.nl/manuals/user/user.pdf`, texto en `external-downloads/professional/UserManual-professional.txt`).
El índice es limpio: cada entrada de la sección 6 "data records" es una keyword real.

### Dimensionamiento
- **~719 keywords de data records del manual 2024 que no están en `database.cc`**.
- Inventario completo: `ProjectDocs/inventario-features-faltantes-2024.txt`.
- Las features se implementan 1 a 1: (1) leer del manual 2024 qué hace, (2) implementar en el
  código, (3) crear un test que la verifique.
- Estimación: ~1-3 features/día según complejidad → **varios meses de trabajo** en total.
  Se prioriza por valor práctico.

### Metodología (por feature)
1. Leer la sección del manual 2024 (`/tmp/user_manual.txt` o `UserManual-professional.txt`).
2. Añadir la keyword al enum (`tochnog.h` + `tochnog-mod.h`), registrarla en `database.cc`.
3. Implementar la lógica en el archivo correspondiente.
4. Crear un test `.dat` mínimo en `validation-suite/test-2014/` que ejercite la feature.
5. Verificar: build limpio (`./scripts/build_safe.sh --clean`), test nuevo PASS,
   y regresión (hypo1-4 + tests de referencia).

### Prioridades

#### P0 — Importación de mallas y archivos (valor alto, tu necesidad directa)
- [x] `include filename` — incluir archivos en el data part (implementado y verificado 2026-08-03).
- [x] `input_gmsh` — importar malla de gmsh (implementado 2026-08-03, lee `tochnog_in.msh` formato 2.2, nodos y elementos lineales/cuadráticos).
- [x] `input_abaqus` — importar malla de abaqus (implementado 2026-08-03, lee `abaqus.inp` y genera `tochnog_abaqus.dat` para `include`).
- [x] `input_abaqus_continue`, `input_abaqus_name`, `input_abaqus_group`, `input_abaqus_set` — sub-opciones (implementadas 2026-08-03: continue dispara la generación, name filtra tipos, group escribe material→group_*, set filtra elementos).
- [ ] `input_abaqus_mesh` — sub-opción pendiente (controla si se escriben timesteps/prints).
- [x] `input_feflow_mesh` — importar malla de FEFLOW (implementado 2026-08-03, lee `feflow.fem` ASCII con secciones `coordinates`/`elements`).
- [x] `input_feflow_fem` — sub-opción (implementado 2026-08-03: -yes usa `.fem`, -no usa `.dac`).
- [ ] `input_feflow_mesh_hydraulic_head` — sub-opción pendiente (requiere leer resultados de presión).

#### P1 — Chequeos y diagnóstico (COMPLETA 2026-08-04)
- [x] `check_target` — suprime el fallo de target (registra nota en tn.log).
- [x] `check_element_node` — detecta elementos con nodos duplicados.
- [x] `check_nan` — chequea NaN en node_dof al final de cada paso.
- [x] `check_solver eps` — avisa si términos de la diagonal del solver son < eps.
- [x] `check_element_shape factor` — avisa si elementos muy distorsionados.
- [x] `check_memory` — reporta pico de uso de RAM.
- [x] `check_memory_usage` — guarda pico de memoria en `check_memory_usage_result`.
- [x] `check_data` — verifica integridad de la base de datos (items requeridos presentes).
- [x] `check_error` — suprime mensajes de error (via `pri`).
- [x] `check_warning` — suprime mensajes de warning (via `pri`).
- [ ] `control_check_data`.

#### P2 — Control de malla (medio)
- [ ] `control_mesh_move`, `control_mesh_mirror`, `control_mesh_copy`,
      `control_mesh_rotate` (+`_angle`), `control_mesh_remove` (+sub-items),
      `control_mesh_keep_element`/`keep_element_group`/`keep_geometry`/`keep_node`,
      `control_mesh_change_element_group`, `control_mesh_duplicate_element_group`,
      `control_mesh_convert*`, `control_mesh_cut_geometry`, `control_mesh_delete_element`.

#### P3 — Condiciones de contorno (medio)
- [x] `bounda_constant` — mantiene dofs prescritos constantes. (2026-08-04)
- [x] `bounda_time_increment` — bounda_time con solo cargas e incremento fijo. (2026-08-04)
- [x] `bounda_time_offset` — offset de tiempo para bounda_time_increment. (2026-08-04)
- [x] `bounda_factor` — factor lineal por coordenada (a0+a1x+...). (2026-08-04)
- [x] `bounda_factor_parabolic_x` — factor cuadrático en x (a0+a1x+a2x²). (2026-08-04)
- [x] `bounda_found` — flag de impresión (indica si el bounda se usó). (2026-08-04, registro)
- [x] `bounda_time_units` — conversión de unidades de tiempo/longitud en bounda_time. (2026-08-04)
- [x] `bounda_geometry_method` — tipo de nodo (NODE/NODE_START_REFINED) para geometría. (2026-08-04)
- [x] `bounda_dof` — prescribe dofs (Dirichlet) por rango o geometría. Alias de `bounda_unknown`. (2026-08-04)
- [x] `bounda_alternate` — omite un bounda_dof rotativo entre iteraciones. (2026-08-04)
- [x] `bounda_normal` — nodos deslizan en un plano (velocidad normal anulada por proyección). (2026-08-04)
- [x] `bounda_water` — presión de poro desde la columna de agua (dens*g*Δz al nivel freático). (2026-08-04)
- [x] `bounda_dof_radial` — velocidad prescrita radial a un punto. (2026-08-04)
- [x] `bounda_dof_cylindrical` — velocidad prescrita cilíndrica a una línea. (2026-08-04)
- [ ] `bounda_baseline_correction` (+`_parameters`, requiere `bounda_time_smc` — señales sísmicas SMC, de nicho).

**P3 COMPLETA en lo implementable (2026-08-04)**: 14 features (constant, time_increment, time_offset, factor, factor_parabolic_x, found, time_units, geometry_method, dof, alternate, normal, water, dof_radial, dof_cylindrical). Solo `bounda_baseline_correction` pendiente (de nicho, requiere SMC).

#### P4 — Materiales geotécnicos (alto esfuerzo)
- [x] `groundflow_pressure_factor` — multiplicador de la presión de poro al calcular el esfuerzo total. (2026-08-04)
- [ ] `group_materi_plasti_hypo_masin` (+ `_clay`, `_clay_advanced_parameters`,
      `_clay_ocr`, `_ocr`, `_structure`), `control_materi_plasti_hypo_masin_ocr_apply`.
- [ ] `group_materi_plasti_tension_direct_normal` (+ `_automatic`) — requiere `group_materi_plasti_tension_direct` (no existe).
      `group_materi_plasti_mohr_coul_direct_normal` (+ `_automatic`) — requiere `group_interface_materi_plasti_mohr_coul_direct`.
- [ ] `group_groundflow_permeability_vertical_stress`, `groundflow_pressure_factor`.
- [ ] `group_materi_damage_mazars`, `group_materi_expansion_linear`,
      `group_materi_expansion_volume`.
- [ ] Modelos hiperelásticos: `group_materi_hyper_besseling`, `_blatz_ko`,
      `_mooney_rivlin`, `_neohookean`, `_reduced_polynomial`, `_volumetric_*`.
- [ ] Viscoelástico/viscoplástico: `group_materi_maxwell_chain`,
      `group_materi_plasti_visco_exponential*`, `_visco_power*`.

#### P5 — Post-proceso y salida (medio)
- [ ] `control_print_history_smooth`, `control_print_gid_*` (varios),
      `control_print_vtk_*`, `control_print_gmsh_*`, `control_print_frd_*`,
      `control_print_materi_stress_force`, `control_print_interface_stress*`.
- [ ] `print_group_data`, `strain_settlement_diagram*`, `post_point_move` (hecho).
- [ ] `force_edge_multi_linear_factor_x`, `force_edge_projected*`, `force_volume*`,
      `force_gravity_geometry`.

#### P6 — Groundflow, contacto, miscelánea
- [ ] `groundflow_pressure_factor`, `groundflow_seepage_*`, `groundflow_phreatic_*`,
      `groundflow_flux_edge_normal*`, `groundflow_total_pressure_limit`.
- [ ] `groundflow_pressure_atmospheric` — relacionado con `bounda_water`: el
      clamp de `groundflow_phreatic_coord()` limita la presión estática a este
      umbral (default 0 → solo succión). Documentado en manual-developer/bounda_water.md.
- [ ] Contacto: `contact_apply`, `contact_heat_generation`, `contact_penalty_*`,
      `contact_plasti_friction`, `contact_target_*`.
- [ ] `control_reset_dof*`, `control_reset_value_*`, `change_dataitem_apply`,
      `change_dataitem_geometry`, `control_data_*`, `control_distribute*`.
- [ ] `control_mesh_generate_interface_geometry`, `slide_axisymmetric`,
      `materi_displacement_relative`, `group_interface_materi_plasti_tension_direct`.

### Fichero de control
- Inventario completo: `ProjectDocs/inventario-features-faltantes-2024.txt` (718 líneas).
- Marcar con `[x]` cada feature al implementarla y verificar su test.
- Al final de cada sesión: guardar el progreso en memoria (Engram) con el plan actualizado.

---

## 7. Ideas a explorar (proyectos futuros, no planificadas)

### 7.1 Exportación a Paraview
- Via `control_print_vtk_*` (features del manual profesional P5).
- Alto valor para post-proceso: visualización de resultados en Paraview.

### 7.2 Librería pytochnog (wrapper Python)
Objetivos:
- **Preproceso**: generar mallados y archivos `.dat` desde Python.
- **Lanzar cálculos**: ejecutar tochnog desde Python, incluso en bucle
  (parametrización) o en paralelo (multiples cálculos simultáneos).
- **Explotación de resultados**: leer las bases de datos de resultados
  (tn.dvd, tn.log, .his), para gráficas, extremos, combinaciones de carga,
  etc.
- Arquitectura sugerida: subprocess para lanzar tochnog + parsing de los
  archivos de salida (o bindings C++ si se decide integrar la librería).
- Relacionado: la feature `include` (ya implementada) permite parametrizar
  inputs desde Python generando archivos incluidos.

### 7.3 Desbloqueo de features bloqueadas de P2
- `control_mesh_rotate`: DESBLOQUEADO (2026-08-04). Implementados los elementos
  3D: HEX8 (vía la fórmula general tensorial en polynom.cc, `npol=2`) y PRISM6
  (shape functions de prisma triangular explícitas en polynom.cc). `mesh_rotate_3d`
  en mesh.cc convierte quad4→hex8 y tria3→prism6 (rota alrededor del eje y).
  Requiere `number_of_integration_points` (≥8 para hex8, ≥6 para prism6) en la
  inicialización. Pendiente: n>1 segmentos rotacionales (multi-capa).
- `control_mesh_keep_node`: DESBLOQUEADO (2026-08-04) usando `delete_node()` +
  `db_delete_index()` que ya existían en delete.cc.

---

## 8. Elementos de Tochnog Professional (inventario y estado)

### 8.1 Estado actual (2026-08-04)

| Elemento | Estado | Notas |
|---|---|---|
| `-bar2/3/4` | ✓ | Truss/barras |
| `-tria3/6` | ✓ | 2D |
| `-quad4/9/16` | ✓ | 2D |
| `-tet4/10` | ✓ | 3D (shape functions explícitas) |
| `-prism6` | ✓ | 3D (nuevo, 2026-08-04) |
| `-hex8` | ✓ | 3D (fórmula general tensorial) |
| `-hex27`, `-hex64` | enum+name | Referenciados por la fórmula general (npol=3/4) pero no verificados |
| `-quad6/8` | **FALTA** | 2D cuadrático/transición |
| `-prism12/15/18` | **FALTA** | 3D prismas de orden superior |
| `-hex18/20` | **FALTA** | 3D hexaedros de orden superior |

### 8.2 Elementos de orden superior que faltan

- **`-quad6`** (6 nodos: quad con 2 lados de 3 nodos, para transiciones
  tria6/quad9) — shape functions de transición.
- **`-quad8`** (8 nodos, serendipity) — fórmula tensorial mixta.
- **`-prism12/15/18`** — prismas con lados cuadráticos (extensión de PRISM6).
- **`-hex18/20`** — hexaedros serendipity/transición (extensión de HEX8).

Implementación: añadir shape functions + puntos de integración en `polynom.cc`
(mismo patrón que PRISM6/la fórmula general). La infraestructura (enum,
database.cc, elem.cc genérico) ya está lista.

### 8.3 Membrane y plate (proyecto futuro — anotado)

- **`membrane`**: elementos 2D SIN rigidez a flexión (solo tensión en el plano).
  Factible: se apoya en QUAD4/TRIA3 existentes. Requiere formulación de membrana
  (estado plano de tensiones, sin dofs rotacionales).
- **`plate`**: elementos 2D CON rigidez a flexión (desplazamientos transversales
  + rotaciones). MÁS COMPLEJO: añade grados de libertad rotacionales (Reissner-
  Mindlin o Kirchhoff), lo que toca el ensamblaje y el solver.
- Recomendación: membrane primero (apoyado en elementos planos), plate después.
- El usuario NO plantea shell por ahora (complejidad de definición y acoplamiento).

### 8.4 Beam

- **`-beam`**: existe en el enum (`BEAM`) y hay `beam.cc` (beam 2D con 3
  dofs/nodo: NDOF=3, NNOL=2, NDIM=2, más `beam_3d`). PERO no está registrado
  como elemento de malla en database.cc con data_length. El manual 2024 lista
  solo `-bar2/3/4` en 1D — beam es una feature de la versión GNU original.
  Verificar cómo se activa (si es por otro mecanismo) antes de decidir.
