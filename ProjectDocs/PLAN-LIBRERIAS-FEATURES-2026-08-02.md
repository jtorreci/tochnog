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

## 3. Priorización propuesta

### Fase 1 (rápido, alto valor, ~cada una 0.5-2 días)
1. `check_used` — verificación de que todos los datos del input se usan (auditoría).
2. `materi_plasti_maximum_iterations` — control de convergencia plástica.
3. `control_print_history_factor` / `control_print_data_versus_data_factor` — multiplicadores en salida.
4. `force_point` — fuerza puntual fuera de nodo.
5. `post_point_move` — seguimiento de partícula.

### Fase 2 (medio, ~2-5 días c/u)
6. `control_mesh_switch` — rotación de ejes para mallar.
7. `control_mesh_move` — mallas oblicuas.
8. `bounda_time_on_off` / `bounda_time_until_*` — control temporal de condiciones de contorno.
9. `group_materi_plasti_tension_direct_normal` + `mohr_coul_direct_normal` (con/auto).
10. `groundflow_pressure_factor`.

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
