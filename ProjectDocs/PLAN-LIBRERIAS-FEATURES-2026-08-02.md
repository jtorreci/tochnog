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
- [ ] `control_check_data`. — **NOTA 2026-08-12**: verificar si existe el keyword `check_data`; `control_check_data` no está en database.cc.

#### P2 — Control de malla (medio)
- [x] `control_mesh_move`, `control_mesh_mirror`, `control_mesh_copy`,
      `control_mesh_rotate` (+`_angle`), `control_mesh_remove` (+sub-items),
      `control_mesh_keep_element`/`keep_element_group`/`keep_geometry`/`keep_node`,
      `control_mesh_change_element_group`, `control_mesh_delete_element`.
      **NOTA 2026-08-12**: presentes en database.cc; el estado exacto
      (implementación y tests) debe verificarse antes de marcarlas como
      completas. `control_mesh_duplicate_element_group`,
      `control_mesh_convert*`, `control_mesh_cut_geometry` pendientes.

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
- [x] **Familias de modelos YA IMPLEMENTADAS (verificadas y documentadas 2026-08-10, P4-A)**:
  - Hiperelásticos: `group_materi_hyper_besseling`, `_blatz_ko`, `_mooney_rivlin`,
    `_neohookean`, `_reduced_polynomial`, `_volumetric_linear`, `_volumetric_murnaghan`,
    `_volumetric_ogden`, `_volumetric_polynomial`, `_volumetric_simotaylor`,
    `_hyper_stiffness` — hyperela.cc, todos funcionales (tests blatz1, ho_mech1).
  - Viscoelástico: `group_materi_maxwell_chain` (viscoela.cc, test viscel1).
    `group_materi_maxwell_chain_nonlinear` REGISTRADO pero `visco_elastiticity_nonlinear()`
    (visconon.cc) es un STUB vacío — no funcional.
  - Viscosidad: `group_materi_viscosity`, `_viscosity_heatgeneration`,
    `_viscosity_user` (viscosit.cc, test viscos1; `user_viscosity()` es stub).
  - Daño: `group_materi_damage_mazars` + `materi_damage` (damage.cc, test damage1).
    **Bug CORREGIDO (2026-08-10)**: acumulación de `epseq` en damage.cc:92 usaba
    `epseq += epseq + x` (duplicante) en vez de `epseq += x`. Presente en el fuente
    GNU original de 2014. Test nuevo `damage2.dat` (biaxial no-equal) lo detecta
    (fixed dam=0.471973, buggy dam=0.789301).
  - Expansión: `group_materi_expansion_linear` (stress.cc:437) y
    `group_materi_expansion_volume` (materi.cc:166, solo densidad) — tests expans1/2.
  - Viscoplástico: `group_materi_plasti_visco_exponential`, `_visco_power`,
    `_visco_always` (stress.cc:100-115,737-755). Falta
    `..._visco_exponential_limit` (hardcode `EPS_VISCO=3`) y variantes
    `_name`/`_values`.
- [ ] `group_materi_plasti_hypo_masin` (+ `_clay`, `_clay_advanced_parameters`,
      `_clay_ocr`, `_ocr`, `_structure`), `control_materi_plasti_hypo_masin_ocr_apply`.
      **P4-B1 PARCIAL (2026-08-10)**: `group_materi_plasti_hypo_masin` (ley básica,
      5 params), `_structure`, `_ocr`, `control_*_ocr_apply` IMPLEMENTADOS via
      `masin.c` (port C del UMAT Fortran autorizado, validado a 5e-7) + dispatch
      en hypoplas.cc. Test hypomasin1.dat end-to-end (sigxx -340.7 vs -334.8,
      e 0.6333 vs 0.6663). Pendiente P4-B2: variante clay anisotrópica
      (alpha_G/alpha_E/alpha_nu/dirección), intergranular strain masin, visco,
      y las variantes strength-reduction/visco del UMAT (paquetes descargados).
      **P4-B2 (2026-08-10)**: `group_materi_plasti_hypo_masin_clay` (anisotrópico,
      5 params) + `_clay_advanced_parameters` (αG αf ay oc) +
      `_clay_avanced_direction` (diri) + `_clay_ocr` + `_clay_structure` +
      `control_*_clay_ocr_apply` + `group_materi_plasti_hypo_strain_intergranular_masin_clay`
      (R Ag ng mrat βr χ θ) IMPLEMENTADOS. Tests: hypomasin2.dat (αG=2: sigxx
      -429.6 vs -418.6 ref), hypomasin3.dat (intergranular: kernel OK, end-to-end
      con discrepancia ~20% por iteraciones del equilibrio — pendiente de fix).
      **P4-B3 (2026-08-10)**: `group_materi_plasti_hypo_masin_clay_visco_jm`
      (ocparam beta_deg ksi gama_deg Dref) IMPLEMENTADO via `masin_visco.c`
      (port C de umat_visco.f, Jerman-Masin 2020, validado: sig11 -82.5684
      idéntico al Fortran). Test hypomasin4.dat (sigxx -82.88 vs -82.57).
      **P4-B4 (2026-08-10)**: `group_materi_plasti_hypo_masin_clay_visco`
      (Dr Iv, Niemunis visco law del manual) IMPLEMENTADO desde la teoría
      (sin UMAT autorizado — el "Niemunis" de soilmodels es el HCA 2005, no
      relacionado). Kernel masin_niemunis_visco_umat con defaults derivados
      del clay (lambda=λ*, ee0=e0, pe0=p0, betaR=1), creep clamp por
      estabilidad. Test hypomasin5.dat (sigxx -104.1). Refactor interfaz:
      Dref movido a _clay_visco_jm para dejar _clay_visco fiel al manual.
- [ ] `group_materi_plasti_tension_direct_normal` (+ `_automatic`) — requiere `group_materi_plasti_tension_direct` (no existe).
      `group_materi_plasti_mohr_coul_direct_normal` (+ `_automatic`) — requiere `group_interface_materi_plasti_mohr_coul_direct`.
- [x] `group_groundflow_permeability_vertical_stress` — kp = a/(sigv/sig0)^b
      con clamp [min,max], combinada con group_groundflow_permeability.
      **IMPLEMENTADO y VALIDADO (2026-08-11)**: groundda.cc groundflow_data.
      Fix clave P4-E2b: el stress nodal se lee via db_dbl(NODE_DOF, inod,
      VERSION_NEW), no indexando new_unknowns manualmente; el h[] pasado es
      un indicador de base nodal, no shape functions (promedio sobre nodos).
      VALIDADO: sigv=-100 -> kp=0.1; sigv=-200 -> kp=0.05 (la ley exacta).
      Test groundperm_vs.dat (pres nodo medio 0.498756 vs 0.5, 0.25%).
- [x] `materi_plasti_hypo_*` variantes del kernel hipoplástico (lowangles, cohesion,
      intergranularstrain, pressuredependentvoidratio, wolfersdorff) — registradas,
      verificar lógica y tests (P4-A1). **VERIFICADO (2026-08-11, P4-A1)**: las 5
      variantes funcionan. Tests: hypo1 (wolfersdorff), hypo2/3/4
      (intergranularstrain), + 3 nuevos: hypo_cohesion.dat (c=5, sigyy=-3945.4),
      hypo_lowangles.dat (rval=2, powxi=2, sigyy=-4263.7), hypo_pdvr.dat
      (-yes, sigyy=-943.0). Documentado en
      manual-user/manual-developer/group_materi_plasti_hypo_wolfersdorff.md.

#### P4-E — Otros modelos de suelo (SoilModels.com) — METODOLOGÍA MASIN
Metodología establecida en P4-B1 (port fiel del UMAT Fortran autorizado → C puro
→ validación numérica a ~1e-6 contra el original vía driver de elemento único →
integración con keywords/dispatch/tests → documentación). Infraestructura
reutilizable: `masin.c` como plantilla de kernel, hypoplas.cc como patrón de
dispatch, drivers de validación en validation-suite/reference-masin/.

El usuario descarga el UMAT de SoilModels a petición (login gratuito). Orden
propuesto por demanda práctica:
- [x] `SANISAND` (Dafalias & Manzari, arena no-cohesiva, UMAT disponible) — P4-E1.
      **PORT FUNCIONAL (2026-08-10)**: `sanisand.c` (~1600 líneas C) integrado en
      tochnog (keyword `group_materi_plasti_sanisand`, 19 params,
      `materi_history_variables 36`). Test hyposanisand1.dat end-to-end pasa.
      **LIMITACIÓN**: valida ~8% vs el Fortran (no 1e-6 como Masin) por el
      substepping adaptativo elasto-plástico simplificado. Constitutivamente
      correcto (endurecimiento de arena densa).
      **P4-E1b (2026-08-10)**: afinado a ~3.5% (driver, paso 20: C=-2873 vs
      F=-2777). Fix: tol_f=1e-6 (era 1e-3, la tol de yield del Fortran es
      independiente de testing) + llamada duplicada de f_plas eliminada.
      **P4-E1d (2026-08-10)**: fix del intersect_DM — (1) el bucle interno de
      Newton acumulaba xi prematuramente (xip1 desde xi fijo, el Fortran usa
      xip1=xi+dxi); (2) la bisection devolvía el punto medio 0.5, el Fortran
      usa el xi del cruce del Newton (~0.018). Ambos bugs se compensaban:
      arreglar solo uno empeoraba el global. Con ambos corregidos: e coincide
      0.007% y pasos 1-11 al 0.1%; el paso 20 acumula ~5.4% (substepping
      plástico fino).
      **P4-E1f (2026-08-10)**: BUG DEFINITIVO del intersect encontrado — la
      bisection del Fortran ESTRECHA el intervalo (if fy05<0 -> y00=y05 else
      y11=y05), el port tenia y00/y11 fijos (convergia a 0.5). Porte la
      bisection correcta: paso 1 EXACTO al Fortran (sig11=-189.9448,
      e=0.69833, a11=0.206539 — identicos), e exacta 0.001% en todos los
      pasos, pasos 1-8 al 0.12%. Verificada la identidad de TODAS las
      funciones fisicas (sin errores de transcripcion). La divergencia
      residual (paso 20: a11 C=0.675 vs F=0.606) es SOLO el deviator.
      **P4-E1h (2026-08-10)**: diagnosticado y cerrado como VALIDADO. De/Gt
      identicos; la diferencia de LDeR viene de los gradientes (dependen de
      alpha). Prueba decisiva: recompilar el MISMO C con -O0 vs -O2 cambia
      sig11 del paso 1 en 1.1e-3, ~20x MAS que la diferencia C-vs-Fortran
      (4.7e-5). INDECIDIBLE cual es el correcto: el Fortran es otra
      compilacion con otro orden de operaciones, no una verdad. La diferencia
      late (6.3%) es el limite de reproducibilidad entre implementaciones
      numericamente equivalentes, amplificado por el substepping adaptativo
      (caotico). e exacta 0.001% y pasos 1-8 al 0.12% validan la fisica. No
      se persigue el bit-exact (espejismo). Documentado en ambos manuales.
- [ ] `PM4Sand` (Boulanger & Ziotopoulou, arenas licuefactables) — P4-E2.
      **SOLO ARTÍCULOS DISPONIBLES (2026-08-10)**: el UMAT NO está descargado
      (solo los papers de Boulanger & Ziotopoulou 2017 y relacionados). Para
      implementarlo con la metodología Masin hace falta el UMAT de referencia
      (SoilModels o UCD) para validar numéricamente. Sin él, implementar desde
      la teoría sería como el visco Dr/Iv (validación débil). Posponer hasta
      obtener el UMAT. LECCIÓN P4-E: la implementación de modelos constitutivos
      se ha revelado como retante hasta conseguir una CALIBRACIÓN CORRECTA
      (SANISAND requirió 6 iteraciones de debug: intersect, bisection,
      tolerancias); el orden de operaciones del compilador impone un límite de
      reproducibilidad (no bit-exact) entre implementaciones equivalentes. Para
      cada nuevo modelo: (1) obtener el UMAT autorizado, (2) driver de
      referencia primero, (3) validar contra resultado publicado o tercera
      implementación, no solo contra el UMAT.
- [ ] `Sand Hypoplasticity` (Gudehus/Bauer, wolfersdorff ya cubierto en hypo.c —
      solo validar) — P4-E3.
- [ ] `EMC` / otros (ISA, barodesy, viscohypoplasticity) según demanda — P4-E4.
- [ ] `hypo.c` refactor a C idiomático (ver P4-F) ANTES de portar más kernels
      f2c: cada nuevo modelo se portaría directo a C limpio.

#### P4 pendientes documentados para el futuro (2026-08-11)
- **Firma f2c de hypo_ (cosmético)**: `hypo_` conserva la declaración de
  parámetros después de la firma (estilo f2c años 90). Funciona y no bloquea;
  modernizarla (mover los parámetros a la firma) es de riesgo sin beneficio
  funcional. POSPUESTO.
- **P4-C (bloqueado)**: `group_materi_plasti_tension_direct_normal` y
  `group_materi_plasti_mohr_coul_direct_normal` requieren
  `group_materi_plasti_tension_direct` y elementos de interfaz
  (`group_interface_*`) que NO existen en esta base. Requiere evaluar primero
  si los elementos de interfaz son implementables. BLOQUEADO.
- **P4-E (modelos nuevos)**: PM4Sand (sin UMAT), Sand Hypoplasticity (solo
  validar wolfersdorff ya cubierto), EMC/otros (ISA, barodesy) según demanda.
  Requieren UMAT de referencia o trabajo de calibración extenso. POSPUESTO.

#### P4-F — Refactor de código adaptado de Fortran (hypo.c) — C idiomático
`hypo.c` es un port f2c→C puro (wolfersdorff, 1380 líneas, `static` locals,
notación de punteros f2c, `f2c.h` con tipos `integer`/`doublereal`). Objetivos:
- [x] Eliminar dependencia de `f2c.h` (tipos propios `long int`/`double`,
      funciones `hypo_*` con firma explícita). El binario ya NO enlaza libf2c
      (0 símbolos); `f2c.h` solo aportaba typedefs. **HECHO (2026-08-11)**: hypo.c
      convertido a C puro (doublereal→double, integer/logical/ftnlen→long int,
      macros TRUE_/FALSE_/min/max/abs definidas localmente), `f2c.h` ELIMINADO
      (de hypo.c y de tochnog.h/tochnog-mod.h). Bug clave: `abs` de f2c.h era una
      macro que funcionaba con double; sin ella, `abs(int)` de stdlib truncaba a
      int y rompía hypo1-3. Validado: hypo1-4 + regresión completa en verde.
- [x] Reemplazar `static` locals y paso por referencia estilo f2c por structs
      de estado por punto de integración (reentrante, sin estado global).
      **HECHO (2026-08-11, paso 2)**: los 46 `static` locals de hypo.c
      convertidos a automáticas (las 3 constantes c__9/c_b45/c__81 quedan como
      static const de solo lectura). Los warnings de "may be uninitialized"
      (tcohesion, c1, c2) son falsos positivos: copy_/power_ son funciones
      externas que inicializan. Validado: hypo1-4 + regresión completa
      (Masin, SANISAND, damage, groundflow) en verde. Esto hace hypo_ REENTRANTE
      (habilita OpenMP futuro). Nota: la regresión confirma que ninguna variable
      local dependía de la persistencia entre llamadas.
- [x] Eliminar macros `min`/`max` de f2c.h que rompen C++ estándar
      (tochnog.h:45 ya documenta el workaround). **HECHO (2026-08-11, paso 3)**:
      sin macros min/max en ningun .cc. math.cc (que definia su propia macro
      max) usa ternarios; hypo.c usa ternarios (no scalar_dmax/dmin porque el
      enlazado C de hypo.c no resuelve el mangling C++ de math.cc). El
      bloque #undef min/max obsoleto de tochnog.h fue removido en el paso 1.
- [ ] Beneficio esperado de velocidad: bajo para el propio wolfersdorff (el
      código numérico es el mismo); el gana está en legibilidad, reentrancia y
      en los kernels nuevos (masin.c ya es C limpio). La validación de
      regresión: hypo1-4 + hypomasin1.

#### P5 — Post-proceso y salida (medio)

#### P5-T — Exportación tabular + SQLite (post-proceso programático) [PLAN APROBADO 2026-08-11]
Objetivo: almacén SQLite + CSV opcional para post-proceso con código (pandas,
análisis estadístico, magnitudes derivadas, gráficas, esfuerzos sobre líneas).
Dependencias OPCIONALES en compilación (patrón SUPERLU/PETSC: `tn_sqlite.h`
con `SQLITE_USE`; si se compila sin soporte, al solicitar la feature se
advierte y se continúa con CSV). En esta máquina se instala `libsqlite3-dev`
(runtime ya presente).
- [x] **P5-T0**: infraestructura condicional — `tn_sqlite.h` (`SQLITE_USE`),
      makefile (`SQLITE_INCLUDE`/`SQLITE_LIB`), `sqlite.cc` (clase `SqliteDB`
      con RAII, esquema normalizado). Instalar `libsqlite3-dev`. (commit 81d0d5d)
- [x] **P5-T1**: `control_print_tabular` — CSV de un instante (`-last`) +
      escritura SQLite tabla `primary`. Reutiliza acceso `db_dbl(NODE_DOF,
      VERSION_PRINT)` + `dof_scal_vec_mat` (patrón print_vt.cc). (commit 1c677bb)
- [x] **P5-T2**: series temporales — CSV multi-incremento con columna `t` +
      SQLite por paso (patrón print_history). Desde el inicio. (commit 7904b4c)
- [x] **P5-T3**: magnitudes derivadas C++ (`template<int D> Tensor`, von Mises,
      Tresca, principales con `matrix_jacobi`) → tabla `derived`. Integración
      con `print_vtk` (`POINT_DATA` sin duplicar lógica).
- [x] **P5-T4**: `tools/postprocess.py` (pandas+sqlite3): estadísticas, gráficas
      tiempo, esfuerzos sobre `geometry_line`, variables de usuario → tabla
      `user`.

Esquema SQLite (normalizado por tablas, clave compuesta `(node,t)`):
```
primary(node, t, ux.., sigxx.., exx..)     -- variables fijas
derived(node, t, vmises, tresca, sig1..3)  -- magnitudes C++
user(node, t, energia, ...)                -- variables Python
meta(key, value)                           -- malla, unidades, convencion
```
Añadir magnitudes NO altera `primary` (tablas separadas por familia; JOIN por
(node,t) en SQL/pandas).


- [ ] `control_print_history_smooth`, `control_print_vtk_*`.
- [ ] ~~`control_print_materi_stress_force`~~ — **DEPENDIENTE (2026-08-12)**.
      Requiere `post_calcul_materi_stress_force` (integración de tensiones
      sobre cortes), que no existe en el GNU. Post-proceso numérico no
      trivial, pendiente de infraestructura.
- [ ] ~~`control_print_interface_stress*`~~ — **DEPENDIENTE (2026-08-12)**.
      Requiere elementos de interfaz/contacto (tensiones de interfaz), que
      el GNU no tiene (ni enums ni keywords). Implica reescribir parte del
      solver; fuera de alcance de post-proceso.
- [ ] ~~`control_print_gid_*` (varios)~~ — **DESCARTADO (2026-08-12)**.
      GiD es el único formato propietario contemplado, el `print_gid_6` del
      GNU es antiguo frente a los cambios recientes, y GiD puede importar
      formatos no nativos (VTK, Gmsh, CSV). Si hay que codificar un formato
      propietario, es preferible un formato estándar. La cobertura estándar
      actual (VTK, GMSH, FRD, CSV/SQLite) se considera suficiente.
- [x] `control_print_gmsh` + `control_print_gmsh_dummy` +
      `control_print_gmsh_element_data` + `control_print_gmsh_node_method`
      (familia GMSH, 2026-08-12).
- [x] `control_print_frd` + `control_print_frd_freecad` +
      `control_print_frd_prepomax` (familia FRD, 2026-08-12).
- [ ] `print_group_data`, `strain_settlement_diagram*`, `post_point_move` (hecho).
- [ ] `force_edge_multi_linear_factor_x`, `force_edge_projected*`, `force_volume*`,
      `force_gravity_geometry`.

#### P6 — Groundflow, contacto, miscelánea
- [ ] `groundflow_seepage_*`, `groundflow_phreatic_*`,
      `groundflow_flux_edge_normal*`, `groundflow_total_pressure_limit`.
- [ ] `groundflow_pressure_atmospheric` — relacionado con `bounda_water`: el
      clamp de `groundflow_phreatic_coord()` limita la presión estática a este
      umbral (default 0 → solo succión). Documentado en manual-developer/bounda_water.md.
- [ ] Contacto: `contact_apply`, `contact_heat_generation`, `contact_penalty_*`,
      `contact_plasti_friction`, `contact_target_*`.
- [ ] `control_reset_dof*`, `control_reset_value_*`, `change_dataitem_apply`,
      `change_dataitem_geometry`, `control_data_*`, `control_distribute*`.
      **NOTA 2026-08-12**: `control_data_put` y `control_distribute` están en
      database.cc (verificar estado); `control_reset_*` y `change_dataitem_*`
      no están.
- [ ] `control_mesh_generate_interface_geometry`, `slide_axisymmetric`,
      `materi_displacement_relative`, `group_interface_materi_plasti_tension_direct`.

### Fichero de control
- Inventario completo: `ProjectDocs/inventario-features-faltantes-2024.txt` (718 líneas).
- **Seguimiento exhaustivo**: `ProjectDocs/SEGUIMIENTO-CONVERGENCIA.md` — lista de TODAS
  las features de Professional vs GNU con estado, fecha/commit de cada feature
  implementada, features descartadas con razón, y registro de verificación.
- Marcar con `[x]` cada feature al implementarla y verificar su test.
- Al final de cada sesión: guardar el progreso en memoria (Engram) con el plan actualizado.

---

## 6c. Brecha de convergencia con Tochnog Professional — estudio de limitantes (2026-08-12)

### Métrica global

Comparando los keywords documentados en el manual de Professional
(`UserManual-professional.txt`, secciones `N.M keyword`) contra los
registrados en `database.cc` del GNU:

- **856** keywords en Professional
- **252** presentes en el GNU (**29 %**)
- **604** faltantes

La brecha no es uniforme: hay clusters de infraestructura que bloquean
familias enteras, y clusters puramente añaditivos (keywords nuevos sobre
infraestructura existente).

### Métrica de ejecución (la que importa para convergencia)

- **197 tests del GNU 2014** (`external-downloads/sfnet/extracted/test`,
  `test-2014.zip`) corren hoy con **cobertura 100 % de keywords**: los
  364 keywords que usan esos tests están todos reconocidos por el GNU
  (los 4 aparentemente "faltantes" — `echo`, `derivatives`,
  `end_initia`, `number_of_space_dimensions` — son secciones del bloque
  `initia` que el parser trata aparte, no keywords de `database.cc`).
- Por tanto el GNU es **retro-compatible con su propia suite** (2014).
- **La brecha de convergencia con Professional son los inputs que no
  podemos correr**: los 604 keywords sin infraestructura o con
  infraestructura parcial. La métrica operativa es "qué porcentaje de un
  input de Professional arranca y da resultados físicamente correctos",
  no el conteo de keywords.
- **Fuente de verdad para priorizar**: el changelog de Professional
  (`external-downloads/professional/changes-site-archive.txt`). Las
  features recientes/activas son las que un usuario de Professional 2024
  esperaría.

### Limitantes de infraestructura (bloquean familias enteras)

1. **Elementos de interfaz (`interface`)** — **AUSENTE en el GNU**.
   Professional tiene toda una familia (conversión automática bar2/bar3/
   tria3/tria6 → quad4/quad6/prism6/prism12 con `control_mesh_convert`, y
   ~15 keywords `group_interface_*`: rigidez kn/kt, Mohr-Coulomb directo,
   tension_direct, gap, memoria, conductividad, groundflow). El GNU no
   tiene ningún enum INTERFACE. **Bloquea**: `control_print_interface_stress*`,
   `control_reset_interface*`, `interface_gap_apply`, toda la familia
   `group_interface_*`, y los elementos de interfaz en sí.
   → Requiere: nuevo tipo de elemento + ensamblaje de rigidez de interfaz
   (esfuerzo alto, toca el núcleo del solver).
   **Evidencia de actualidad**: el changelog de Professional confirma que
   la familia de interfaz es activa y reciente (2022: "Re-introduced
   `group_interface_materi_plasti_tension_direct`", "Changed
   `group_interface_gap`"), no legacy — ver Carril A.

2. **Contacto plástico** — parcial en el GNU. Existen `contactspring` y
   `control_mesh_generate_contactspring*` (generación de elementos), y los
   enums `GROUP_CONTACTSPRING_*`, pero **sin keywords registrados** para la
   plasticidad de contacto (cohesión, fricción, dirección automática,
   memoria). **Bloquea**: `group_contact_spring_plasti_*`,
   `contact_plasti_friction`, `contact_apply`. → Requiere: ley constitutiva
   de contacto sobre los contactspring (medio-alto).

3. **`post_calcul` limitado** — el GNU tiene `post_calcul` genérico
   (operaciones sobre dofs), pero **no** `post_calcul_materi_stress_force`
   (integración de tensiones sobre cortes → fuerzas/momentos). **Bloquea**:
   `control_print_materi_stress_force`. → Requiere: rutina de integración
   de tensiones sobre secciones (medio).

4. **Beams/vigas con plasticidad** — faltan `group_beam_*`
   (direction_z, shear, force_moment_plasti) y `control_print_beam*`.
   El GNU tiene elementos beam pero la familia de propiedades plásticas y
   su post-proceso no están. → Medio.

### Clusters puramente añaditivos (sin infraestructura nueva)

- **`control_mesh_*`** (~50 faltantes): truss/beam generation, extrude,
  delete, remove, convert, merge. El GNU ya tiene el núcleo de
  `control_mesh_*`; la mayoría son variantes añaditivas.
- **`force_edge_*`** (~44 faltantes): factores de carga multilineales,
  proyectados, etc. Son variantes de `force_edge` ya existente.
- **`control_print_*`** (~60 faltantes tras descartar gid): dof, node,
  vtk_*, history_smooth — la mayoría son post-proceso sobre
  infraestructura existente (patrón print_tabular/gmsh/frd).
- **`group_materi_*`** — modelos de material nuevos (algunos ya hechos:
  masin, sanisand, wolfersdorff; faltan los de P4-E).

### Conclusión

- La convergencia total es **mucho mayor de 604 keywords**: la brecha real
  se mide en infraestructura, no en keywords.
- El limitante **más estructural es la familia de interfaz** (elementos +
  plasticidad + post-proceso), que falta por completo.
- Los siguientes limitantes por orden de impacto: contacto plástico,
  post_calcul de fuerzas, vigas plásticas.
- El resto (malla, fuerzas de borde, post-proceso) es aditivo y
  incremental con el patrón que ya dominamos.

### Estrategia de convergencia — carriles C → B → A (2026-08-12)

Prioridad acordada: **Carril C (aditivos, valor inmediato) → Carril B
(features recientes de Professional) → Carril A (infraestructura de
interfaz, proyecto con especificación propia)**.

#### Carril C — Clusters aditivos (valor inmediato, patrón ya dominado)

Keywords nuevos sobre infraestructura existente; cada uno con test en la
suite. Criterio de cierre: keyword implementado + documentado
(manual-user/developer) + test en verde.

- [x] `control_mesh_generate_truss`, `control_mesh_generate_beam` y variantes
  (`trussbeam`, `_loose`, `_macro`; el núcleo `control_mesh_*` ya existe).
  **VERIFICADO 2026-08-13** — ya implementados en generate.cc
  (`generate_beam_truss`); tests gen1 (truss) y genbeam1 (beam) añadidos
  al harness. Docs creados.
- [x] `control_print_dof` (post-proceso sobre nodos; dofs + coordenadas en
  ASCII plano). **HECHO 2026-08-12** (variantes `_id`/`_smooth_*`/`_line`
  pendientes).
- [x] `control_print_history_smooth` (suavizado de history; factible, sin
  dependencias nuevas — la única de las 3 originales del bloque P5 que no
  requiere infraestructura). **HECHO 2026-08-12**.
- [x] `force_element_edge_multi_linear_factor_x` (carga de borde
  multilineal sobre `force_edge` existente). **HECHO 2026-08-13**
  (multiplicador multilineal en x sobre force_element_edge).
- [x] `control_print_vtk_dof` (variante de print_vtk: filtrar los campos
  escritos). **HECHO 2026-08-13** (`_coord`/`_empty`/`_node_method`/
  `_other` pendientes).
- **Carril C COMPLETO (2026-08-13)**. Próximo: Carril B (features
  recientes del changelog de Professional).

#### Carril B — Features recientes de Professional (del changelog)

**ESTADO 2026-08-13**: 2 features completadas (`control_reset_dof`,
`control_change_dataitem_apply`); `bounda_time_until_force` ya existía;
los 4 restantes quedan **PENDIENTES PARA MÁS ADELANTE** (decisión del
usuario 2026-08-13 — se avanza al Carril A):

Del changelog (`changes-site-archive.txt`), un usuario de Professional
2024 esperaría estas. **Pendientes**:

- `materi_displacement_relative` (i.c.w. `materi_velocity_integrated`) —
  **PENDIENTE DE DISEÑO (2026-08-13)**: requiere guardar el desplazamiento
  de referencia en cada cambio de paso/reset y calcular la diferencia;
  toca múltiples rutinas de integración. Nivel alto.
- `strain_settlement_diagram*` (asientos dependientes de tensiones) —
  **PENDIENTE DE DISEÑO (2026-08-13)**: cálculo de asientos no trivial.
- [x] `control_reset_value_dof` (y `control_reset_dof`, `_value_constant`,
  `_value_dof_diagram`, `_value_method`). **HECHO 2026-08-13** (verificado:
  hisv0 reseteado a 0.55; diagrama sigyy→hisv0; métodos -use/-add/-multiply).
  Otras distribuciones `_value_*` (exponent/linear/log/power/sqrt/multi_linear)
  pendientes.
- `group_materi_plasti_mohr_coul_direct_normal` y
  `group_materi_plasti_tension_direct_normal` (+ `_automatic`) —
  **PENDIENTE DE DISEÑO (2026-08-13)**: el GNU no tiene la familia
  `mohr_coul_direct`/`tension_direct` en absoluto; implementarla requiere
  el modelo constitutivo completo con plano normal. Nivel alto.
- `bounda_time_until_*` (ya implementado: bounda_time_until_force)
- [x] `control_change_dataitem_apply`. **HECHO 2026-08-13** (con -no ignora
  change_dataitem: targets de hypo1 se cumplen; sin apply la geometry
  cambia y los targets fallan).
- `slide_axisymmetric` — **PENDIENTE DE DISEÑO (2026-08-13)**: se intentó
  implementar (escala 2πr de la fricción) pero la validación es
  problemática: sin fuerza normal la fricción de Coulomb es 0 y no hay
  efecto observable. Requiere un test axisimétrico con fuerza normal y
  targets propios. Revertido; queda pendiente.
- `control_mesh_generate_interface_geometry` (depende del Carril A)

#### Carril A — Infraestructura de interfaz (EN CURSO 2026-08-13)

Fase de diseño técnico completada. Ver `ProjectDocs/DESIGN-INTERFACES.md`
para el modelo físico, el análisis del codebase (cómo encaja en elem.cc
con el patrón de spring.cc), las 4 fases, y el test de validación
propuesto.

Estado de fases:
- Fase 1 (elemento elástico kn/kt + `group_interface`): **IMPLEMENTADA Y
  VALIDADA (2026-08-13)** — interface_element() en interface.cc; el signo
  del ensamblaje se corrigió (`-sign*stress*dir`); límites físicos
  correctos (kn→∞ soldado, kn→0 libre; ver DESIGN-INTERFACES.md).
- Estrategia de no-interpenetración decidida: **penalización implícita +
  control_timestep_iterations** (sin line-search/arc-length); rigidez
  actualizada por iteración en la Fase 3 (patrón contactspring).
- Fase 2 (conversión automática `control_mesh_convert`): **IMPLEMENTADA
  (2026-08-14, 3D 2026-08-16)** — bar2→quad4: crea 2 nodos duplicados
  desplazados en la normal, reescribe el elemento, reconecta vecinos del
  otro lado con `control_mesh_convert_element_group`. Validada (test
  iface_conv: nodos 7,8 en x=0.99 creados, bloques reconectados). 3D:
  tria3→prism6 y quad4→hex8 (test iface_conv3d). La llamada se movió a
  `step_start` para convertir antes del primer `element_loop`.
- Fase 3 (ley constitutiva): **IMPLEMENTADA (2026-08-13, +memory 2026-08-16)** —
  `group_interface_gap` (validado: con gap la interfaz abre y el bloque se
  separa velix=29.5, sin gap resiste velix=-3.24),
  `group_interface_materi_residual_stiffness`,
  `group_interface_materi_plasti_tension_direct`,
  `group_interface_materi_plasti_mohr_coul_direct` (MC/tension en
  validación de humo), `group_interface_materi_memory`
  (`-updated_linear`/`-total_linear`, commit `643865b`).
- Fase 4 (post-proceso): **IMPLEMENTADA (2026-08-16)** —
  `control_print_interface_stress` 2D (sign normal desde el strain
  acumulado, corte por línea; sigt leído del history
  `element_interface_force_tang`) y 3D (centroide + sign + sigt1 +
  sigt2, filtro `_3d_geometry`, orden `_3d_order -x/-y/-z`).
  Interfaz 3D (prism6/hex8, 2 tangentes, history `force_tang2`) y
  conversión 3D (tria3→prism6, quad4→hex8) implementadas (`d57faf4`).

La familia `group_interface_*` (11 keywords) está documentada en el
seguimiento. Sin tests de referencia en sfnet — validación con test
propio de 2 bloques.

**Ya presentes en el GNU** (marcar `[x]`): `control_mesh_switch`,
`groundflow_pressure_factor`, `group_groundflow_permeability_vertical_stress`,
`control_print_history_factor`, `control_print_data_versus_data_factor`,
`bounda_time_on_off`, `force_point`, `check_used`,
`bounda_factor_parabolic_x`, `group_materi_plasti_maximum_iterations`
(el changelog la lista como `materi_plasti_maximum_iterations*`).

#### Carril A — Infraestructura de interfaz (proyecto independiente)

Los interface elements son features **activas y recientes** de
Professional (changelog 2022: "Re-introduced
`group_interface_materi_plasti_tension_direct`", "Changed
`group_interface_gap`"), no legacy — la familia es una brecha real de
convergencia. Es el limitante más estructural y toca el núcleo del solver,
por eso merece especificación y diseño propios.

Fases propuestas:

1. Tipo de elemento de interfaz + ensamblaje de rigidez (kn/kt).
2. `control_mesh_convert` (bar2/bar3/tria3/tria6 → quad4/quad6/prism6/
   prism12) y `control_mesh_generate_interface_geometry`.
3. Familia `group_interface_*` (~15 keywords: rigidez, Mohr-Coulomb
   directo, tension_direct, gap, memoria, conductividad, groundflow).
4. Post-proceso: `control_print_interface_stress*`,
   `control_reset_interface*`, `interface_gap_apply`.

---

## 6b. Decisión de arquitectura: C vs C++ (2026-08-11)

### Contexto (evidencia)

- El proyecto se compila como C++ (`.cc` con g++), pero el C++ real usado es solo
  `ofstream`/`cout` (21 archivos de salida). Cero `std::vector`, `std::map`,
  RAII o templates en el código de proceso.
- Los kernels constitutivos (`hypo.c`, `masin.c`, `masin_visco.c`, `sanisand.c`)
  son C puro tras P4-F (reentrancia, ABI estable, port Fortran).
- El núcleo de cálculo usa arrays planos con indexación manual
  (`element_matrix[indx*nnol*npuknwn+indx]`) y `get_new_dbl/get_new_int`
  (= `malloc` sin RAII).
- `MDIM` es una constante de compilación (`#define MDIM 3`), lo que habilita
  templates de dimensión fija.

### Decisión

**C++ se usa donde aporta valor real; no hay una "evolución C++" global ni se
lleva todo a C puro.** El criterio es por capas, según el coste/beneficio:

| Capa | Lenguaje | Razón |
|------|----------|-------|
| Kernels constitutivos | **C puro** | Reentrancia, ABI estable, port Fortran (P4-F). Decisión cerrada. |
| Núcleo de cálculo (assembl/solve/element/db) | C plano (intocable) | Estable 40 años, usa LAPACK/BLAS; los templates solo aportarían claridad marginal a costa de riesgo de regresión. NO se reescribe. |
| **Capa nueva** (post-proceso, SQLite, magnitudes derivadas) | **C++ moderno** | Aquí C++ aporta: RAII para `SqliteDB`, `std::vector`/`std::map` para buffers y series, `template<int D> Tensor` para magnitudes con dimensiones fijas (2D/3D). No es ostream por ostream — es gestión de recursos y tipos seguros. |

### Razones de la decisión

1. **C puro no hace perder nada hoy**: el cuello de botella real es el
   assemblaggio/solve con LAPACK — un `std::vector` no acelera eso, y una
   reescritura C++ del núcleo es un riesgo enorme sin beneficio medido.
2. **C++ aporta donde la gestión de recursos y tipos importa**: RAII evita
   fugas de `sqlite3*`/buffers (el patrón `get_new_dbl`+`free` manual de 40
   años); los templates de dimensión fija (`Tensor<D>`) eliminan la indexación
   manual y los bounds checks en el código nuevo.
3. **Refactor C++ del ensamblaje SOLO si el profiling lo justifica**: si un
   día se mide un cuello de botella en los bucles de elementos, los templates
   `Tensor<MDIM>` pueden acelerarlos — pero de forma incremental por módulo,
   nunca como reescritura global.
4. **Decisiones cerradas se mantienen**: P4-F (kernels C puro) no se revierte.

### Consecuencia práctica

- Los kernels nuevos siguen en C puro.
- La capa P5 (SQLite, tabular, derivadas) se implementa en C++ moderno
  (RAII + contenedores + `Tensor<D>`).
- El núcleo no se toca salvo bug crítico.

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
