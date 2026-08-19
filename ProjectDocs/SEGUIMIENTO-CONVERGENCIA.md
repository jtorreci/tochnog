# Tochnog Professional vs GNU — Seguimiento de convergencia


Documento de control del proceso de convergencia entre Tochnog Professional
(manual `UserManual-professional.txt`) y el GNU (fork de sfnet 2014).

## Métrica global

- **1094** keywords documentados en Professional (tras limpiar 1 keyword
  espurio del manual)
- **304** presentes en el GNU (**28 %** de cobertura nominal)
- **790** faltantes
- Además: 197 tests del GNU 2014 corren con cobertura 100 % de keywords (ver plan, sección 6c)

## Features implementadas por nosotros — registro de verificación

Cada feature implementada en este proceso (posterior al fork sfnet 2014)
se verifica contra un ejemplo del manual de Professional, un test de la
suite sfnet, o un test propio. El registro completo:

| Feature | Commit | Fecha | Verificación |
|---------|--------|-------|--------------|
| `check_used` | `3793892` | 2026-08-04 | test propio |
| `check_element_shape`, `check_memory` | `762f2ca` | 2026-08-04 | test propio |
| `check_solver` | `e42e12d` | 2026-08-04 | test propio |
| `check_data`, `check_error`, `check_warning` | `5a2b3d2` | 2026-08-04 | test propio |
| `control_mesh_move/mirror/copy` | `3e94dac` | 2026-08-05 | test propio |
| `control_mesh_rotate`, elementos 3D | `d8d3a72` | 2026-08-05 | tests 3D |
| `control_mesh_keep_*`, `change_element_group`, `delete_element` | `afc1dad` | 2026-08-05 | test propio |
| `control_mesh_remove` | `1e5b4b9` | 2026-08-05 | test propio |
| `control_mesh_extrude` | `5dec3c2` | 2026-08-05 | test propio |
| `bounda_unknown` (P3) | `99d6cd6` | 2026-08-05 | test propio |
| `bounda_dof` (alias unknown) | `b5e02d5` | 2026-08-05 | test propio |
| `bounda_alternate` | `0ba382a` | 2026-08-05 | test propio |
| `bounda_normal` | `d20ca47` | 2026-08-05 | test propio |
| `bounda_water` | `e7861f6` | 2026-08-05 | test propio |
| `bounda_dof_radial`, `_cylindrical` | `dd53575` | 2026-08-05 | test propio |
| `bounda_time_units`, `geometry_method`, `factor_parabolic_x` | `ba20c45` | 2026-08-05 | test propio |
| `groundflow_pressure_factor` | `75d1ac1` | 2026-08-05 | test propio |
| `group_groundflow_permeability_vertical_stress` | `07f77ac` | 2026-08-05 | test propio |
| `group_materi_plasti_hypo_masin` | `974b284` | 2026-08-06 | tests hypomasin1-5 contra `reference-masin` (paths de tensión) |
| `..._masin_clay_visco` | `a0f55a6` | 2026-08-06 | test propio (visco Niemunis Dr/Iv) |
| `..._masin_clay_anisotropic` | `4225318` | 2026-08-06 | test propio (anisotropía + intergranular strain) |
| `group_materi_plasti_sanisand` | `b6034cc` | 2026-08-07 | contra UMAT de Dafalias & Manzari (P4-E1, hasta 5.4% de deviator residual) |
| `control_print_tabular` (P5-T1/T2) | `1c677bb`/`7904b4c` | 2026-08-11 | contra valores de referencia del test |
| `control_print_history_smooth` | `8f8ad67` | 2026-08-12 | media móvil verificada numéricamente (test smooth1) |
| `control_print_dof` | `ca2a9e4` | 2026-08-12 | contra SQLite de tabular (test dof1) |
| `control_print_gmsh` (familia) | `90798f9` | 2026-08-12 | validado con gmsh 4.12 (test gmsh) |
| `control_print_frd` (familia) | `19a8d08` | 2026-08-12 | formato byte-a-byte contra CalculiX (test frd) |
| `control_print_vtk_dof` | `b3e58a4` | 2026-08-13 | campos filtrados verificados (test vtk_dof1) |
| `force_element_edge_multi_linear_factor_x` | `86f6d03` | 2026-08-13 | efecto del factor 0/1 en target (test mlx1) |
| `control_mesh_generate_truss/beam` | `bd7b3ab` | 2026-08-13 | contra tests sfnet genera1/ho_othr1 (gen1/genbeam1) |
| `control_reset_dof` + `_value_constant` + `_value_dof` + `_value_dof_diagram` + `_value_method` | `b6eaee4` | 2026-08-13 | hisv0 reseteado a 0.55 (0.5986 sin reset); diagrama sigyy→hisv0; métodos -use/-add/-multiply |
| `control_change_dataitem_apply` | `bdcdeaf` | 2026-08-13 | con -no ignora change_dataitem (targets hypo1 se cumplen); sin apply la geometry cambia y los targets fallan |
| Carril A (diseño) | — | 2026-08-13 | diseño técnico en DESIGN-INTERFACES.md: modelo físico (strain = dif. de desplazamiento entre lados), análisis del codebase (patrón spring.cc en elem.cc), 4 fases, test de validación propuesto. Sin tests de referencia sfnet. |
| Carril A Fase 1 (`group_interface` + `_elasti_stiffness`) | `a82cbc7` | 2026-08-13 | interface_element() en interface.cc. **Validada** (test 2 bloques): signo corregido (-sign*stress*dir); kn=100→0.044, kn=1e6→-0.003≈soldado, kn=0.001→≈1.0 libre. Límites físicos correctos. Estrategia de no-interpenetración: penalización implícita + control_timestep_iterations (sin line-search/arc-length). |
| Carril A Fase 3 (gap, residual_stiffness, tension_direct, mohr_coul_direct — fix y validación) | `9c2f4c8` (código T1-T5; docs T8 + build_safe T7 en el commit del orquestador) | 2026-08-14 | VALIDADO (familia `iface_mc`, 13 tests / 18 runs en verde). Mohr-Coulomb acumulativo (history `element_interface_force_tang`, activación por presencia del record, phi_flow = dilatancia); gap: cerrada si strain > gap, default -1e20 (siempre cerrada), hueco físico = gap negativo; tension_direct abre en tracción sobre la fuerza normal TOTAL. a/a' y b/b' invariantes en nº de pasos. Corrige los bugs de `60bf78c` (2026-08-13). |
| Carril A Fase 2 (`control_mesh_convert` bar2→quad4) | `490545b` | 2026-08-14 | validado (test iface_conv): nodos 7,8 en x=0.99 creados, elemento reescrito quad4, bloques reconectados. |
| Carril A Fase 4 (`control_print_interface_stress` 2D) | `01f6c3e` | 2026-08-14 | validado (test iface_stress): interface_stress.0 generado con distancia+sign (strain acumulado * kn); sign crece con la compresion. sigt=0 y 3D pendientes. |

**Carril A Fase 4 — sigt implementado** | `8ef45c4` | 2026-08-14 | `interface_sigt` ya no es 0: se lee del history `element_interface_force_tang` (fuerza tangencial total acumulada, Fase 3), consistente con `sign`. Verificado con probe: sigt=23.28 == F_t acumulada del último paso. Documentado en ambos manuales. 3D pendiente. |

**Carril A — `group_interface_materi_memory`** | `643865b` | 2026-08-16 | modelo de memoria de la ley de interfaz: `-updated_linear` (default, normal/tangente de la config actual) o `-total_linear` (geometría de referencia tiempo 0, `NODE_START_REFINED`). Valores inválidos → `db_error`. Test `iface_mc_mem` (familia iface_mc, 7º run → 19 runs). |

**Carril A — validación numérica MC/tension completa** | `8022af3` | 2026-08-16 | cierre de RF-4: `iface_mc_dil` (phi_flow=30° → strain_normal acumulado −2.8065 vs +0.0802 sin dilatancia; apertura por deslizamiento > compresión), `iface_mc_dil_1step` (invarianza paso, −2.82675, dif ~0.7%) y `iface_mc_num` (clamp MC: c=50, phi=0, vely=1000 → `element_interface_force_tang`=50.0 exacto = max_fric=c). Familia iface_mc → 10 runs (22 totales). |

**Carril A — interfaz 3D + conversión 3D (items C+D)** | `d57faf4` | 2026-08-16 | `interface_element` soporta 3D (prism6/hex8): normal real por cross product orientada de lado2 hacia lado1 (compresión positiva), 2 tangentes ortonormales, 2 fuerzas tangenciales acumuladas (nueva history `element_interface_force_tang2`). MC 3D clampa la magnitud de la fuerza tangencial total preservando la razón; dilatancia usa |du_tang| en ambos planos. `interface_convert` ahora convierte tria3→prism6 y quad4→hex8; la llamada se movió de `step_close` a `step_start` (antes del primer `element_loop`). Tests: `iface_3d` (fricción alta sostiene carga tangencial 3D), `iface_3d_slip` (fricción nula desliza), `iface_conv3d` (tria3→prism6). |

**Carril A — print 3D `control_print_interface_stress` (item E)** | `b8804c6` | 2026-08-16 | 3D implementado: imprime centroide (x y z) + `sign` + `sigt1` + `sigt2` por elemento de interfaz. `_3d_geometry` filtra por geometría (patrón adjust.cc); `_3d_order -x/-y/-z` ordena las líneas por la coordenada del centroide (recolección+sort). Validado con `iface_3d_stress` (centroide 1.5 0.5 0.5, sign crece con compresión, sigt2≈0) e `iface_3d_order` (2 interfaces ordenadas por x). |

**Carril B — `control_mesh_generate_interface` (+ `_geometry` + `_method` + subdivisión cuadrática)** | `947ac83` + `87f157d` + `3d9ed64` + `ced696f` + `9a2b422` + `14a94d3` | 2026-08-17 | genera elementos de interfaz entre dos grupos de elementos que comparten una cara (nodos duplicados espacialmente). Sintaxis `index eg_i eg_a eg_b eg_j eg_c eg_d...`; el elemento generado se asigna a `eg_i`. `_geometry` filtra por geometría; `_method` (`-element_geometry`) selecciona por `element_geometry` y/o genera un registro `element_geometry`. **Subdivisión de caras cuadráticas** (`interface_face_subdivide`): quad9/bar3 (arista 3 nodos → 2 `quad4`), tet10 (tria6 → 4 `prism6`), hex27 (quad9 → 4 `hex8`) — acopla todos los nodos incl. medios/centro; clasificación por coordenadas. Validado con `iface_gen*` (lineal, geometría, método) e `iface_gen_quad9` (2 quad9 → 2 interfaces quad4, fricción MC sostiene). |

**Carril B — `group_materi_plasti_mohr_coul_direct`/`tension_direct` (+ `_normal`/`_normal_automatic`/`_visco`/`_wall`)** | `b6d4be6` + `7f51812` + `e553b6d` + `d3ba0cb` | 2026-08-17 | leyes plásticas de **cut-off directo de tensiones** en un plano con normal `n` (patrón del manual: sin deformaciones plásticas). `tension_direct sigy` capa la tracción normal; `mohr_coul_direct phi c phi_flow` capa el cortante a `max_fric = max(c - sig_n*tan(phi),0)`. `_normal` da la normal explícita; `_normal_automatic -yes` la toma del elemento; `_visco tm` relaja el cut-off (`factor=1-exp(-dt/tm)`); `_wall` usa valores alternativos si el elemento está pegado a una pared (`plasti_on_boundary`). `materi_direct_cutoff()` en stress.cc. Validado con `materi_direct`, `materi_direct_mc`, `materi_direct_auto`, `materi_direct_visco` y `materi_direct_wall`. |

**Carril B — `materi_displacement_relative`** | `96eabda` + `841c39e` | 2026-08-17 | opción del initia que añade un dof de **desplazamiento relativo** (`disr*`): acumula el desplazamiento desde un punto de referencia. La referencia se re-sincroniza (reset a 0) en dos eventos: cambio de timestep en `control_timestep` (dt nuevo ≠ persistido en `MATERI_DISPLACEMENT_RELATIVE_REF`, en `top.cc`) y reset de desplazamiento en `control_reset_dof` (cuando se resetea `disx`, en `data.cc`). Integración del dof junto a `dis_indx` en `dof.cc`. Requiere `materi_displacement` + `materi_velocity` + `materi_velocity_integrated`. Validado con `mat_rel` (dt cambia: disy=2.0, disry=1.0) y `mat_rel_reset` (reset de disx: disy=1.0, disry=0.1). |

**Carril B — `slide_axisymmetric` + variantes `control_reset_value_*` espaciales** | `3024d24` + `650ed26` | 2026-08-18 | `slide_axisymmetric -yes` escala la fricción de `slide_geometry` por `2*pi*r` (r = coordenada radial del nodo, `slide.cc`; validado con `slide_axi`). Variantes espaciales de `control_reset_value` (distribuciones en x/y/z): `_linear` (`ax x + ay y + az z`), `_power` (`ax x^bx + ...`), `_square_root`, `_exponent`, `_logarithmic`, `_logarithmic_second`, `_multi_linear` (tabla vs coordenada vertical) — implementadas en `data.cc` (bloque de reset); validada `_linear` con `reset_value_linear` (disx → x). |

**Carril B — `mesh_activate_gravity_time*` + `strain_settlement_*`** | `f6b8935` + `b58fd8f` + `4fc3419` + `6f44549` | 2026-08-18 | **Activación gradual de la gravedad** (`mesh_activate_gravity_time`, `_element`, `_element_group`, `_geometry`, `_time_initial`, `_method` 1 y 2, `_stiffness_factor`, `control_mesh_activate_gravity_apply`): `mesh_activate_gravity_factor()` en mesh.cc, aplicada en materi.cc (el vector de gravedad se multiplica por el factor de activación; el método 2 escala la matriz del elemento por el factor de rigidez). **Creep de asentamiento** (`strain_settlement_parameters`, `_element_group`, `strain_settlement_diagram`+`_dof`+`_number`): `strain_settlement_creep()` en materi.cc añade el creep vertical a `inc_ept` en `set_deften_etc`. Ley con saturación (decisión 2026-08-18, el OCR del manual es ambiguo): `eps_zz = Ar*(t/t_ref)^n/(t_plus+(t/t_ref)^n)`. Validado con `mesh_act_grav`, `mesh_act_grav2` (método 2), `strain_settle` (creep comprime la columna) y `strain_settle_diag` (Ar duplicado vía diagrama). |

**Carril B — familia `contact_*` (apply, plasti_friction, targets)** | `(commit en curso)` | 2026-08-19 | el GNU ya tenía el algoritmo de contacto (`parallel_contact` en contact.cc: geometría target, penalties, stick, relaxation, heatgeneration) y `contact_friction` simple. Se completan: `contact_apply` (gate `-yes`/`-no` por timestep), `contact_plasti_friction` (Mohr-Coulomb `max(c + Fn*tan(phi),0)` en el slip, reemplaza el `mu*Fn`), `contact_target_element_group` (filtro de targets por grupo), `contact_target_geometry`/`_switch` (alias de `contact_geometry`/`_switch`). Validado con `contact` (la familia se parsea y el modelo corre estable; el algoritmo de contacto es experimental y la detección de penetración depende del caso). |

Nota: las features marcadas solo "GNU" en el detalle por familia (sin
fila en esta tabla) ya existían en el fork sfnet 2014 y no fueron
implementadas por nosotros.

## Convenciones del seguimiento

- `ESTADO`: `[x]` = implementada (commit abajo); `GNU` = presente en el codebase GNU
  (puede faltar doc); `PENDIENTE` = falta; `DESCARTADO` = no se implementará (razón abajo).
- `Fecha/Commit`: fecha ISO del commit que implementó la feature (no hay números de versión).
- `Verificación`: test de la suite, ejemplo del manual, o test propio.

> **Regla de trabajo**: cada feature nueva se investiga primero — ¿hay ejemplo de uso en el
> manual de Professional, en su repo/Google Drive, o en la suite sfnet? Si lo hay, se verifica
> contra ese ejemplo; si no, se monta un test propio. La verificación queda anotada aquí y en
> el manual-developer de la feature.

## Descarte de features

Las features que decidimos NO implementar se anotan aquí con su razón, y además se reflejan
en el apartado 'Diferencias con la versión Professional' del manual correspondiente.

| Feature | Razón del descarte |
|---------|--------------------|
| `control_print_gid_*` (familia ~20) | Único formato propietario contemplado; `print_gid_6` del GNU es de 1998 (anterior a cambios recientes); GiD puede importar formatos no nativos (VTK, Gmsh, CSV). Ver plan sección 6c. |
| `control_print_interface_stress*` | Requiere elementos de interfaz/contacto, que el GNU no tiene (ni enums). Depende del Carril A. |
| `control_print_materi_stress_force` | Requiere `post_calcul_materi_stress_force` (integración de tensiones sobre cortes), no existente. |

## Equivalencias de nomenclatura Professional vs GNU

Algunas familias tienen nombres distintos entre Professional y el GNU
(fork sfnet 2014). En el detalle por familia, un keyword de Professional
marcado `PENDIENTE` puede estar cubierto en el GNU bajo otro nombre:

| Professional | GNU | Notas |
|--------------|-----|-------|
| `force_edge*` | `force_element_edge*` | Fuerzas de borde distribuidas. El GNU usa `force_element_edge`, `force_element_edge_factor`, etc. |
| `materi_plasti_maximum_iterations` | `group_materi_plasti_maximum_iterations` | Prefijo `group_` en el GNU. |
| `control_mesh_convert` (interfaces) | — | Solo conversión de interfaz; el GNU no tiene elementos de interfaz (Carril A). |

## Detalle por familia

**Leyenda**: `[x]` = implementada por nosotros (commit abajo); `GNU` = presente en el GNU
(posiblemente sin doc); `PENDIENTE` = falta.

<details><summary>Ver detalle completo</summary>




### area_element (5/13)

- [x] `area_element_group` — presente en el GNU
- [ ] `area_element_group_element` — PENDIENTE
- [ ] `area_element_group_interface` — PENDIENTE
- [ ] `area_element_group_method` — PENDIENTE
- [ ] `area_element_group_node` — PENDIENTE
- [x] `area_element_group_sequence` — presente en el GNU
- [x] `area_element_group_sequence_element` — presente en el GNU
- [ ] `area_element_group_sequence_element_group` — PENDIENTE
- [x] `area_element_group_sequence_geometry` — presente en el GNU
- [ ] `area_element_group_sequence_geometry_method` — PENDIENTE
- [ ] `area_element_group_sequence_interface` — PENDIENTE
- [x] `area_element_group_sequence_time` — presente en el GNU
- [ ] `area_element_group_time` — PENDIENTE

### area_node (3/3)

- [x] `area_node_dataitem` — presente en el GNU
- [x] `area_node_dataitem_double` — presente en el GNU
- [x] `area_node_dataitem_integer` — presente en el GNU

### axisymmetric (0/1)

- [ ] `axisymmetric` — PENDIENTE

### beam (1/1)

- [x] `beam_rotation` — presente en el GNU

### bounda_alternate (1/1)

- [x] `bounda_alternate` — implementada (commit `0ba382a`, 2026-08-05)

### bounda_baseline (0/2)

- [ ] `bounda_baseline_correction` — PENDIENTE
- [ ] `bounda_baseline_correction_parameters` — PENDIENTE

### bounda_constant (1/1)

- [x] `bounda_constant` — presente en el GNU

### bounda_dof (3/3)

- [x] `bounda_dof` — implementada (commit `b5e02d5`, 2026-08-05)
- [x] `bounda_dof_cylindrical` — implementada (commit `dd53575`, 2026-08-05)
- [x] `bounda_dof_radial` — implementada (commit `5058e46`, 2026-08-05)

### bounda_factor (2/2)

- [x] `bounda_factor` — presente en el GNU
- [x] `bounda_factor_parabolic_x` — implementada (commit `ba20c45`, 2026-08-05)

### bounda_force (1/1)

- [x] `bounda_force` — presente en el GNU

### bounda_found (1/1)

- [x] `bounda_found` — presente en el GNU

### bounda_geometry (1/1)

- [x] `bounda_geometry_method` — implementada (commit `ba20c45`, 2026-08-05)

### bounda_normal (1/1)

- [x] `bounda_normal` — implementada (commit `d20ca47`, 2026-08-05)

### bounda_print (0/3)

- [ ] `bounda_print_mesh_dof` — PENDIENTE
- [ ] `bounda_print_mesh_dof_geometry` — PENDIENTE
- [ ] `bounda_print_mesh_dof_values` — PENDIENTE

### bounda_sine (1/1)

- [x] `bounda_sine` — presente en el GNU

### bounda_time (4/11)

- [x] `bounda_time` — presente en el GNU
- [ ] `bounda_time_factor` — PENDIENTE
- [x] `bounda_time_increment` — presente en el GNU
- [ ] `bounda_time_o` — PENDIENTE
- [ ] `bounda_time_smc` — PENDIENTE
- [ ] `bounda_time_smc_o` — PENDIENTE
- [ ] `bounda_time_smc_units` — PENDIENTE
- [x] `bounda_time_units` — implementada (commit `ba20c45`, 2026-08-05)
- [ ] `bounda_time_until_data` — PENDIENTE
- [ ] `bounda_time_until_value_minimum` — PENDIENTE
- [x] `bounda_time_user` — presente en el GNU

### bounda_water (1/1)

- [x] `bounda_water` — implementada (commit `e7861f6`, 2026-08-05)

### change (4/6)

- [x] `change_dataitem` — presente en el GNU
- [ ] `change_dataitem_geometry` — PENDIENTE
- [x] `change_dataitem_time` — presente en el GNU
- [x] `change_dataitem_time_discrete` — presente en el GNU
- [ ] `change_dataitem_time_method` — PENDIENTE
- [x] `change_dataitem_time_user` — presente en el GNU

### check_data (1/1)

- [x] `check_data` — implementada (commit `5a2b3d2`, 2026-08-04)

### check_element (2/2)

- [x] `check_element_node` — implementada (commit `3793892`, 2026-08-04)
- [x] `check_element_shape` — implementada (commit `762f2ca`, 2026-08-04)

### check_error (1/1)

- [x] `check_error` — implementada (commit `5a2b3d2`, 2026-08-04)

### check_memory (3/3)

- [x] `check_memory` — implementada (commit `762f2ca`, 2026-08-04)
- [x] `check_memory_usage` — presente en el GNU
- [x] `check_memory_usage_result` — presente en el GNU

### check_nan (1/1)

- [x] `check_nan` — implementada (commit `3793892`, 2026-08-04)

### check_solver (1/1)

- [x] `check_solver` — implementada (commit `e42e12d`, 2026-08-04)

### check_target (1/1)

- [x] `check_target` — implementada (commit `3793892`, 2026-08-04)

### check_used (1/1)

- [x] `check_used` — implementada (commit `3793892`, 2026-08-04)

### check_warning (1/1)

- [x] `check_warning` — implementada (commit `5a2b3d2`, 2026-08-04)

### condif_convection (0/7)

- [ ] `condif_convection_edge_normal` — PENDIENTE
- [ ] `condif_convection_edge_normal_element` — PENDIENTE
- [ ] `condif_convection_edge_normal_element_group` — PENDIENTE
- [ ] `condif_convection_edge_normal_element_node` — PENDIENTE
- [ ] `condif_convection_edge_normal_element_side` — PENDIENTE
- [ ] `condif_convection_edge_normal_geometry` — PENDIENTE
- [ ] `condif_convection_edge_normal_node` — PENDIENTE

### condif_heat (0/20)

- [ ] `condif_heat_edge_normal` — PENDIENTE
- [ ] `condif_heat_edge_normal_element` — PENDIENTE
- [ ] `condif_heat_edge_normal_element_group` — PENDIENTE
- [ ] `condif_heat_edge_normal_element_node` — PENDIENTE
- [ ] `condif_heat_edge_normal_element_node_factor` — PENDIENTE
- [ ] `condif_heat_edge_normal_element_side` — PENDIENTE
- [ ] `condif_heat_edge_normal_factor` — PENDIENTE
- [ ] `condif_heat_edge_normal_geometry` — PENDIENTE
- [ ] `condif_heat_edge_normal_node` — PENDIENTE
- [ ] `condif_heat_edge_normal_sine` — PENDIENTE
- [ ] `condif_heat_edge_normal_time` — PENDIENTE
- [ ] `condif_heat_volume` — PENDIENTE
- [ ] `condif_heat_volume_element` — PENDIENTE
- [ ] `condif_heat_volume_element_group` — PENDIENTE
- [ ] `condif_heat_volume_factor` — PENDIENTE
- [ ] `condif_heat_volume_geometry` — PENDIENTE
- [ ] `condif_heat_volume_sine` — PENDIENTE
- [ ] `condif_heat_volume_time` — PENDIENTE
- [ ] `condif_heat_volume_user` — PENDIENTE
- [ ] `condif_heat_volume_user_parameters` — PENDIENTE

### condif_radiation (0/7)

- [ ] `condif_radiation_edge_normal` — PENDIENTE
- [ ] `condif_radiation_edge_normal_element` — PENDIENTE
- [ ] `condif_radiation_edge_normal_element_group` — PENDIENTE
- [ ] `condif_radiation_edge_normal_element_node` — PENDIENTE
- [ ] `condif_radiation_edge_normal_element_side` — PENDIENTE
- [ ] `condif_radiation_edge_normal_geometry` — PENDIENTE
- [ ] `condif_radiation_edge_normal_node` — PENDIENTE

### condif_temperature (1/1)

- [x] `condif_temperature` — presente en el GNU

### contact_apply (1/1)

- [x] `contact_apply` — implementada (gate `-yes`/`-no` en `parallel_contact`; validada con `contact`)

### contact_heat (0/1)

- [ ] `contact_heat_generation` — PENDIENTE

### contact_penalty (3/3)

- [x] `contact_penalty_pressure` — presente en el GNU
- [x] `contact_penalty_temperature` — presente en el GNU
- [x] `contact_penalty_velocity` — presente en el GNU

### contact_plasti (1/1)

- [x] `contact_plasti_friction` — implementada (Mohr-Coulomb `max(c + Fn*tan(phi), 0)` en el slip)

### contact_target (3/3)

- [x] `contact_target_element_group` — implementada (filtro de targets por grupo)
- [x] `contact_target_geometry` — implementada (alias de `contact_geometry`)
- [x] `contact_target_geometry_switch` — implementada (alias de `contact_geometry_switch`)

### control_bounda (0/2)

- [ ] `control_bounda_relax` — PENDIENTE
- [ ] `control_bounda_relax_geometry` — PENDIENTE

### control_change (0/1)

- [x] `control_change_dataitem_apply` — implementada (commit `bdcdeaf`, 2026-08-13)

### control_check (0/1)

- [ ] `control_check_data` — PENDIENTE

### control_contact (0/1)

- [ ] `control_contact_apply` — PENDIENTE

### control_convection (0/1)

- [ ] `control_convection_apply` — PENDIENTE

### control_data (4/12)

- [ ] `control_data_activate` — PENDIENTE
- [ ] `control_data_arithmetic` — PENDIENTE
- [ ] `control_data_arithmetic_double` — PENDIENTE
- [ ] `control_data_copy` — PENDIENTE
- [ ] `control_data_copy_factor` — PENDIENTE
- [ ] `control_data_copy_index` — PENDIENTE
- [ ] `control_data_copy_index_factor` — PENDIENTE
- [x] `control_data_delete` — presente en el GNU
- [x] `control_data_put` — presente en el GNU
- [x] `control_data_put_double` — presente en el GNU
- [x] `control_data_put_integer` — presente en el GNU
- [ ] `control_data_save` — PENDIENTE

### control_dependency (0/1)

- [ ] `control_dependency_apply` — PENDIENTE

### control_distribute (1/6)

- [x] `control_distribute` — presente en el GNU
- [ ] `control_distribute_correlation_distance` — PENDIENTE
- [ ] `control_distribute_correlation_length` — PENDIENTE
- [ ] `control_distribute_minimum_maximum` — PENDIENTE
- [ ] `control_distribute_parameters` — PENDIENTE
- [ ] `control_distribute_seed` — PENDIENTE

### control_element (0/2)

- [ ] `control_element_group` — PENDIENTE
- [ ] `control_element_group_apply` — PENDIENTE

### control_geometry (0/1)

- [ ] `control_geometry_moving` — PENDIENTE

### control_ground (0/1)

- [ ] `control_ground` — PENDIENTE

### control_inertia (0/1)

- [ ] `control_inertia_apply` — PENDIENTE

### control_input (0/1)

- [ ] `control_input` — PENDIENTE

### control_materi (2/15)

- [ ] `control_materi_damage_apply` — PENDIENTE
- [ ] `control_materi_dynamic` — PENDIENTE
- [ ] `control_materi_elasti_k0` — PENDIENTE
- [ ] `control_materi_failure_apply` — PENDIENTE
- [ ] `control_materi_plasti_hardsoil_gammap_initial` — PENDIENTE
- [x] `control_materi_plasti_hypo_masin_clay_ocr_apply` — presente en el GNU
- [x] `control_materi_plasti_hypo_masin_ocr_apply` — presente en el GNU
- [ ] `control_materi_plasti_hypo_niemunis_visco_ocr_apply` — PENDIENTE
- [ ] `control_materi_plasti_hypo_pressure_dependent_void_ratio` — PENDIENTE
- [ ] `control_materi_plasti_hypo_substepping` — PENDIENTE
- [ ] `control_materi_plasti_tension_apply` — PENDIENTE
- [ ] `control_materi_plasti_visco_apply` — PENDIENTE
- [ ] `control_materi_undrained_apply` — PENDIENTE
- [ ] `control_materi_updated_apply` — PENDIENTE
- [ ] `control_materi_viscosity_apply` — PENDIENTE

### control_mesh (35/91)

- [x] `control_mesh_activate_gravity_apply` — implementada (registrado; selección de records a aplicar)
- [x] `control_mesh_adjust_geometry` — presente en el GNU
- [x] `control_mesh_change_element_group` — implementada (commit `afc1dad`, 2026-08-05)
- [x] `control_mesh_convert` — implementada (commit `490545b`, 2026-08-14; Carril A Fase 2, bar2→quad4)
- [x] `control_mesh_convert_element_group` — implementada (commit `490545b`, 2026-08-14)
- [ ] `control_mesh_convert_quad9_quad6` — PENDIENTE
- [ ] `control_mesh_convert_tria6_tria3` — PENDIENTE
- [x] `control_mesh_copy` — implementada (commit `3e94dac`, 2026-08-05)
- [ ] `control_mesh_cut_force` — PENDIENTE
- [ ] `control_mesh_cut_geometry` — PENDIENTE
- [x] `control_mesh_delete_element` — implementada (commit `afc1dad`, 2026-08-05)
- [x] `control_mesh_delete_geometry` — presente en el GNU
- [ ] `control_mesh_delete_geometry_direct` — PENDIENTE
- [x] `control_mesh_delete_geometry_element` — presente en el GNU
- [ ] `control_mesh_delete_geometry_element_group` — PENDIENTE
- [x] `control_mesh_delete_geometry_factor` — presente en el GNU
- [ ] `control_mesh_delete_geometry_method` — PENDIENTE
- [ ] `control_mesh_delete_geometry_move_node` — PENDIENTE
- [ ] `control_mesh_delete_geometry_projection_type` — PENDIENTE
- [ ] `control_mesh_delete_geometry_stop` — PENDIENTE
- [ ] `control_mesh_delete_geometry_stop_geometry` — PENDIENTE
- [x] `control_mesh_delete_small` — presente en el GNU
- [ ] `control_mesh_duplicate_element_group` — PENDIENTE
- [ ] `control_mesh_element_group_apply` — PENDIENTE
- [x] `control_mesh_extrude` — implementada (commit `5dec3c2`, 2026-08-05)
- [ ] `control_mesh_extrude_contact_spring_element_group` — PENDIENTE
- [ ] `control_mesh_extrude_contact_spring_element_group_new` — PENDIENTE
- [ ] `control_mesh_extrude_direction` — PENDIENTE
- [ ] `control_mesh_extrude_element` — PENDIENTE
- [ ] `control_mesh_extrude_element_group_new` — PENDIENTE
- [x] `control_mesh_extrude_n` — presente en el GNU
- [x] `control_mesh_generate_beam` — implementada (commit `bd7b3ab`, 2026-08-13)
- [ ] `control_mesh_generate_contact_spring` — PENDIENTE
- [ ] `control_mesh_generate_contact_spring_element` — PENDIENTE
- [ ] `control_mesh_generate_contact_spring_element_group` — PENDIENTE
- [x] `control_mesh_generate_interface` — implementada (genera interfaces entre grupos con cara compartida; validada con `iface_gen`, `iface_gen_geom`, `iface_gen_geom_off`)
- [x] `control_mesh_generate_interface_geometry` — implementada (filtro por geometría del par de elementos)
- [x] `control_mesh_generate_interface_method` — implementada (`-element_geometry` para selección y/o generación; validada con `iface_gen_method` e `iface_gen_method_gen`)
- [x] `control_mesh_generate_spring1` — presente en el GNU
- [x] `control_mesh_generate_spring2` — presente en el GNU
- [x] `control_mesh_generate_truss` — implementada (commit `bd7b3ab`, 2026-08-13)
- [ ] `control_mesh_generate_truss_beam` — PENDIENTE
- [x] `control_mesh_generate_truss_beam_loose` — implementada (commit `bd7b3ab`, 2026-08-13)
- [x] `control_mesh_generate_truss_beam_macro` — implementada (commit `bd7b3ab`, 2026-08-13)
- [ ] `control_mesh_generate_truss_beam_separate` — PENDIENTE
- [ ] `control_mesh_gid_batch` — PENDIENTE
- [ ] `control_mesh_interface_triangle` — PENDIENTE
- [x] `control_mesh_keep_element` — implementada (commit `afc1dad`, 2026-08-05)
- [x] `control_mesh_keep_element_group` — implementada (commit `afc1dad`, 2026-08-05)
- [x] `control_mesh_keep_geometry` — presente en el GNU
- [x] `control_mesh_keep_node` — implementada (commit `415d26f`, 2026-08-05)
- [x] `control_mesh_macro` — presente en el GNU
- [ ] `control_mesh_macro_concentrate` — PENDIENTE
- [x] `control_mesh_macro_element` — presente en el GNU
- [x] `control_mesh_macro_parameters` — presente en el GNU
- [ ] `control_mesh_map` — PENDIENTE
- [x] `control_mesh_merge` — presente en el GNU
- [ ] `control_mesh_merge_eps_coord` — PENDIENTE
- [ ] `control_mesh_merge_geometry` — PENDIENTE
- [ ] `control_mesh_merge_geometry_not` — PENDIENTE
- [x] `control_mesh_merge_macro_generate` — presente en el GNU
- [x] `control_mesh_mirror` — implementada (commit `3e94dac`, 2026-08-05)
- [x] `control_mesh_move` — implementada (commit `3e94dac`, 2026-08-05)
- [x] `control_mesh_multiply` — presente en el GNU
- [ ] `control_mesh_re` — PENDIENTE
- [x] `control_mesh_remove` — implementada (commit `1e5b4b9`, 2026-08-05)
- [ ] `control_mesh_remove_geometry` — PENDIENTE
- [ ] `control_mesh_remove_keep_geometry` — PENDIENTE
- [ ] `control_mesh_remove_really` — PENDIENTE
- [ ] `control_mesh_remove_really_activate_all` — PENDIENTE
- [ ] `control_mesh_remove_really_activate_factor` — PENDIENTE
- [x] `control_mesh_renumber` — presente en el GNU
- [ ] `control_mesh_renumber_element_geometry_o` — PENDIENTE
- [ ] `control_mesh_renumber_element_group_o` — PENDIENTE
- [x] `control_mesh_rotate` — implementada (commit `d8d3a72`, 2026-08-05)
- [x] `control_mesh_rotate_angle` — implementada (commit `415d26f`, 2026-08-05)
- [x] `control_mesh_split` — presente en el GNU
- [ ] `control_mesh_split_element_from` — PENDIENTE
- [ ] `control_mesh_split_element_to` — PENDIENTE
- [x] `control_mesh_split_only` — presente en el GNU
- [x] `control_mesh_switch` — presente en el GNU
- [ ] `control_mesh_truss_distribute_mpc` — PENDIENTE
- [ ] `control_mesh_truss_distribute_mpc_air` — PENDIENTE
- [ ] `control_mesh_truss_distribute_mpc_dof` — PENDIENTE
- [ ] `control_mesh_truss_distribute_mpc_element_group_isoparametric` — PENDIENTE
- [ ] `control_mesh_truss_distribute_mpc_element_group_truss` — PENDIENTE
- [ ] `control_mesh_truss_distribute_mpc_exact` — PENDIENTE
- [ ] `control_mesh_truss_distribute_mpc_exact_minimal_length` — PENDIENTE
- [ ] `control_mesh_truss_distribute_mpc_exact_minimal_length_connect` — PENDIENTE
- [ ] `control_mesh_truss_distribute_mpc_geometry_isoparametric` — PENDIENTE
- [ ] `control_mesh_truss_distribute_mpc_geometry_truss` — PENDIENTE

### control_mpc (0/2)

- [ ] `control_mpc_apply` — PENDIENTE
- [ ] `control_mpc_element_group` — PENDIENTE

### control_plasti (0/1)

- [ ] `control_plasti_apply` — PENDIENTE

### control_post (0/3)

- [ ] `control_post` — PENDIENTE
- [ ] `control_post_apply` — PENDIENTE
- [ ] `control_post_element_force` — PENDIENTE

### control_print (21/85)

- [x] `control_print` — presente en el GNU
- [ ] `control_print_` — PENDIENTE
- [ ] `control_print_beam_force_moment` — PENDIENTE
- [ ] `control_print_beam_force_moment_coordinates` — PENDIENTE
- [ ] `control_print_beam_force_moment_switch` — PENDIENTE
- [x] `control_print_data_versus_data` — presente en el GNU
- [x] `control_print_data_versus_data_factor` — implementada (commit `3793892`, 2026-08-04)
- [x] `control_print_database` — presente en el GNU
- [ ] `control_print_database_method` — PENDIENTE
- [x] `control_print_dof` — implementada (commit `ca2a9e4`, 2026-08-12)
- [ ] `control_print_dof_id` — PENDIENTE
- [ ] `control_print_dof_line` — PENDIENTE
- [ ] `control_print_dof_line_coordinates` — PENDIENTE
- [ ] `control_print_dof_line_element_group` — PENDIENTE
- [ ] `control_print_dof_line_eps_iso` — PENDIENTE
- [ ] `control_print_dof_line_method` — PENDIENTE
- [ ] `control_print_dof_line_move` — PENDIENTE
- [ ] `control_print_dof_line_n` — PENDIENTE
- [ ] `control_print_dof_line_time` — PENDIENTE
- [ ] `control_print_dof_point` — PENDIENTE
- [ ] `control_print_dof_point_coordinates` — PENDIENTE
- [ ] `control_print_dof_point_time` — PENDIENTE
- [ ] `control_print_dof_rhside` — PENDIENTE
- [ ] `control_print_dof_smooth_dof` — PENDIENTE
- [ ] `control_print_dof_smooth_n` — PENDIENTE
- [x] `control_print_element` — presente en el GNU
- [ ] `control_print_element_method` — PENDIENTE
- [x] `control_print_frd` — implementada (commit `19a8d08`, 2026-08-12)
- [x] `control_print_frd_freecad` — implementada (commit `19a8d08`, 2026-08-12)
- [x] `control_print_frd_prepomax` — implementada (commit `19a8d08`, 2026-08-12)
- [ ] `control_print_frequency_timeinterval` — PENDIENTE
- [ ] `control_print_frequency_timestep` — PENDIENTE
- [x] `control_print_gid` — presente en el GNU
- [ ] `control_print_gid_batch` — DESCARTADO
- [ ] `control_print_gid_beam_vectors` — DESCARTADO
- [ ] `control_print_gid_beam_vectors_normal` — DESCARTADO
- [ ] `control_print_gid_contact_spring2` — DESCARTADO
- [ ] `control_print_gid_coord` — DESCARTADO
- [ ] `control_print_gid_dof` — DESCARTADO
- [ ] `control_print_gid_dof_calcul` — DESCARTADO
- [ ] `control_print_gid_element_group` — DESCARTADO
- [ ] `control_print_gid_element_mpc` — DESCARTADO
- [x] `control_print_gid_empty` — presente en el GNU
- [ ] `control_print_gid_group` — DESCARTADO
- [ ] `control_print_gid_mesh_activate_gravity` — DESCARTADO
- [ ] `control_print_gid_method` — DESCARTADO
- [ ] `control_print_gid_node_method` — DESCARTADO
- [ ] `control_print_gid_other` — DESCARTADO
- [ ] `control_print_gid_safety_slip_critical` — DESCARTADO
- [ ] `control_print_gid_save_di` — PENDIENTE
- [ ] `control_print_gid_smooth_dof` — DESCARTADO
- [ ] `control_print_gid_smooth_n` — DESCARTADO
- [ ] `control_print_gid_spring2` — DESCARTADO
- [ ] `control_print_gid_truss_vector` — DESCARTADO
- [ ] `control_print_gid_truss_vector_normal` — DESCARTADO
- [x] `control_print_gmsh` — implementada (commit `90798f9`, 2026-08-12)
- [x] `control_print_gmsh_dummy` — implementada (commit `90798f9`, 2026-08-12)
- [x] `control_print_gmsh_element_data` — implementada (commit `90798f9`, 2026-08-12)
- [x] `control_print_gmsh_node_method` — implementada (commit `90798f9`, 2026-08-12)
- [x] `control_print_history` — presente en el GNU
- [x] `control_print_history_factor` — implementada (commit `3793892`, 2026-08-04)
- [ ] `control_print_history_relative_time` — PENDIENTE
- [x] `control_print_history_smooth` — implementada (commit `8f8ad67`, 2026-08-12)
- [x] `control_print_interface_stress` — implementada (commit `01f6c3e`, 2026-08-14; Carril A Fase 4, 2D)
- [x] `control_print_interface_stress_2d_coordinates` — implementada (commit `01f6c3e`, 2026-08-14)
- [ ] `control_print_interface_stress_3d_geometry` — PENDIENTE
- [ ] `control_print_interface_stress_3d_order` — PENDIENTE
- [ ] `control_print_materi_stress_force` — PENDIENTE
- [ ] `control_print_mesh_dof` — PENDIENTE
- [ ] `control_print_node` — PENDIENTE
- [ ] `control_print_node_angular` — PENDIENTE
- [ ] `control_print_node_angular_middle` — PENDIENTE
- [ ] `control_print_node_geometry` — PENDIENTE
- [ ] `control_print_node_sort` — PENDIENTE
- [ ] `control_print_node_zero` — PENDIENTE
- [ ] `control_print_number_iterations` — PENDIENTE
- [ ] `control_print_partialname` — PENDIENTE
- [x] `control_print_tecplot` — presente en el GNU
- [x] `control_print_vtk` — presente en el GNU
- [ ] `control_print_vtk_coord` — PENDIENTE
- [x] `control_print_vtk_dof` — implementada (commit `b3e58a4`, 2026-08-13)
- [ ] `control_print_vtk_dof_calcul` — PENDIENTE
- [ ] `control_print_vtk_empty` — PENDIENTE
- [ ] `control_print_vtk_node_method` — PENDIENTE
- [ ] `control_print_vtk_other` — PENDIENTE

### control_repeat (1/3)

- [x] `control_repeat` — presente en el GNU
- [ ] `control_repeat_save` — PENDIENTE
- [ ] `control_repeat_save_calculate` — PENDIENTE

### control_reset (5/18)

- [x] `control_reset_dof` — implementada (commit `b6eaee4`, 2026-08-13)
- [ ] `control_reset_element_dof` — PENDIENTE
- [ ] `control_reset_element_group` — PENDIENTE
- [ ] `control_reset_geometry` — PENDIENTE
- [ ] `control_reset_interface` — PENDIENTE
- [ ] `control_reset_interface_strain` — PENDIENTE
- [ ] `control_reset_node` — PENDIENTE
- [x] `control_reset_value_constant` — implementada (commit `b6eaee4`, 2026-08-13)
- [x] `control_reset_value_dof` — implementada (commit `b6eaee4`, 2026-08-13)
- [x] `control_reset_value_dof_diagram` — implementada (commit `b6eaee4`, 2026-08-13)
- [x] `control_reset_value_exponent` — implementada (distribución espacial exponencial; validada con `reset_value_linear`)
- [x] `control_reset_value_linear` — implementada (distribución espacial lineal `ax x + ay y + az z`; validada con `reset_value_linear`)
- [x] `control_reset_value_logarithmic_` — implementada (distribución espacial logarítmica)
- [x] `control_reset_value_logarithmic_second` — implementada (distribución espacial logarítmica de segundo tipo)
- [x] `control_reset_value_method` — implementada (commit `b6eaee4`, 2026-08-13)
- [x] `control_reset_value_multi_linear` — implementada (tabla multilineal vs coordenada vertical)
- [x] `control_reset_value_power` — implementada (distribución espacial potencial `ax x^bx + ...`)
- [x] `control_reset_value_square_root` — implementada (distribución espacial raíz cuadrada)

### control_restart (1/1)

- [x] `control_restart` — presente en el GNU

### control_safety (0/1)

- [ ] `control_safety_slip` — PENDIENTE

### control_slide (0/3)

- [ ] `control_slide_damping_apply` — PENDIENTE
- [ ] `control_slide_plasti_apply` — PENDIENTE
- [ ] `control_slide_sti` — PENDIENTE

### control_solver (0/7)

- [ ] `control_solver` — PENDIENTE
- [ ] `control_solver_bicg_error` — PENDIENTE
- [ ] `control_solver_bicg_restart` — PENDIENTE
- [ ] `control_solver_bicg_stop` — PENDIENTE
- [ ] `control_solver_matrix_save` — PENDIENTE
- [ ] `control_solver_pardiso_ordering` — PENDIENTE
- [ ] `control_solver_pardiso_out_of_core` — PENDIENTE

### control_support (0/2)

- [ ] `control_support_edge_normal_damping_apply` — PENDIENTE
- [ ] `control_support_edge_normal_sti` — PENDIENTE

### control_system (0/1)

- [ ] `control_system_call` — PENDIENTE

### control_timestep (5/10)

- [x] `control_timestep` — presente en el GNU
- [ ] `control_timestep_adjust_minimum_iterations` — PENDIENTE
- [x] `control_timestep_iterations` — presente en el GNU
- [x] `control_timestep_iterations_automatic` — presente en el GNU
- [ ] `control_timestep_iterations_automatic_minimum_maximum_wished` — PENDIENTE
- [x] `control_timestep_iterations_automatic_stop` — presente en el GNU
- [x] `control_timestep_multiplier` — presente en el GNU
- [ ] `control_timestep_until_data` — PENDIENTE
- [ ] `control_timestep_until_maximum` — PENDIENTE
- [ ] `control_timestep_until_minimum` — PENDIENTE

### control_truss (0/1)

- [ ] `control_truss_rope_apply` — PENDIENTE

### control_zip (0/1)

- [ ] `control_zip` — PENDIENTE

### convection (0/2)

- [ ] `convection_apply` — PENDIENTE
- [ ] `convection_stabilization` — PENDIENTE

### data (0/5)

- [ ] `data_activate` — PENDIENTE
- [ ] `data_activate_time` — PENDIENTE
- [ ] `data_delete` — PENDIENTE
- [ ] `data_delete_time` — PENDIENTE
- [ ] `data_ignore` — PENDIENTE

### dependency_apply (0/1)

- [ ] `dependency_apply` — PENDIENTE

### dependency_diagram (1/1)

- [x] `dependency_diagram` — presente en el GNU

### dependency_geometry (0/1)

- [ ] `dependency_geometry` — PENDIENTE

### dependency_item (1/1)

- [x] `dependency_item` — presente en el GNU

### dependency_method (0/1)

- [ ] `dependency_method` — PENDIENTE

### dependency_number (0/1)

- [ ] `dependency_number` — PENDIENTE

### dependency_type (0/1)

- [ ] `dependency_type` — PENDIENTE

### derivatives (0/1)

- [ ] `derivatives` — PENDIENTE

### dof (1/3)

- [ ] `dof_element_dof` — PENDIENTE
- [x] `dof_label` — presente en el GNU
- [ ] `dof_limit` — PENDIENTE

### dtime (1/1)

- [x] `dtime` — presente en el GNU

### echo (0/1)

- [ ] `echo` — PENDIENTE

### element (1/1)

- [x] `element` — presente en el GNU

### element_beam (1/3)

- [x] `element_beam_direction` — presente en el GNU
- [ ] `element_beam_direction_z` — PENDIENTE
- [ ] `element_beam_force_moment` — PENDIENTE

### element_boundary (0/1)

- [ ] `element_boundary` — PENDIENTE

### element_contact (0/3)

- [ ] `element_contact_spring_direction` — PENDIENTE
- [ ] `element_contact_spring_force` — PENDIENTE
- [ ] `element_contact_spring_strain` — PENDIENTE

### element_dof (1/3)

- [x] `element_dof` — presente en el GNU
- [ ] `element_dof_initial` — PENDIENTE
- [ ] `element_dof_initial_speci` — PENDIENTE

### element_empty (1/1)

- [x] `element_empty` — presente en el GNU

### element_geometry (0/2)

- [ ] `element_geometry` — PENDIENTE
- [ ] `element_geometry_present` — PENDIENTE

### element_group (1/2)

- [x] `element_group` — presente en el GNU
- [ ] `element_group_apply` — PENDIENTE

### element_interface (0/7)

- [ ] `element_interface_intpnt_direction` — PENDIENTE
- [ ] `element_interface_intpnt_gap_status` — PENDIENTE
- [ ] `element_interface_intpnt_materi_tension_status` — PENDIENTE
- [ ] `element_interface_intpnt_strain` — PENDIENTE
- [ ] `element_interface_intpnt_strain_average` — PENDIENTE
- [ ] `element_interface_intpnt_stress` — PENDIENTE
- [ ] `element_interface_intpnt_stress_average` — PENDIENTE

### element_intpnt (0/7)

- [ ] `element_intpnt_dof` — PENDIENTE
- [ ] `element_intpnt_h` — PENDIENTE
- [ ] `element_intpnt_iso_coord` — PENDIENTE
- [ ] `element_intpnt_materi_plasti_hardsoil_gammap_initial` — PENDIENTE
- [ ] `element_intpnt_materi_undrained_pressure` — PENDIENTE
- [ ] `element_intpnt_method` — PENDIENTE
- [ ] `element_intpnt_npoint` — PENDIENTE

### element_middle (1/1)

- [x] `element_middle` — presente en el GNU

### element_normal (0/1)

- [ ] `element_normal` — PENDIENTE

### element_print (0/1)

- [ ] `element_print_group_data_values` — PENDIENTE

### element_spring (1/2)

- [x] `element_spring_force` — presente en el GNU
- [ ] `element_spring_strain` — PENDIENTE

### element_truss (2/4)

- [x] `element_truss_direction` — presente en el GNU
- [x] `element_truss_force` — presente en el GNU
- [ ] `element_truss_strain` — PENDIENTE
- [ ] `element_truss_strain_temperature` — PENDIENTE

### element_volume (1/1)

- [x] `element_volume` — presente en el GNU

### end (0/2)

- [ ] `end_data` — PENDIENTE
- [ ] `end_initia` — PENDIENTE

### force_edge (0/44)

- [ ] `force_edge` — PENDIENTE
- [ ] `force_edge_diagram` — PENDIENTE
- [ ] `force_edge_element` — PENDIENTE
- [ ] `force_edge_element_group` — PENDIENTE
- [ ] `force_edge_element_node` — PENDIENTE
- [ ] `force_edge_element_side` — PENDIENTE
- [ ] `force_edge_factor` — PENDIENTE
- [ ] `force_edge_geometry` — PENDIENTE
- [ ] `force_edge_multi_linear_factor_x` — PENDIENTE
- [ ] `force_edge_node` — PENDIENTE
- [ ] `force_edge_node_factor` — PENDIENTE
- [ ] `force_edge_normal` — PENDIENTE
- [ ] `force_edge_normal_element` — PENDIENTE
- [ ] `force_edge_normal_element_group` — PENDIENTE
- [ ] `force_edge_normal_element_node` — PENDIENTE
- [ ] `force_edge_normal_element_side` — PENDIENTE
- [ ] `force_edge_normal_factor` — PENDIENTE
- [ ] `force_edge_normal_geometry` — PENDIENTE
- [ ] `force_edge_normal_node` — PENDIENTE
- [ ] `force_edge_normal_node_factor` — PENDIENTE
- [ ] `force_edge_normal_sine` — PENDIENTE
- [ ] `force_edge_normal_time` — PENDIENTE
- [ ] `force_edge_projected` — PENDIENTE
- [ ] `force_edge_projected_element` — PENDIENTE
- [ ] `force_edge_projected_element_group` — PENDIENTE
- [ ] `force_edge_projected_element_node` — PENDIENTE
- [ ] `force_edge_projected_element_side` — PENDIENTE
- [ ] `force_edge_projected_factor` — PENDIENTE
- [ ] `force_edge_projected_geometry` — PENDIENTE
- [ ] `force_edge_projected_node` — PENDIENTE
- [ ] `force_edge_projected_node_factor` — PENDIENTE
- [ ] `force_edge_projected_sine` — PENDIENTE
- [ ] `force_edge_projected_time` — PENDIENTE
- [ ] `force_edge_sine` — PENDIENTE
- [ ] `force_edge_time` — PENDIENTE
- [ ] `force_edge_water` — PENDIENTE
- [ ] `force_edge_water_element` — PENDIENTE
- [ ] `force_edge_water_element_group` — PENDIENTE
- [ ] `force_edge_water_element_node` — PENDIENTE
- [ ] `force_edge_water_element_side` — PENDIENTE
- [ ] `force_edge_water_factor` — PENDIENTE
- [ ] `force_edge_water_geometry` — PENDIENTE
- [ ] `force_edge_water_node` — PENDIENTE
- [ ] `force_edge_water_time` — PENDIENTE

### force_gravity (2/3)

- [x] `force_gravity` — presente en el GNU
- [ ] `force_gravity_geometry` — PENDIENTE
- [x] `force_gravity_time` — presente en el GNU

### force_point (1/1)

- [x] `force_point` — implementada (commit `3793892`, 2026-08-04)

### force_volume (0/7)

- [ ] `force_volume` — PENDIENTE
- [ ] `force_volume_element` — PENDIENTE
- [ ] `force_volume_element_group_0` — PENDIENTE
- [ ] `force_volume_factor` — PENDIENTE
- [ ] `force_volume_geometry` — PENDIENTE
- [ ] `force_volume_sine` — PENDIENTE
- [ ] `force_volume_time` — PENDIENTE

### geometry_bounda (3/3)

- [x] `geometry_bounda_sine_x` — presente en el GNU
- [x] `geometry_bounda_sine_y` — presente en el GNU
- [x] `geometry_bounda_sine_z` — presente en el GNU

### geometry_boundary (0/1)

- [ ] `geometry_boundary` — PENDIENTE

### geometry_brick (1/1)

- [x] `geometry_brick` — presente en el GNU

### geometry_circle (2/3)

- [x] `geometry_circle` — presente en el GNU
- [ ] `geometry_circle_part` — PENDIENTE
- [x] `geometry_circle_segment` — presente en el GNU

### geometry_cylinder (2/4)

- [x] `geometry_cylinder` — presente en el GNU
- [ ] `geometry_cylinder_part` — PENDIENTE
- [ ] `geometry_cylinder_part_start_vector` — PENDIENTE
- [x] `geometry_cylinder_segment` — presente en el GNU

### geometry_element (0/4)

- [ ] `geometry_element_geometry` — PENDIENTE
- [ ] `geometry_element_geometry_method` — PENDIENTE
- [ ] `geometry_element_group` — PENDIENTE
- [ ] `geometry_element_group_method` — PENDIENTE

### geometry_ellipse (1/1)

- [x] `geometry_ellipse` — presente en el GNU

### geometry_exclude (0/1)

- [ ] `geometry_exclude` — PENDIENTE

### geometry_factor (0/1)

- [ ] `geometry_factor` — PENDIENTE

### geometry_hexahedral (0/1)

- [ ] `geometry_hexahedral` — PENDIENTE

### geometry_line (1/2)

- [x] `geometry_line` — presente en el GNU
- [ ] `geometry_line_eps_iso` — PENDIENTE

### geometry_list (0/1)

- [ ] `geometry_list` — PENDIENTE

### geometry_method (0/1)

- [ ] `geometry_method` — PENDIENTE

### geometry_moving (0/6)

- [ ] `geometry_moving` — PENDIENTE
- [ ] `geometry_moving_n` — PENDIENTE
- [ ] `geometry_moving_operat` — PENDIENTE
- [ ] `geometry_moving_operat_parameter` — PENDIENTE
- [ ] `geometry_moving_operat_time` — PENDIENTE
- [ ] `geometry_moving_parameter` — PENDIENTE

### geometry_mpc (0/1)

- [ ] `geometry_mpc` — PENDIENTE

### geometry_node (0/1)

- [ ] `geometry_node_type` — PENDIENTE

### geometry_point (1/1)

- [x] `geometry_point` — presente en el GNU

### geometry_polynomial (1/1)

- [x] `geometry_polynomial` — presente en el GNU

### geometry_projection (0/1)

- [ ] `geometry_projection_type` — PENDIENTE

### geometry_quadrilateral (1/2)

- [x] `geometry_quadrilateral` — presente en el GNU
- [ ] `geometry_quadrilateral_eps_iso` — PENDIENTE

### geometry_set (1/1)

- [x] `geometry_set` — presente en el GNU

### geometry_sphere (2/2)

- [x] `geometry_sphere` — presente en el GNU
- [x] `geometry_sphere_segment` — presente en el GNU

### geometry_tetrahedral (0/1)

- [ ] `geometry_tetrahedral` — PENDIENTE

### geometry_triangle (1/2)

- [x] `geometry_triangle` — presente en el GNU
- [ ] `geometry_triangle_eps_iso` — PENDIENTE

### global (0/4)

- [ ] `global_element_dof_apply` — PENDIENTE
- [ ] `global_element_dof_from_node_dof` — PENDIENTE
- [ ] `global_node_dof_empty` — PENDIENTE
- [ ] `global_post_point_node_type` — PENDIENTE

### ground (1/1)

- [x] `ground` — presente en el GNU

### group_axisymmetric (1/1)

- [x] `group_axisymmetric` — presente en el GNU

### group_beam (3/7)

- [ ] `group_beam_direction_z` — PENDIENTE
- [ ] `group_beam_direction_z_reference_point` — PENDIENTE
- [ ] `group_beam_force_moment_plasti` — PENDIENTE
- [x] `group_beam_inertia` — presente en el GNU
- [x] `group_beam_memory` — presente en el GNU
- [ ] `group_beam_shear` — PENDIENTE
- [x] `group_beam_young` — presente en el GNU

### group_condif (4/5)

- [ ] `group_condif_` — PENDIENTE
- [x] `group_condif_absorption` — presente en el GNU
- [x] `group_condif_capacity` — presente en el GNU
- [x] `group_condif_conductivity` — presente en el GNU
- [x] `group_condif_density` — presente en el GNU

### group_contact (0/8)

- [ ] `group_contact_spring_direction` — PENDIENTE
- [ ] `group_contact_spring_direction_automatic` — PENDIENTE
- [ ] `group_contact_spring_direction_automatic_planes` — PENDIENTE
- [ ] `group_contact_spring_memory` — PENDIENTE
- [ ] `group_contact_spring_plasti_cohesion` — PENDIENTE
- [ ] `group_contact_spring_plasti_friction` — PENDIENTE
- [ ] `group_contact_spring_plasti_friction_automatic` — PENDIENTE
- [ ] `group_contact_spring_sti` — PENDIENTE

### group_dof (0/2)

- [ ] `group_dof_initial` — PENDIENTE
- [ ] `group_dof_initial_speci` — PENDIENTE

### group_ground (0/1)

- [ ] `group_ground` — PENDIENTE

### group_integration (2/3)

- [x] `group_integration_method` — presente en el GNU
- [ ] `group_integration_method_reduced_factor` — PENDIENTE
- [x] `group_integration_points` — presente en el GNU

### group_interface (6/11)

- [x] `group_interface` — Fase 1 implementada (commit `a82cbc7`, 2026-08-13; memory/condif pendientes)
- [ ] `group_interface_condif_conductivity` — PENDIENTE
- [x] `group_interface_gap` — Fase 3 implementada (commit `9c2f4c8`, 2026-08-14; cerrada si strain > gap, default -1e20, hueco físico = gap negativo; validada con iface_mc_gap)
- [ ] `group_interface_ground` — PENDIENTE
- [x] `group_interface_materi_elasti_sti` — Fase 1 implementada (commit `a82cbc7`, 2026-08-13)
- [ ] `group_interface_materi_expansion_normal` — PENDIENTE
- [ ] `group_interface_materi_memory` — PENDIENTE
- [x] `group_interface_materi_plasti_mohr_coul_direct` — Fase 3 implementada (commit `9c2f4c8`, 2026-08-14; MC acumulativo con history `element_interface_force_tang`, activación por presencia del record; validada con iface_mc_slip/iface_mc)
- [x] `group_interface_materi_plasti_tension_direct` — Fase 3 implementada (commit `9c2f4c8`, 2026-08-14; abre en tracción sobre la fuerza normal TOTAL; validada con iface_mc_tension)
- [x] `group_interface_materi_residual_sti` — Fase 3 implementada (commit `60bf78c`, 2026-08-13)
- [ ] `group_interface_tangential_reference_point` — PENDIENTE

### group_materi (56/124)

- [x] `group_materi_damage_mazars` — presente en el GNU
- [x] `group_materi_damping` — presente en el GNU
- [ ] `group_materi_damping_method` — PENDIENTE
- [x] `group_materi_density` — presente en el GNU
- [ ] `group_materi_density_ground` — PENDIENTE
- [ ] `group_materi_elasti_borja_tamagnini` — PENDIENTE
- [ ] `group_materi_elasti_c` — PENDIENTE
- [ ] `group_materi_elasti_c_direction` — PENDIENTE
- [x] `group_materi_elasti_camclay_g` — presente en el GNU
- [x] `group_materi_elasti_camclay_poisson` — presente en el GNU
- [ ] `group_materi_elasti_camclay_pressure_min` — PENDIENTE
- [x] `group_materi_elasti_compressibility` — presente en el GNU
- [ ] `group_materi_elasti_hardsoil` — PENDIENTE
- [ ] `group_materi_elasti_k0` — PENDIENTE
- [x] `group_materi_elasti_lade` — presente en el GNU
- [x] `group_materi_elasti_poisson` — presente en el GNU
- [ ] `group_materi_elasti_poisson_power` — PENDIENTE
- [ ] `group_materi_elasti_shear_factor` — PENDIENTE
- [ ] `group_materi_elasti_stress_pressure_history_factor` — PENDIENTE
- [x] `group_materi_elasti_transverse_isotropy` — presente en el GNU
- [x] `group_materi_elasti_volumetric_poisson` — presente en el GNU
- [x] `group_materi_elasti_volumetric_young_order` — presente en el GNU
- [x] `group_materi_elasti_volumetric_young_values` — presente en el GNU
- [x] `group_materi_elasti_young` — presente en el GNU
- [x] `group_materi_elasti_young_polynomial` — presente en el GNU
- [x] `group_materi_elasti_young_power` — presente en el GNU
- [ ] `group_materi_elasti_young_user` — PENDIENTE
- [x] `group_materi_expansion_linear` — presente en el GNU
- [x] `group_materi_expansion_volume` — presente en el GNU
- [ ] `group_materi_factor` — PENDIENTE
- [ ] `group_materi_failure_crunching` — PENDIENTE
- [x] `group_materi_failure_damage` — presente en el GNU
- [x] `group_materi_failure_plasti_kappa` — presente en el GNU
- [x] `group_materi_failure_rupture` — presente en el GNU
- [ ] `group_materi_failure_void_fraction` — PENDIENTE
- [ ] `group_materi_history_variable_user` — PENDIENTE
- [ ] `group_materi_history_variable_user_parameters` — PENDIENTE
- [x] `group_materi_hyper_besseling` — presente en el GNU
- [x] `group_materi_hyper_blatz_ko` — presente en el GNU
- [x] `group_materi_hyper_mooney_rivlin` — presente en el GNU
- [x] `group_materi_hyper_neohookean` — presente en el GNU
- [ ] `group_materi_hyper_reduced_polynomial` — PENDIENTE
- [x] `group_materi_hyper_volumetric_linear` — presente en el GNU
- [x] `group_materi_hyper_volumetric_murnaghan` — presente en el GNU
- [x] `group_materi_hyper_volumetric_ogden` — presente en el GNU
- [x] `group_materi_hyper_volumetric_polynomial` — presente en el GNU
- [ ] `group_materi_hyper_volumetric_simo_taylor` — PENDIENTE
- [x] `group_materi_maxwell_chain` — presente en el GNU
- [x] `group_materi_membrane` — presente en el GNU
- [x] `group_materi_memory` — presente en el GNU
- [ ] `group_materi_plasti_bounda` — PENDIENTE
- [ ] `group_materi_plasti_bounda_factor` — PENDIENTE
- [x] `group_materi_plasti_camclay` — presente en el GNU
- [ ] `group_materi_plasti_cap1` — PENDIENTE
- [ ] `group_materi_plasti_cap2` — PENDIENTE
- [x] `group_materi_plasti_compression` — presente en el GNU
- [ ] `group_materi_plasti_compression_direct` — PENDIENTE
- [ ] `group_materi_plasti_compression_direct_visco` — PENDIENTE
- [ ] `group_materi_plasti_coord_limit` — PENDIENTE
- [x] `group_materi_plasti_diprisco` — presente en el GNU
- [ ] `group_materi_plasti_diprisco_density` — PENDIENTE
- [ ] `group_materi_plasti_druck_prag` — PENDIENTE
- [ ] `group_materi_plasti_element_group` — PENDIENTE
- [ ] `group_materi_plasti_element_group_factor` — PENDIENTE
- [ ] `group_materi_plasti_generalised_non_associate_cam_clay_for_bonded_soils` — PENDIENTE
- [x] `group_materi_plasti_gurson` — presente en el GNU
- [ ] `group_materi_plasti_hardsoil` — PENDIENTE
- [ ] `group_materi_plasti_heat_generation` — PENDIENTE
- [x] `group_materi_plasti_hypo_cohesion` — presente en el GNU
- [x] `group_materi_plasti_hypo_masin` — implementada (commit `974b284`, 2026-08-06)
- [x] `group_materi_plasti_hypo_masin_clay` — presente en el GNU
- [x] `group_materi_plasti_hypo_masin_clay_advanced_parameters` — presente en el GNU
- [x] `group_materi_plasti_hypo_masin_clay_avanced_direction` — presente en el GNU
- [x] `group_materi_plasti_hypo_masin_clay_ocr` — presente en el GNU
- [x] `group_materi_plasti_hypo_masin_clay_structure` — presente en el GNU
- [x] `group_materi_plasti_hypo_masin_clay_visco` — implementada (commit `a0f55a6`, 2026-08-06)
- [x] `group_materi_plasti_hypo_masin_ocr` — presente en el GNU
- [x] `group_materi_plasti_hypo_masin_structure` — presente en el GNU
- [ ] `group_materi_plasti_hypo_minimum_void_ratio` — PENDIENTE
- [ ] `group_materi_plasti_hypo_niemunis_visco` — PENDIENTE
- [ ] `group_materi_plasti_hypo_niemunis_visco_ocr` — PENDIENTE
- [ ] `group_materi_plasti_hypo_strain_intergranular` — PENDIENTE
- [x] `group_materi_plasti_hypo_strain_intergranular_masin_clay` — presente en el GNU
- [ ] `group_materi_plasti_hypo_strain_isa` — PENDIENTE
- [ ] `group_materi_plasti_hypo_void_ratio_linear` — PENDIENTE
- [ ] `group_materi_plasti_hypo_wol` — PENDIENTE
- [x] `group_materi_plasti_kinematic_hardening` — presente en el GNU
- [ ] `group_materi_plasti_mohr_coul` — PENDIENTE
- [x] `group_materi_plasti_mohr_coul_direct` — implementada (cut-off directo de cortante en un plano; `materi_direct_cutoff` en stress.cc; validada con `materi_direct_mc`)
- [x] `group_materi_plasti_mohr_coul_direct_normal` — implementada (normal explícita del plano; validada con `materi_direct_mc`)
- [x] `group_materi_plasti_mohr_coul_direct_normal_automatic` — implementada (normal del elemento; validada con `materi_direct_auto`)
- [x] `group_materi_plasti_mohr_coul_direct_visco` — implementada (relajación visco: `factor = 1-exp(-dt/tm)`; validada con `materi_direct_visco`)
- [x] `group_materi_plasti_mohr_coul_direct_wall` — implementada (valores alternativos si el elemento está pegado a una pared, vía `plasti_on_boundary`; validada con `materi_direct_wall`)
- [ ] `group_materi_plasti_mohr_coul_hardening_softening` — PENDIENTE
- [ ] `group_materi_plasti_mpc` — PENDIENTE
- [ ] `group_materi_plasti_mpc_factor` — PENDIENTE
- [ ] `group_materi_plasti_pressure_limit` — PENDIENTE
- [ ] `group_materi_plasti_residual_sti` — PENDIENTE
- [x] `group_materi_plasti_tension` — presente en el GNU
- [x] `group_materi_plasti_tension_direct` — implementada (cut-off directo de tracción normal en un plano; validada con `materi_direct`)
- [x] `group_materi_plasti_tension_direct_normal` — implementada (normal explícita; validada con `materi_direct`)
- [x] `group_materi_plasti_tension_direct_normal_automatic` — implementada (normal del elemento; validada con `materi_direct_auto`)
- [x] `group_materi_plasti_tension_direct_visco` — implementada (relajación visco; validada con `materi_direct_visco`)
- [x] `group_materi_plasti_tension_direct_wall` — implementada (valores alternativos en pared; validada con `materi_direct_wall`)
- [x] `group_materi_plasti_user` — presente en el GNU
- [x] `group_materi_plasti_visco_exponential` — presente en el GNU
- [ ] `group_materi_plasti_visco_exponential_limit` — PENDIENTE
- [ ] `group_materi_plasti_visco_exponential_name` — PENDIENTE
- [ ] `group_materi_plasti_visco_exponential_values` — PENDIENTE
- [x] `group_materi_plasti_visco_power` — presente en el GNU
- [ ] `group_materi_plasti_visco_power_name` — PENDIENTE
- [ ] `group_materi_plasti_visco_power_value` — PENDIENTE
- [x] `group_materi_plasti_vonmises` — presente en el GNU
- [x] `group_materi_plasti_vonmises_nadai` — presente en el GNU
- [x] `group_materi_stokes` — presente en el GNU
- [ ] `group_materi_stress_null` — PENDIENTE
- [ ] `group_materi_stress_null_direction` — PENDIENTE
- [ ] `group_materi_stress_null_direction_automatic` — PENDIENTE
- [ ] `group_materi_umat` — PENDIENTE
- [ ] `group_materi_umat_parameters` — PENDIENTE
- [ ] `group_materi_umat_pardiso_decompose` — PENDIENTE
- [ ] `group_materi_undrained_capacity` — PENDIENTE
- [x] `group_materi_viscosity` — presente en el GNU
- [x] `group_materi_viscosity_heatgeneration` — presente en el GNU

### group_plasti (0/1)

- [ ] `group_plasti_apply` — PENDIENTE

### group_porosity (0/1)

- [ ] `group_porosity` — PENDIENTE

### group_spring (2/4)

- [x] `group_spring_direction` — presente en el GNU
- [ ] `group_spring_memory` — PENDIENTE
- [x] `group_spring_plasti` — presente en el GNU
- [ ] `group_spring_sti` — PENDIENTE

### group_time (1/2)

- [x] `group_time` — presente en el GNU
- [ ] `group_time_` — PENDIENTE

### group_truss (5/9)

- [x] `group_truss_area` — presente en el GNU
- [x] `group_truss_density` — presente en el GNU
- [ ] `group_truss_elasti_elongation_force_diagram` — PENDIENTE
- [ ] `group_truss_elasti_young` — PENDIENTE
- [ ] `group_truss_expansion` — PENDIENTE
- [ ] `group_truss_initial_force` — PENDIENTE
- [x] `group_truss_memory` — presente en el GNU
- [x] `group_truss_plasti` — presente en el GNU
- [x] `group_truss_rope` — presente en el GNU

### group_type (1/1)

- [x] `group_type` — presente en el GNU

### group_volume (1/1)

- [x] `group_volume_factor` — presente en el GNU

### group_wave (1/1)

- [x] `group_wave_speed_of_sound` — presente en el GNU

### icontrol (1/1)

- [x] `icontrol` — presente en el GNU

### incremental (0/1)

- [ ] `incremental_driver` — PENDIENTE

### inertia (0/1)

- [ ] `inertia_apply` — PENDIENTE

### input_abaqus (6/6)

- [x] `input_abaqus` — presente en el GNU
- [x] `input_abaqus_continue` — presente en el GNU
- [x] `input_abaqus_group` — presente en el GNU
- [x] `input_abaqus_mesh` — presente en el GNU
- [x] `input_abaqus_name` — presente en el GNU
- [x] `input_abaqus_set` — presente en el GNU

### input_fe (0/1)

- [ ] `input_fe` — PENDIENTE

### input_gmsh (1/1)

- [x] `input_gmsh` — presente en el GNU

### interface (0/1)

- [ ] `interface_gap_apply` — PENDIENTE

### linear (0/1)

- [ ] `linear_calculation_apply` — PENDIENTE

### materi_acceleration (0/1)

- [ ] `materi_acceleration` — PENDIENTE

### materi_damage (1/2)

- [x] `materi_damage` — presente en el GNU
- [ ] `materi_damage_apply` — PENDIENTE

### materi_displacement (1/2)

- [x] `materi_displacement` — presente en el GNU
- [x] `materi_displacement_relative` — implementada (dof relativo `disr*`; referencia re-sincronizada en cambio de timestep en `control_timestep` y en reset de desplazamiento en `control_reset_dof`; integrada en `dof.cc`; validada con `mat_rel` y `mat_rel_reset`)

### materi_dynamic (0/1)

- [ ] `materi_dynamic` — PENDIENTE

### materi_elasti (0/1)

- [ ] `materi_elasti_young_power_apply` — PENDIENTE

### materi_failure (0/1)

- [ ] `materi_failure_apply` — PENDIENTE

### materi_history (0/1)

- [ ] `materi_history_variable` — PENDIENTE

### materi_maxwell (1/1)

- [x] `materi_maxwell_stress` — presente en el GNU

### materi_plasti (4/16)

- [ ] `materi_plasti_camclay_history` — PENDIENTE
- [ ] `materi_plasti_cap1_history` — PENDIENTE
- [ ] `materi_plasti_diprisco_history` — PENDIENTE
- [x] `materi_plasti_f` — presente en el GNU
- [x] `materi_plasti_f_nonlocal` — presente en el GNU
- [ ] `materi_plasti_generalised_non_associate_cam_clay_for_bonded_soils_history` — PENDIENTE
- [ ] `materi_plasti_hardsoil_history` — PENDIENTE
- [ ] `materi_plasti_hypo_history` — PENDIENTE
- [ ] `materi_plasti_hypo_substepping` — PENDIENTE
- [x] `materi_plasti_kappa` — presente en el GNU
- [ ] `materi_plasti_kappa_shear` — PENDIENTE
- [ ] `materi_plasti_maximum_iterations` — PENDIENTE
- [ ] `materi_plasti_phimob` — PENDIENTE
- [x] `materi_plasti_rho` — presente en el GNU
- [ ] `materi_plasti_tension_apply` — PENDIENTE
- [ ] `materi_plasti_visco_apply` — PENDIENTE

### materi_strain (4/21)

- [x] `materi_strain_elasti` — presente en el GNU
- [ ] `materi_strain_energy` — PENDIENTE
- [x] `materi_strain_intergranular` — presente en el GNU
- [ ] `materi_strain_isa_c` — PENDIENTE
- [ ] `materi_strain_isa_eacc` — PENDIENTE
- [x] `materi_strain_plasti` — presente en el GNU
- [ ] `materi_strain_plasti_camclay` — PENDIENTE
- [ ] `materi_strain_plasti_cap` — PENDIENTE
- [ ] `materi_strain_plasti_compression` — PENDIENTE
- [ ] `materi_strain_plasti_diprisco` — PENDIENTE
- [ ] `materi_strain_plasti_druckprag` — PENDIENTE
- [ ] `materi_strain_plasti_generalised_non_associate_cam_clay_for_bonded_soils` — PENDIENTE
- [ ] `materi_strain_plasti_hardsoil` — PENDIENTE
- [ ] `materi_strain_plasti_mohr_coul` — PENDIENTE
- [ ] `materi_strain_plasti_tension` — PENDIENTE
- [ ] `materi_strain_plasti_vonmises` — PENDIENTE
- [x] `materi_strain_total` — presente en el GNU
- [ ] `materi_strain_total_compression_kappa` — PENDIENTE
- [ ] `materi_strain_total_kappa` — PENDIENTE
- [ ] `materi_strain_total_shear_kappa` — PENDIENTE
- [ ] `materi_strain_total_tension_kappa` — PENDIENTE

### materi_stress (1/2)

- [x] `materi_stress` — presente en el GNU
- [ ] `materi_stress_pressure_history` — PENDIENTE

### materi_velocity (2/2)

- [x] `materi_velocity` — presente en el GNU
- [x] `materi_velocity_integrated` — presente en el GNU

### materi_void (1/1)

- [x] `materi_void_fraction` — presente en el GNU

### materi_work (1/1)

- [x] `materi_work` — presente en el GNU

### mesh (1/51)

- [x] `mesh` — presente en el GNU
- [x] `mesh_activate_gravity_element` — implementada (selección de elementos por rango; `mesh_activate_gravity_factor` en mesh.cc)
- [x] `mesh_activate_gravity_element_group` — implementada (selección por grupos)
- [x] `mesh_activate_gravity_geometry` — implementada (selección por geometría)
- [x] `mesh_activate_gravity_method` — implementada (método 1 y 2; método 2 mantiene el elemento activo con rigidez reducida; validada con `mesh_act_grav2`)
- [x] `mesh_activate_gravity_sti` — implementada (factor de rigidez aplicado a la matriz del elemento en materi.cc)
- [x] `mesh_activate_gravity_time` — implementada (activación gradual de la gravedad; validada con `mesh_act_grav`)
- [x] `mesh_activate_gravity_time_initial` — implementada (time of birth)
- [x] `mesh_activate_gravity_time_strain_settlement` — implementada (registrado)
- [ ] `mesh_boundary` — PENDIENTE
- [ ] `mesh_correct` — PENDIENTE
- [ ] `mesh_delete_geometry_moving` — PENDIENTE
- [ ] `mesh_element_group_apply` — PENDIENTE
- [ ] `mesh_gid_arc_coord` — PENDIENTE
- [ ] `mesh_gid_assign_conditions_line` — PENDIENTE
- [ ] `mesh_gid_assign_conditions_point` — PENDIENTE
- [ ] `mesh_gid_assign_conditions_surface` — PENDIENTE
- [ ] `mesh_gid_circle_coord` — PENDIENTE
- [ ] `mesh_gid_circle_element_group` — PENDIENTE
- [ ] `mesh_gid_circle_hollow` — PENDIENTE
- [ ] `mesh_gid_circle_radius` — PENDIENTE
- [ ] `mesh_gid_cylinder_coord` — PENDIENTE
- [ ] `mesh_gid_cylinder_element_group` — PENDIENTE
- [ ] `mesh_gid_cylinder_height` — PENDIENTE
- [ ] `mesh_gid_cylinder_hollow` — PENDIENTE
- [ ] `mesh_gid_cylinder_normal` — PENDIENTE
- [ ] `mesh_gid_cylinder_radius` — PENDIENTE
- [ ] `mesh_gid_line_element_group` — PENDIENTE
- [ ] `mesh_gid_line_point` — PENDIENTE
- [ ] `mesh_gid_line_size` — PENDIENTE
- [ ] `mesh_gid_line_structured_concentrate` — PENDIENTE
- [ ] `mesh_gid_line_structured_nel` — PENDIENTE
- [ ] `mesh_gid_line_structured_size` — PENDIENTE
- [ ] `mesh_gid_point_coord` — PENDIENTE
- [ ] `mesh_gid_rectangle_coord` — PENDIENTE
- [ ] `mesh_gid_rectangle_element_group` — PENDIENTE
- [ ] `mesh_gid_rectangle_hollow` — PENDIENTE
- [ ] `mesh_gid_size` — PENDIENTE
- [ ] `mesh_gid_sphere_coord` — PENDIENTE
- [ ] `mesh_gid_sphere_element_group` — PENDIENTE
- [ ] `mesh_gid_sphere_hollow` — PENDIENTE
- [ ] `mesh_gid_sphere_radius` — PENDIENTE
- [ ] `mesh_gid_surface_element` — PENDIENTE
- [ ] `mesh_gid_surface_element_group` — PENDIENTE
- [ ] `mesh_gid_surface_line` — PENDIENTE
- [ ] `mesh_gid_surface_structured_nel` — PENDIENTE
- [ ] `mesh_gid_surface_structured_size` — PENDIENTE
- [ ] `mesh_gid_volume_element_group` — PENDIENTE
- [ ] `mesh_gid_volume_surface` — PENDIENTE
- [ ] `mesh_interface_triangle_coordinates` — PENDIENTE
- [ ] `mesh_interface_triangle_element_group` — PENDIENTE

### message (0/1)

- [ ] `message` — PENDIENTE

### mpc (0/17)

- [ ] `mpc_apply` — PENDIENTE
- [ ] `mpc_element_group` — PENDIENTE
- [ ] `mpc_element_group_always` — PENDIENTE
- [ ] `mpc_element_group_closest` — PENDIENTE
- [ ] `mpc_element_group_coord_geometry` — PENDIENTE
- [ ] `mpc_element_group_dof` — PENDIENTE
- [ ] `mpc_element_group_eps_iso` — PENDIENTE
- [ ] `mpc_element_group_geometry` — PENDIENTE
- [ ] `mpc_element_group_keep` — PENDIENTE
- [ ] `mpc_geometry` — PENDIENTE
- [ ] `mpc_geometry_dof` — PENDIENTE
- [ ] `mpc_geometry_method` — PENDIENTE
- [ ] `mpc_geometry_switch` — PENDIENTE
- [ ] `mpc_geometry_tolerance` — PENDIENTE
- [ ] `mpc_linear_quadratic` — PENDIENTE
- [ ] `mpc_node_factor` — PENDIENTE
- [ ] `mpc_node_number` — PENDIENTE

### mrange (0/1)

- [ ] `mrange` — PENDIENTE

### mstring (0/1)

- [ ] `mstring` — PENDIENTE

### node (9/23)

- [x] `node` — presente en el GNU
- [x] `node_boundary` — presente en el GNU
- [x] `node_bounded` — presente en el GNU
- [ ] `node_bounded_index` — PENDIENTE
- [ ] `node_convection_apply` — PENDIENTE
- [x] `node_damping` — presente en el GNU
- [x] `node_deformed_mesh` — presente en el GNU
- [x] `node_dof` — presente en el GNU
- [x] `node_dof_calcul` — presente en el GNU
- [ ] `node_dof_start_re` — PENDIENTE
- [ ] `node_dynamic_pressure` — PENDIENTE
- [ ] `node_force` — PENDIENTE
- [ ] `node_geometry_present` — PENDIENTE
- [ ] `node_inertia` — PENDIENTE
- [x] `node_mass` — presente en el GNU
- [ ] `node_mesh` — PENDIENTE
- [x] `node_rhside` — presente en el GNU
- [ ] `node_slide` — PENDIENTE
- [ ] `node_start_re` — PENDIENTE
- [ ] `node_static_pressure` — PENDIENTE
- [ ] `node_sti` — PENDIENTE
- [ ] `node_support_edge_normal_plasti_tension_status` — PENDIENTE
- [ ] `node_total_pressure` — PENDIENTE

### nonlocal (0/2)

- [ ] `nonlocal` — PENDIENTE
- [ ] `nonlocal_name` — PENDIENTE

### number (0/1)

- [ ] `number_of_space_dimensions` — PENDIENTE

### plasti (0/1)

- [ ] `plasti_apply` — PENDIENTE

### post_apply (0/1)

- [ ] `post_apply` — PENDIENTE

### post_calcul (1/22)

- [x] `post_calcul` — presente en el GNU
- [ ] `post_calcul_absolute` — PENDIENTE
- [ ] `post_calcul_apparent_total` — PENDIENTE
- [ ] `post_calcul_label` — PENDIENTE
- [ ] `post_calcul_limit` — PENDIENTE
- [ ] `post_calcul_materi_stress_force_average` — PENDIENTE
- [ ] `post_calcul_materi_stress_force_direction_exclude` — PENDIENTE
- [ ] `post_calcul_materi_stress_force_direction_exclude_epsilon` — PENDIENTE
- [ ] `post_calcul_materi_stress_force_direction_include` — PENDIENTE
- [ ] `post_calcul_materi_stress_force_direction_include_epsilon` — PENDIENTE
- [ ] `post_calcul_materi_stress_force_element_group` — PENDIENTE
- [ ] `post_calcul_materi_stress_force_outer` — PENDIENTE
- [ ] `post_calcul_materi_stress_force_plot_switch` — PENDIENTE
- [ ] `post_calcul_materi_stress_force_reference_point` — PENDIENTE
- [ ] `post_calcul_materi_stress_force_thickness_switch` — PENDIENTE
- [ ] `post_calcul_multiply` — PENDIENTE
- [ ] `post_calcul_safety_default` — PENDIENTE
- [ ] `post_calcul_safety_maximum` — PENDIENTE
- [ ] `post_calcul_safety_method` — PENDIENTE
- [ ] `post_calcul_static_pressure` — PENDIENTE
- [ ] `post_calcul_static_pressure_height` — PENDIENTE
- [ ] `post_calcul_static_pressure_height_element_group` — PENDIENTE

### post_count (0/1)

- [ ] `post_count` — PENDIENTE

### post_data (0/3)

- [ ] `post_data` — PENDIENTE
- [ ] `post_data_factor` — PENDIENTE
- [ ] `post_data_result` — PENDIENTE

### post_element (0/9)

- [ ] `post_element_force` — PENDIENTE
- [ ] `post_element_force_force` — PENDIENTE
- [ ] `post_element_force_geometry` — PENDIENTE
- [ ] `post_element_force_group` — PENDIENTE
- [ ] `post_element_force_inertia` — PENDIENTE
- [ ] `post_element_force_multiply_factor` — PENDIENTE
- [ ] `post_element_force_normal` — PENDIENTE
- [ ] `post_element_force_number` — PENDIENTE
- [ ] `post_element_force_result` — PENDIENTE

### post_global (1/1)

- [x] `post_global` — presente en el GNU

### post_group (0/1)

- [ ] `post_group_volume_summed` — PENDIENTE

### post_integrate (2/2)

- [x] `post_integrate` — presente en el GNU
- [x] `post_integrate_result` — presente en el GNU

### post_line (5/5)

- [x] `post_line` — presente en el GNU
- [x] `post_line_dof` — presente en el GNU
- [x] `post_line_dof_calcul` — presente en el GNU
- [x] `post_line_n` — presente en el GNU
- [x] `post_line_operat` — presente en el GNU

### post_node (4/8)

- [x] `post_node` — presente en el GNU
- [ ] `post_node_factor` — PENDIENTE
- [x] `post_node_result` — presente en el GNU
- [ ] `post_node_rhside_` — PENDIENTE
- [x] `post_node_rhside_free` — presente en el GNU
- [x] `post_node_rhside_ratio` — presente en el GNU
- [ ] `post_node_rhside_ratio_dof_type` — PENDIENTE
- [ ] `post_node_rhside_ratio_method` — PENDIENTE

### post_point (4/6)

- [x] `post_point` — presente en el GNU
- [x] `post_point_dof` — presente en el GNU
- [x] `post_point_dof_calcul` — presente en el GNU
- [ ] `post_point_element_group` — PENDIENTE
- [ ] `post_point_eps_iso` — PENDIENTE
- [x] `post_point_move` — implementada (commit `3793892`, 2026-08-04)

### post_quadrilateral (4/5)

- [x] `post_quadrilateral` — presente en el GNU
- [x] `post_quadrilateral_dof` — presente en el GNU
- [x] `post_quadrilateral_dof_calcul` — presente en el GNU
- [ ] `post_quadrilateral_element_group` — PENDIENTE
- [x] `post_quadrilateral_n` — presente en el GNU

### post_strain (0/3)

- [ ] `post_strain_volume_absolute` — PENDIENTE
- [ ] `post_strain_volume_initial` — PENDIENTE
- [ ] `post_strain_volume_relative` — PENDIENTE

### print_ (0/1)

- [ ] `print_` — PENDIENTE

### print_apply (0/1)

- [ ] `print_apply` — PENDIENTE

### print_arithmetic (1/1)

- [x] `print_arithmetic` — presente en el GNU

### print_control (1/1)

- [x] `print_control` — presente en el GNU

### print_data (0/1)

- [ ] `print_data_name` — PENDIENTE

### print_database (0/1)

- [ ] `print_database_calculation` — PENDIENTE

### print_de (0/1)

- [ ] `print_de` — PENDIENTE

### print_element (0/2)

- [ ] `print_element_geometry_present` — PENDIENTE
- [ ] `print_element_geometry_present_node_type` — PENDIENTE

### print_failure (1/1)

- [x] `print_failure` — presente en el GNU

### print_frd (0/2)

- [ ] `print_frd_freecad` — PENDIENTE
- [ ] `print_frd_prepomax` — PENDIENTE

### print_gid (0/10)

- [ ] `print_gid_calculation` — PENDIENTE
- [ ] `print_gid_contact_spring2` — PENDIENTE
- [ ] `print_gid_coord` — PENDIENTE
- [ ] `print_gid_de` — PENDIENTE
- [ ] `print_gid_group` — PENDIENTE
- [ ] `print_gid_mesh_activate_gravity` — PENDIENTE
- [ ] `print_gid_node_method` — PENDIENTE
- [ ] `print_gid_smooth_dof` — PENDIENTE
- [ ] `print_gid_smooth_n` — PENDIENTE
- [ ] `print_gid_spring2` — PENDIENTE

### print_gmsh (0/3)

- [ ] `print_gmsh_calculation` — PENDIENTE
- [ ] `print_gmsh_dummy` — PENDIENTE
- [ ] `print_gmsh_node_method` — PENDIENTE

### print_group (0/1)

- [ ] `print_group_data` — PENDIENTE

### print_mesh (0/1)

- [ ] `print_mesh_dof` — PENDIENTE

### print_node (0/2)

- [ ] `print_node_geometry_present` — PENDIENTE
- [ ] `print_node_geometry_present_node_type` — PENDIENTE

### print_precision (0/1)

- [ ] `print_precision` — PENDIENTE

### print_tecplot (0/1)

- [ ] `print_tecplot_calculation` — PENDIENTE

### print_vtk (0/4)

- [ ] `print_vtk_calculation` — PENDIENTE
- [ ] `print_vtk_coord` — PENDIENTE
- [ ] `print_vtk_group` — PENDIENTE
- [ ] `print_vtk_node_method` — PENDIENTE

### print_where (1/1)

- [x] `print_where` — presente en el GNU

### processors (0/3)

- [ ] `processors` — PENDIENTE
- [ ] `processors_maximum` — PENDIENTE
- [ ] `processors_partition` — PENDIENTE

### repeat (0/2)

- [ ] `repeat_save_calculate_result` — PENDIENTE
- [ ] `repeat_save_result` — PENDIENTE

### safety (0/31)

- [ ] `safety_slip_circle_grid_middle` — PENDIENTE
- [ ] `safety_slip_circle_grid_middle_n` — PENDIENTE
- [ ] `safety_slip_circle_grid_radius` — PENDIENTE
- [ ] `safety_slip_circle_grid_radius_n` — PENDIENTE
- [ ] `safety_slip_circle_grid_result` — PENDIENTE
- [ ] `safety_slip_circle_grid_segment_n` — PENDIENTE
- [ ] `safety_slip_circle_line_middle` — PENDIENTE
- [ ] `safety_slip_circle_line_middle_n` — PENDIENTE
- [ ] `safety_slip_circle_line_radius` — PENDIENTE
- [ ] `safety_slip_circle_line_radius_n` — PENDIENTE
- [ ] `safety_slip_circle_line_result` — PENDIENTE
- [ ] `safety_slip_circle_line_segment_n` — PENDIENTE
- [ ] `safety_slip_combined_linear` — PENDIENTE
- [ ] `safety_slip_combined_linear_n` — PENDIENTE
- [ ] `safety_slip_combined_linear_result` — PENDIENTE
- [ ] `safety_slip_combined_linear_segment_n` — PENDIENTE
- [ ] `safety_slip_ellipsoide` — PENDIENTE
- [ ] `safety_slip_ellipsoide_method` — PENDIENTE
- [ ] `safety_slip_ellipsoide_n` — PENDIENTE
- [ ] `safety_slip_ellipsoide_result` — PENDIENTE
- [ ] `safety_slip_ellipsoide_segment_n` — PENDIENTE
- [ ] `safety_slip_grd` — PENDIENTE
- [ ] `safety_slip_grd_method` — PENDIENTE
- [ ] `safety_slip_grd_method_direction` — PENDIENTE
- [ ] `safety_slip_grd_segment_n` — PENDIENTE
- [ ] `safety_slip_multi_linear` — PENDIENTE
- [ ] `safety_slip_multi_linear_n` — PENDIENTE
- [ ] `safety_slip_multi_linear_result` — PENDIENTE
- [ ] `safety_slip_multi_linear_segment_n` — PENDIENTE
- [ ] `safety_slip_set` — PENDIENTE
- [ ] `safety_slip_set_result` — PENDIENTE

### slide (1/7)

- [x] `slide_axisymmetric` — implementada (escala la fricción de slide por `2*pi*r`, r = coordenada radial del nodo en axisimétrico; validada con `slide_axi`)
- [ ] `slide_damping` — PENDIENTE
- [x] `slide_geometry` — presente en el GNU
- [ ] `slide_plasti_friction` — PENDIENTE
- [ ] `slide_plasti_tension` — PENDIENTE
- [ ] `slide_sti` — PENDIENTE
- [ ] `slide_user` — PENDIENTE

### solver (0/10)

- [ ] `solver` — PENDIENTE
- [ ] `solver_bicg_error` — PENDIENTE
- [ ] `solver_bicg_restart` — PENDIENTE
- [ ] `solver_bicg_stop` — PENDIENTE
- [ ] `solver_matrix_save` — PENDIENTE
- [ ] `solver_matrix_symmetric` — PENDIENTE
- [ ] `solver_pardiso_ordering` — PENDIENTE
- [ ] `solver_pardiso_out_of_core` — PENDIENTE
- [ ] `solver_pardiso_processors` — PENDIENTE
- [ ] `solver_pardiso_processors_maximum` — PENDIENTE

### strain (0/10)

- [x] `strain_settlement_diagram` — implementada (dependencia de parámetro del creep en un dof vía tabla; validada con `strain_settle_diag`)
- [x] `strain_settlement_diagram_dof` — implementada
- [x] `strain_settlement_diagram_number` — implementada
- [x] `strain_settlement_element_group` — implementada
- [x] `strain_settlement_parameters` — implementada (creep de asentamiento con saturación; `strain_settlement_creep` en materi.cc; validada con `strain_settle`)
- [ ] `strain_volume_absolute_time` — PENDIENTE
- [ ] `strain_volume_element` — PENDIENTE
- [ ] `strain_volume_element_group` — PENDIENTE
- [ ] `strain_volume_geometry` — PENDIENTE
- [ ] `strain_volume_relative_time` — PENDIENTE

### support (0/18)

- [ ] `support_edge_normal` — PENDIENTE
- [ ] `support_edge_normal_damping` — PENDIENTE
- [ ] `support_edge_normal_damping_automatic` — PENDIENTE
- [ ] `support_edge_normal_damping_automatic_apparent` — PENDIENTE
- [ ] `support_edge_normal_density` — PENDIENTE
- [ ] `support_edge_normal_element_group` — PENDIENTE
- [ ] `support_edge_normal_element_node` — PENDIENTE
- [ ] `support_edge_normal_element_side` — PENDIENTE
- [ ] `support_edge_normal_factor` — PENDIENTE
- [ ] `support_edge_normal_force_initial` — PENDIENTE
- [ ] `support_edge_normal_geometry` — PENDIENTE
- [ ] `support_edge_normal_node` — PENDIENTE
- [ ] `support_edge_normal_plasti_compression` — PENDIENTE
- [ ] `support_edge_normal_plasti_friction` — PENDIENTE
- [ ] `support_edge_normal_plasti_residual_sti` — PENDIENTE
- [ ] `support_edge_normal_plasti_tension` — PENDIENTE
- [ ] `support_edge_normal_plasti_tension_double` — PENDIENTE
- [ ] `support_edge_normal_time` — PENDIENTE

### target (2/2)

- [x] `target_item` — presente en el GNU
- [x] `target_value` — presente en el GNU

### time (2/2)

- [x] `time_calculation` — presente en el GNU
- [x] `time_current` — presente en el GNU

### timestep (0/2)

- [ ] `timestep_iterations_automatic_apply` — PENDIENTE
- [ ] `timestep_predict_velocity` — PENDIENTE

### tochnog (0/1)

- [ ] `tochnog_version` — PENDIENTE

### truss (0/1)

- [ ] `truss_rope_apply` — PENDIENTE

### up (0/1)

- [ ] `up` — PENDIENTE

### volume (1/2)

- [x] `volume_factor` — presente en el GNU
- [ ] `volume_factor_x` — PENDIENTE

### wave (2/2)

- [x] `wave_fscalar` — presente en el GNU
- [x] `wave_scalar` — presente en el GNU

### zip (0/1)

- [ ] `zip` — PENDIENTE

</details>
