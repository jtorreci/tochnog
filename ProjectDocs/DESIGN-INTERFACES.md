# Carril A — Elementos de interfaz (diseño técnico)

Fecha: 2026-08-13. Estado: DISEÑO (fase 1 pendiente de implementar).

## Objetivo

Implementar la familia de elementos de interfaz de Tochnog Professional:
elementos que modelan juntas/discontinuidades entre bloques de material
(p.ej. entre un pilote y el suelo). La familia cubre el elemento en sí,
su ley constitutiva (`group_interface_*`), la conversión automática
(`control_mesh_convert`), y el post-proceso
(`control_print_interface_stress*`).

## Modelo físico

En un elemento de interfaz las **deformaciones son diferencias de
desplazamiento entre los dos lados opuestos** del elemento (no gradientes
de campo). Con coordenadas locales (n, t1, t2):

- `strain_normal = (u_n^side2 - u_n^side1) / h` (o simplemente la
  diferencia de desplazamiento normal si h = 1).
- `strain_shear1 = (u_t1^side2 - u_t1^side1) / h`.
- `strain_shear2 = (u_t2^side2 - u_t2^side1) / h` (3D).

El elemento `-bar2` (2 nodos) es el "lado 1"; al convertir a `-quad4`
(4 nodos), los nodos 3-4 forman el "lado 2". La conversión añade nodos y
los conecta a los elementos isoparamétricos vecinos.

## Fases de implementación

### Fase 1 — Elemento de interfaz + rigidez elástica (kn, kt)

- Nuevo tipo de elemento tratado en `elem.cc`: cuando el grupo tiene
  `group_interface -yes`, el elemento se procesa con la ley de interfaz.
- `group_interface_materi_elasti_stiffness kn kt,first kt,second`:
  - `stress_normal = kn * strain_normal`
  - `stress_shear1 = kt,first * shear_gamma1` (shear_gamma = 2*strain)
  - `stress_shear2 = kt,second * shear_gamma2`
- La matriz del elemento y el rhs se ensamblan con estas tensiones.
- Elementos 2D soportados primero: `-quad4` (de bar2) y `-quad6`
  (de bar3/quad8/quad9).

### Fase 2 — Conversión automática

- `control_mesh_convert index switch`: bar2→quad4, bar3→quad6, tria3→
  prism6, quad4→hex8, etc. cuando el grupo es interfaz.
- `control_mesh_convert_element_group`: ayuda a conectar los lados
  correctos (grupos a un lado de la interfaz).
- `control_mesh_generate_interface_geometry`.

### Fase 3 — Ley constitutiva completa de interfaz

- `group_interface_materi_elasti_stiffness` (rigideces kn/kt).
- `group_interface_materi_plasti_mohr_coul_direct phi c phi_flow`:
  max fricción = c + Fn*tan(phi).
- `group_interface_materi_plasti_tension_direct tension_limit`.
- `group_interface_gap gap` (espacio inicial; cerrar solo si la
  deformación normal < gap).
- `group_interface_materi_residual_stiffness factor`.
- `group_interface_materi_memory -updated_linear | -total_linear`.
- `group_interface_condif_conductivity k`.
- `group_interface_groundflow_capacity C`, `_permeability pe`,
  `_total_pressure_tension`.
- `group_interface_materi_expansion_normal`.
- `group_interface_tangential_reference_point`.

### Fase 4 — Post-proceso

- `control_print_interface_stress*`: tensiones de interfaz.
- `control_reset_interface*`, `control_reset_interface_strain`.

## Consideraciones de validación

- **Sin tests de referencia**: la suite sfnet no tiene tests de interfaz
  (solo contac*, que es contact spring). No hay ejemplo numérico completo
  en el manual (el ejemplo bar2→hex8 está incompleto con "...").
- **Validación propuesta**: construir un test propio 2D con un quad4 de
  interfaz entre dos bloques, con `group_interface -yes` + `_elasti_
  stiffness`, y verificar que (a) con kn alta la interfaz transmite el
  desplazamiento (como si fuera continua), (b) con kn baja hay salto de
  desplazamiento, (c) las tensiones siguen stress=kn*strain.
- Verificar contra el comportamiento del mismo modelo SIN interfaz
  (bloques conectados) para el límite kn→∞.

## Dependencias y decisiones

- El elemento de interfaz se integra en `elem.cc` (núcleo del solver),
  no como módulo aparte — es donde se ensambla la rigidez.
- Requiere entender cómo `elem()` construye `element_lhside` /
  `element_matrix` para isoparamétricos y duplicar el patrón para la
  ley de interfaz.
- `group_interface` es un keyword nuevo (data_class GROUP, tipo INTEGER).

## Estrategia numérica de no-interpenetración (decidida 2026-08-13)

**Opción elegida: penalización implícita + iteraciones** (patrón de
tochnog, validado en la Fase 1).

- La interfaz elástica usa `kn`/`kt` como penalización. Con `kn` alta los
  lados quedan acoplados (límite = bloques soldados); con `kn` baja
  deslizan.
- El solver resuelve **implícito** (Newton) con
  `control_timestep_iterations`: la interpenetración se corrige
  iterativamente dentro del paso, no por subpasos.
- En la Fase 3 (unilateralidad: `group_interface_gap`,
  `_tension_direct`, `_mohr_coul_direct`), la rigidez se **actualiza
  dentro de las iteraciones**: si la interfaz se abre (tracción o
  deformación normal > gap), `kn` → rigidez residual
  (`_residual_stiffness`); si se cierra, `kn` se reactiva. Es el mismo
  patrón que `group_contactspring` en conspr.cc.
- **No se usa** line-search ni arc-length: el Newton implícito maneja la
  no-linealidad de la transición cerrado/abierto. Line-search se
  reconsideraría solo si la convergencia falla en la Fase 3.

### Validación de la Fase 1 (2026-08-13)

Test de 2 bloques (`iface2.dat`): bloque izquierdo fijo, derecho empujado,
conectados por la interfaz.

| kn | velix(nodo6) | Comportamiento |
|----|-------------|----------------|
| 100 | 0.044 | acopla parcialmente |
| 1000 | 0.0018 | casi soldado |
| 1e6 | -0.0032 | ≈ soldado |
| 0.001 | ≈1.0 | deslizamiento libre |
| (sin interfaz) | -0.0066 | referencia soldado |

Límites correctos: kn→∞ → soldado, kn→0 → libre. La fuerza nodal usa
`-sign*stress*dir` (principio de trabajos virtuales); el signo invertido
hacía que la interfaz empujara en vez de resistir.

## Fase 2 — Conversión automática (diseño 2026-08-13)

### Objetivo

`control_mesh_convert` convierte automáticamente los elementos de interfaz
de baja dimensión a su equivalente isoparamétrico, creando los nodos del
lado opuesto de la interfaz. Caso principal 2D: `-bar2` → `-quad4`.

### Flujo de uso (del manual de Professional)

```
element 1 -bar2 101 102
element_group 1 10
group_interface 10 -yes
control_mesh_convert 110 -yes
control_mesh_convert_element_group 110 0 1   (grupos a un lado de la interfaz)
```

El usuario genera con GID una malla con elementos `-bar2` en la interfaz
(entre pile y soil), y `control_mesh_convert` crea los nodos duplicados y
reconecta los vecinos.

### Algoritmo (caso 2D bar2 → quad4)

Para cada elemento `-bar2` cuyo grupo tiene `group_interface -yes`:

1. **Leer el bar2**: nodos {a, b} — forman el lado 1 de la interfaz.
2. **Calcular la normal** de la interfaz: perpendicular a la línea a-b,
   en el plano de la malla (2D).
3. **Crear 2 nodos nuevos** {a', b'} = copias de {a, b} desplazadas en la
   normal (dirección hacia el otro lado). Copiar NODE, NODE_START_REFINED,
   NODE_DOF, NODE_DOF_START_REFINED (patrón generate.cc:340-390).
4. **Reescribir el elemento** como `-quad4` con nodos {a, b, a', b'}.
5. **Reconectar los vecinos**: para cada elemento vecino que comparte los
   nodos {a, b}:
   - Si el vecino está en un grupo de `control_mesh_convert_element_group`
     (un lado de la interfaz): se queda con {a, b}.
   - Si el vecino está en el OTRO lado: reemplazar {a, b} por {a', b'} en
     su conectividad.
6. **Reconstruir la malla**: `mesh_has_changed(VERSION_NORMAL)` para
   actualizar NODE_NODE/NODE_ELEMENT.

### Identificación de lados

- `control_mesh_convert_element_group index g0 g1 ...` define los grupos a
  UN lado de la interfaz. Los vecinos en esos grupos usan los nodos
  originales; los demás vecinos (que comparten el lado) usan los nuevos.
- La normal apunta desde el lado especificado hacia el otro.

### Condiciones de validez (del manual)

- Cada interfaz debe tener vecinos isoparamétricos que compartan un lado
  completo.
- Las superficies con interfaces no deben intersectarse.

### Invocación en top.cc

La conversión se inserta en `mesh_changed()` (o tras `generate_spring`),
igual que las otras generaciones: comprobar `CONTROL_MESH_CONVERT` activo
y llamar a `interface_convert()`.

### Casos 3D (futuro)

- `-tria3` → `-prism6`, `-quad4` → `-hex8` (3D), etc. Siguen el mismo
  patrón con más nodos nuevos.

### Riesgos

- La reconexión de vecinos es la parte delicada: identificar correctamente
  qué vecinos están a cada lado sin doble-conectar.
- Añadir nodos cambia la numeración; requiere renumbering consistente.
- Los dofs de los nodos nuevos deben inicializarse (copiar del original).

## Riesgos

- `elem.cc` es complejo (1007 líneas) y toca el ensamblaje global; un
  error en la interfaz puede corromper la solución de TODOS los elementos.
- La conversión (`control_mesh_convert`) implica re-mallado (añadir
  nodos) — infraestructura de renumbering sensible.

## Análisis del codebase (2026-08-13) — cómo encaja la interfaz

### Estructura del ensamblaje

`elem()` (elem.cc, 1007 líneas) procesa cada elemento:
1. Lee `new_dof`/`old_dof` (nuknwn valores por nodo) y construye
   `coord`/`new_coord`.
2. Para elementos discretos (spring/truss/beam) llama a funciones
   dedicadas: `spring()`, `truss()`, `beam_3d()` (mismo patrón).
3. Para isoparamétricos: `pol()` calcula funciones de forma `h`,
   gradientes `old_d`/`new_d`, volumen; luego por cada punto de
   integración llama a `general()` (general.cc) que calcula tensiones
   vía `materi_*` y ensambla `element_lhside`/`element_rhside`/
   `element_matrix`.

### Patrón de referencia: `spring()` (spring.cc)

La interfaz es un **spring2 generalizado**: su deformación es la
diferencia de desplazamiento entre los 2 lados. `spring()` muestra el
patrón exacto:

- Firma: `spring(element, name, element_group, coord, old_dof, new_dof,
  element_lhside, element_matrix, element_rhside)`.
- Lee `GROUP_SPRING_STIFFNESS` del grupo; calcula la longitud
  incremental; ensambla fuerzas y matriz con `fac=±1` por nodo.
- Escribe `ELEMENT_SPRING_FORCE` (history) con `PUT` en VERSION_NEW.

### Estrategia de implementación de la Fase 1

1. **Nuevo archivo `interface.cc`** con `interface_element(element, name,
   element_group, coord, old_dof, new_dof, element_lhside, element_matrix,
   element_rhside)` — copia el patrón de `spring()`.

2. **Geometría de la interfaz (2D, cuadrilátero)**: el elemento tiene 4
   nodos (2 por lado, convertidos de bar2). Lados: nodos {0,1} = lado 1,
   nodos {2,3} = lado 2.
   - Vector normal `n` (perpendicular a la línea media).
   - Vector tangencial `t` (a lo largo de la interfaz).
   - Diferencia de desplazamiento: `du = u^lado2 - u^lado1` por nodo.

3. **Deformaciones** (por punto de integración, usando `h` de `pol()`):
   - `strain_normal = (du · n) / h` (h = espesor, o 1 si h se absorbe en
     kn).
   - `strain_shear = (du · t) / h`.

4. **Tensiones** (elástico, `GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS`):
   - `stress_normal = kn * strain_normal`.
   - `stress_shear = kt * 2 * strain_shear` (shear gamma = 2*strain).
   - Aplicar tensión de tracción solo si la interfaz está cerrada
     (`group_interface_gap`), con rigidez residual (`_residual_stiffness`)
     si está abierta (fase 3).

5. **Ensamblaje**: para cada nodo, `element_rhside[i] -= stress *
   du/ddoflow` y la matriz con la derivada (kn*dtime en el patrón de
   spring). La matriz es `[K -K; -K K]` en las DOF de desplazamiento
   (vel_indx para velocity_integrated, o dis_indx).

6. **Conexión en elem.cc**: en `elem()`, antes del procesamiento
   isoparamétrico estándar, si `GROUP_INTERFACE` del grupo es `-YES`,
   llamar a `interface_element()` y saltar el resto (igual que
   spring/truss).

### Enums y keywords de la Fase 1

- Enums nuevos: `GROUP_INTERFACE`, `GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS`
  (en la familia GROUP_*).
- Keywords: `group_interface index switch`, `group_interface_materi_elasti_
  stiffness index kn kt,first kt,second` (data_class GROUP, DOUBLE).
- El elemento convertido (quad4 de interfaz) es un `-quad4` normal con el
  flag de grupo; no necesita un nuevo tipo de elemento.

### Siguientes pasos de diseño pendientes

- **Representación en pol()**: el quad4 de interfaz (de bar2) tiene 4
  nodos pero su deformación es la diferencia entre 2 lados. Opciones:
  (a) usar `pol()` con `-QUAD4` y en `interface_element()` ignorar los
  gradientes, usando solo `h` y las coordenadas de los 4 nodos para
  calcular la normal media y las diferencias de desplazamiento por lado;
  (b) tratar la interfaz como un `-BAR2` "engrosado". La opción (a) es
  más limpia porque reutiliza las funciones de forma de pol y la
  integración.
- **Dof de deformación**: usar `vel_indx` (velocity_integrated) o
  `dis_indx` (displacement), como spring.cc (veli si velocity_integrated,
  dis si displacement).
- **Matriz de rigidez**: `[K -K; -K K]` bloqueada en las DOF de
  desplazamiento de los 2 lados, con K = kn (normal) y kt (tangencial)
  proyectadas sobre la normal/tangente.
- **Test de validación** (diseñar antes de codificar):
  - Modelo: 2 bloques (quad4) conectados por una interfaz (quad4 de
    interfaz entre ellos), cargado.
  - Verificación 1 (kn alta, kt alta): el comportamiento tiende al de los
    bloques conectados directamente (sin interfaz) — continua.
  - Verificación 2 (kn baja): salto de desplazamiento normal entre lados.
  - Verificación 3 (kt baja): deslizamiento tangencial.
  - Verificación 4: stress_normal = kn * strain_normal numéricamente.
  - Este test propio es la única validación (no hay referencia sfnet).
