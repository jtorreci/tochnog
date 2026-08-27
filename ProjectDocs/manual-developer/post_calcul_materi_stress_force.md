# post_calcul_materi_stress_force

## Implementación

- **Lot 1 (infrastructure, commit cc0aba1)**: registration of the
  `post_calcul_materi_stress_force_*` records + the FORCE branch in
  `calculate()`/`calculate_operat()` (`calcul.cc`) + the item
  generation (9 items 2D / 16 items 3D in `NODE_DOF_CALCUL`) + the
  print (`print_materi_stress_force.cc`). See the per-record pages.
- **Lot 2 (this lot, 2D numerical integration)** in `calcul_force.cc`:
  `post_calcul_materi_stress_force()` (the per-node entry, called from
  `calculate_operat` in the parallel node loop) → `msf_calculate_node_2d()`
  (scans the target element groups, accumulates the per-node
  contributions) → `msf_element_contribution_2d()` (per element: end
  faces, integration, node role, plot components) →
  `msf_integrate_side_2d()` (the 1D stress integration over one face).
  The per-node averaged flag
  `POST_CALCUL_MATERI_STRESS_FORCE_AVERAGED_NODE` (INTEGER, NODE
  class, version_all=1, registered in `database.cc`, allocated/deleted
  in `calculate()` like `NODE_DOF_CALCUL`) feeds the `-primary` filter
  (`msf_node_is_averaged()` in `print_materi_stress_force.cc`).
  Enum additions in `tochnog.h`/`tochnog-mod.h` in sync (clean build
  required after enum changes).

## Diseño / decisiones (con evidencia; el test analítico es el árbitro)

1. **Caras extremas (end faces)**: la dirección de espesor de la
   sección es t̂ = (centroide del elemento − reference_point)
   normalizada EN EL PLANO (manual 6.914: el reference_point define
   fuera/dentro en dirección de espesor). Las dos caras extremas = las
   2 aristas cuya normal exterior es MÁS PERPENDICULAR a t̂
   (|n·t̂| mínima): las secciones donde actúan las fuerzas. Las otras 2
   aristas (superficies de la estructura) no producen valores
   primarios. Selección ambigua (elemento distorsionado o
   reference_point sobre la diagonal: |n·t̂| iguales) → aviso único y
   el elemento se omite (valores 0).
2. **Fuente de σ**: la tensión NODAL (`node_dof[stres_indx +
   stress_indx(i,j)*nder]`), la incógnita resuelta. Directa y exacta
   para el campo lineal/cuadrático a lo largo de la arista; la
   alternativa `ELEMENT_DOF` (tensiones de los puntos de integración)
   requeriría extrapolación IP→arista y se descartó. Verificado:
   `materi_stress` debe estar en initia (validado).
3. **Cuadratura de arista (1D)**: `integration_gauss(2)` para npol=2
   (quad4) y `integration_lobatto(3)` para npol=3 (quad9). El
   integrando del MOMENTO σ_nn·dt es cuadrático (quad4) / cúbico
   (quad9): Gauss(2) es exacto para grado ≤3 y Lobatto(3) (Simpson)
   para grado ≤3 → el momento se integra EXACTAMENTE para campos
   lineales/cuadráticos. NOTA de precisión: la nota del plan "npol=2 ≡
   Gauss ±1/√3" es correcta para la arista (el trapezio Lobatto(2)
   sobreestimaría el momento 3/2, verificado analíticamente); los
   pesos de Tochnog suman 1 (integral = longitud·Σ w·f).
4. **Momento**: mom = ∫σ_nn·dt ds con dt = (x_q − centroide)·t̂ — la
   distancia EN LA DIRECCIÓN DE ESPESOR respecto al centro del
   elemento (manual 6.913: "a distance in thickness direction dt
   relative to the middle of the element"). La parametrización de la
   propia arista s daría signos OPUESTOS en las dos caras de la misma
   sección (una recorre de abajo a arriba y la otra al revés) —
   verificado: con (s−s_mid) los nodos del plano medio salían
   (f1−f2)/2 en vez del promedio. dt = (x−C)·t̂ hace el momento
   consistente (ambas caras −M para la ménsula; moms = P·(8−x) medido
   dentro del 1-3%).
5. **Valores** (por unidad de longitud l; 2D plano l=1, axisimétrico
   l=2π·r con r = coordenada radial del centroide, convención de
   area.cc): nor = ∫σ_nn ds (CON SIGNO, tracción +), she = |∫σ_nt ds|
   (solo tamaño, manual), mom = ∫σ_nn·dt ds (con signo). Componentes
   de plot: norx/nory = nor·t̂, nors = |nor|; shex/shey = |she|·t̂,
   shes = |she|; momx/momy = mom·t̂, moms = |mom|. El signo
   tracción/compresión de nor y mom viaja en la DIRECCIÓN del vector
   de plot (t̂ hacia/desde el reference_point); la componente s es el
   tamaño físico (manual 6.913: "the size of the vector formed by
   these components indeed is the real physical size").
6. **Asignación por nodo**: los nodos de las caras extremas son
   PRIMARIOS (reciben el valor de su cara; en una malla conforme los
   nodos compartidos reciben la media de contribuciones idénticas).
   Nodos del plano medio (quad9, los 3 que NO están en las caras):
   con average -yes (default) reciben el PROMEDIO de las dos caras
   (verificado EXACTO: msf_beam2d/msf_quad9) y se marcan en
   `POST_CALCUL_MATERI_STRESS_FORCE_AVERAGED_NODE`; con -no no reciben
   nada (0) y no se marcan (msf_quad9_noavg: -all == -primary).
7. **outer -yes**: solo los nodos PRIMARIOS a máxima distancia del
   reference_point reciben valores (manual 6.915); los promediados no
   reciben nada (decisión documentada). Implementado en 2D.
8. **plot_switch**: -yes invierte las componentes x/y del item
   (dirección de dibujo del vector, manual 6.916); la componente s
   (tamaño) no cambia. Implementado en 2D.

## Evidencia / precisión FE (medida en los tests)

| modelo | magnitud | esperado | medido |
|--------|----------|----------|--------|
| msf_beam2d (8×quad9, L=8, P=1e-2) | moms en x=0..8 | P·(8−x) | 0.0820, 0.0708, 0.0600, 0.0500, 0.0400, 0.0300, 0.0200, 0.0099, 0.0003 (0-2.5%) |
| msf_beam2d | promediado x=0.5 | (f0+f1)/2 | 0.0764 = (0.0820+0.0708)/2 EXACTO |
| msf_beam2d | shes | P (banda ±50%) | 0.0109..0.0133 (polución FE) |
| msf_beam2d_pure (4 puntos) | moms tramo central | 1e-2 cte | 0.01015 (1.5%) |
| msf_beam2d_pure | shes tramo central | 0 | 3e-7 (≈0 EXACTO) |
| msf_quad9 | -all vs -primary | 27 vs 15 | 27 vs 15 |
| msf_nor (8×quad4 axial, N=2e-2) | nors | 2e-2 | 0.0209 (x=0, Poisson) y 0.0200 (0-5%) |
| msf_shear (corte simple) | shes | G·γ·h = 0.384615 | 0.384615 EXACTO |

## Detalles / gotchas

- **GOTCHA (cazado con AddressSanitizer)**: leer las coordenadas
  nodales con `i<MDIM` (3) en vez de `i<ndim` (2) lee 8 bytes FUERA
  del record NODE (heap buffer overflow READ) en el último nodo; el
  valor basura corrompía el solve posterior (el campo de tensiones
  salía LOCALIZADO en las últimas columnas en el axial) por
  interacción con el asignador. Corregido (`i<ndim`). Lección:
  los records NODE tienen `data_length = ndim`, NO MDIM.
- **GOTCHA (solver mixto iterativo)**: con `nu=0.0` y mallas
  multi-elemento el solve BiCGStab con las iteraciones por defecto (2)
  LOCALIZA la deformación axial (σ_xx = 0 en la mitad izquierda);
  con `nu=0.3` (o más iteraciones) el campo es uniforme. No es del
  integrador: los tests usan nu=0.3 (documentado en msf_nor).
- **GOTCHA (par de fuerzas nodales)**: un PAR de fuerzas en la punta
  NO excita flexión pura en el FE discreto: su trabajo sobre el modo
  de flexión es CERO (u_top = u_bottom en la misma sección) y el
  solver responde con cizalla pura (σ_xy cte, σ_xx≈0) — verificado en
  quad4 Y quad9. Por eso el test de flexión pura usa la configuración
  estándar de 4 puntos (viga apoyada + 2 cargas → tramo central con M
  cte y V=0).
- **GOTCHA (malla 1-en-espesor)**: la cizalla integrada de la
  formulación mixta tiene polución en las superficies libres
  (condición de superficie libre débil): |∫σ_xy dy| oscila ±30% de P
  según la sección (el valor exacto se verifica con msf_shear, campo
  uniforme). El quad4 1-en-espesor sufre shear locking (mom ≈ 0.23× y
  she ≈ 0.74× del valor estático, medido) → los tests de flexión usan
  quad9.
- El flag averaged se pre-aloca en `calculate()` (db_allocate) como
  NODE_DOF_CALCUL: el PUT dentro del bucle paralelo no puede alocar.
  version_all=1 → db_version_copy + renumbering lo llevan a
  VERSION_PRINT para el print (como NODE_DOF_CALCUL).
- `validate()` añade en L2: `materi_stress` obligatorio en initia
  (stres_indx ≥ 0) y, en 2D, los grupos objetivo solo pueden contener
  quad4/quad9 (error claro para otros tipos). El aviso "not yet
  implemented" queda SOLO para 3D (lote 3).

## Pendiente

- 3D (hex8/hex27): selección de caras por direction_exclude/include,
  mom1/mom2 (distancias dt y dl), thickness_switch — lote 3.
- axisimétrico: implementado (l=2πr) pero sin test dedicado en este
  lote (pendiente).
- Los tests `.dat` viven en `validation-suite/test-2014/` (gitignored;
  solo el bucle y los checks de `scripts/build_safe.sh` se versionan).
