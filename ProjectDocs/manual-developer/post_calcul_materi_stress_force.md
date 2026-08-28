# post_calcul_materi_stress_force

## Implementación

- **Lot 1 (infrastructure, commit cc0aba1)**: registration of the
  `post_calcul_materi_stress_force_*` records + the FORCE branch in
  `calculate()`/`calculate_operat()` (`calcul.cc`) + the item
  generation (9 items 2D / 16 items 3D in `NODE_DOF_CALCUL`) + the
  print (`print_materi_stress_force.cc`). See the per-record pages.
- **Lot 2 (commit 1babf2d, 2D numerical integration)** in `calcul_force.cc`:
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
- **Lot 3 (this lot, 3D numerical integration)** in `calcul_force.cc`:
  `msf_calculate_node_3d()` (per-node scan, 16 items) →
  `msf_element_contribution_3d()` (face selection from
  direction_exclude/include, per-face frame, node roles) →
  `msf_integrate_face_3d()` (the 2D stress integration over one face).
  No new enums: the averaged flag exists since lot 2. The face tables
  `msf_border_nodes_hex8/hex27` are static duplicates of the
  `border_nodes_*` of `area.cc` (S Y N C comment in the source; the
  originals are static there). `validate()` now also restricts the 3D
  target groups to hex8/hex27 (the manual's isoparametric section
  elements) and the "not yet implemented" notice is gone.

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
   reciben nada (decisión documentada). Implementado en 2D y 3D.
8. **plot_switch**: -yes invierte las componentes x/y del item
   (dirección de dibujo del vector, manual 6.916); la componente s
   (tamaño) no cambia. Implementado en 2D (3 switches) y 3D (4
   switches, manual 6.916).

## Diseño 3D (lote 3; decisiones con evidencia; el test analítico es el árbitro)

1. **Selección de caras** (manual 6.909/6.911): para cada una de las 6
   caras del hex8/hex27 (tablas `msf_border_nodes_*`, SYNC con
   area.cc) la normal física n = producto vectorial de las 2 aristas
   de esquina, orientada hacia FUERA del centroide (patrón
   interface.cc). Con `direction_exclude`: se descartan las caras con
   |n·dir| > 1−eps (los extremos de la estructura, p.ej. las caras
   axiales del túnel); con `direction_include`: las caras con
   |n·dir| < eps (las superficies de la estructura). El chequeo de
   cordura del manual 6.913 "exactly 4 sides should be consistent
   with the specified direction" = exactamente 4 caras con |n·dir| <
   eps (el elemento de sección tiene 4 caras ⊥ dir y 2 ∥ dir); si no,
   aviso único y el elemento se omite. Las CANDIDATAS (las caras no
   descartadas) son: exclude → las 4 caras ⊥ dir (exterior/interior/
   ±circunferenciales del túnel); include → las 2 caras ∥ dir (las
   secciones ±altura del sheet pile).
2. **Caras extremas** (6.908 "two opposing end faces"): las 2
   candidatas con |n·t̂| MÍNIMA, t̂ = normalizado(centroide −
   reference_point) (decisión 2D-1 generalizada; el centroide = la
   media de las 8 ESQUINAS — para el hex27 las posiciones
   0,2,6,8,18,20,24,26, NO los nodos 0-7 que son la cara z=−1 —
   GOTCHA verificado en msf_tunnel3d). Túnel: entre las 4 candidatas
   el par ⊥ radial = las 2 ±circunferenciales (donde actúa la hoop);
   sheet pile: las únicas 2. Las otras candidatas (superficies
   exterior/interior del túnel) son superficies de la estructura y no
   producen valores primarios (patrón 2D). Guardas: (a) las 2 caras
   extremas deben apuntar a LADOS OPUESTOS del elemento (producto
   escalar de las normales < 0.1 — NOTA: para los elementos curvos
   del anillo las caras ±θ de un elemento de 45° NO son paralelas
   (dot = −cos 22.5° = −0.92), el chequeo |n_a·n_b| ≈ 1 sería
   incorrecto — GOTCHA verificado); (b) la 2ª cara debe estar
   claramente separada de la 3ª candidata (elemento distorsionado →
   aviso único y omisión).
3. **Frame por cara**: n (normal saliente), t = la dirección en la
   cara de MENOR extensión física (manual 6.917; thickness_switch
   -yes → la mayor), orientada hacia FUERA del reference_point (la
   dirección de dibujo); l̂ = n×t̂ (dirección de longitud). l = el
   tamaño del elemento en dirección de longitud = la extensión de la
   cara proyectada sobre l̂ (manual 6.913: "the length of an element
   is determined from the nodal coordinates differences in length
   direction"). Para el túnel l = el tamaño axial (la hoop por unidad
   de longitud axial); para el sheet pile l = el ancho.
4. **Cuadratura de cara (2D)** — consistente con la decisión 2D-3:
   Gauss(2)×Gauss(2) para hex8 (4 puntos; el integrando del momento es
   cuadrático) y Lobatto(3)×Lobatto(3) para hex27 (9 puntos, Simpson;
   el integrando del momento es cúbico). Los puntos de la cara usan
   la interpolación del elemento EN la cara (la tensión σ se
   interpola de los valores NODALES con las funciones de forma de la
   cara — decisión 2D-2). El factor 4 = los dos intervalos [−1,1] de
   la integral 2D (análogo del side_len implícito de la integral 1D
   del lote 2) — GOTCHA verificado: sin el factor los valores salían
   4× pequeños (nor = 0.025 en vez de p·R = 0.1).
5. **Valores** (por unidad de longitud l): nor = ∫∫σ_nn dA / l (con
   signo, tracción +), she = |∫∫σ_nt dA| / l (solo tamaño), mom1 =
   ∫∫σ_nn·dt dA / l y mom2 = ∫∫σ_nn·dl dA / l. dt y dl se miden desde
   el PUNTO MEDIO DE LA CARA (la media de sus nodos), no desde el
   centroide del elemento: para los elementos CURVOS del anillo el
   centroide por esquinas queda en la CUERDA (no en el arco medio) y
   añade un momento espurio (verificado en msf_tunnel3d: mom1 ≈ 0.015
   con el centroide vs 3.6e-12 con el punto medio de la cara); para
   los elementos rectos ambas definiciones coinciden con la decisión
   2D-4.
6. **Asignación por nodo** (3D): los nodos de las 2 caras extremas
   (18 en hex27, 8 en hex8 — caras opuestas disjuntas) son PRIMARIOS;
   los 9 nodos del hex27 que NO están en ninguna cara extrema (el
   plano medio entre las caras, 6.908) reciben el PROMEDIO de las dos
   caras con average -yes (default) y se marcan averaged (-primary los
   omite: hex27_avg 45 vs 27, tunnel 144 vs 72). El vector de plot de
   los nodos promediados se dibuja a lo largo de t̂ del elemento (los
   t̂ de las dos caras curvas difieren).
7. **Hex8**: sin promedio (average solo para quad9/hex27) → -all ==
   -primary. El hex8 1-en-espesor NO puede representar el campo de
   flexión cuadrático (u_x ~ y²) → shear locking en flexión (el mismo
   fenómeno del quad4 2D del lote 2) — documentado en
   msf_sheet3d_hex8 (el test usa CORTE PURO, campo lineal EXACTO
   para el hex8: she = G·γ·t = 0.5 EXACTO). El fix (SRI/hex8
   cuadrático) es un lote aparte.

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
| msf_sheet3d (hex27, flexión prescrita) | mom1s | E·κ/12 = 0.0833333 | 0.0833333 EXACTO (81/81 nodos) |
| msf_sheet3d | nors / shes | 0 | < 1e-8 |
| msf_sheet3d_hex8 (corte puro) | shes | G·γ·t = 0.5 | 0.5 EXACTO (20/20 nodos) |
| msf_tunnel3d (anillo hex27) | nors | E·u0·t/R = p·R = 0.1 | 0.1 EXACTO (144/144 nodos) |
| msf_tunnel3d | mom1s / shes | 0 | 3.6e-12 / < 1e-8 |
| msf_hex27_avg | -all vs -primary | 45 vs 27 | 45 vs 27 (18 promediados) |

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
  quad4/quad9 (error claro para otros tipos). En L3: en 3D solo
  hex8/hex27 (el manual 6.913). El aviso "not yet implemented" se
  elimina en L3.
- **GOTCHA MAYOR del GNU (3D, documentado — fuera del alcance del
  integrador)**: el solve mixto 3D materi_stress + BiCG por defecto es
  DEGENERADO. La matriz v-v del hex8/hex27 con la integración default
  (Lobatto 2×2×2 / 3×3×3) tiene una estructura incompleta: para
  muchas configuraciones de carga A·b ≈ 0 y el BiCG se detiene en la
  iteración 0 con x = 0 (medido: dAd = 3.9e-29 con la carga en la
  cara superior del cantilever) o diverge (hex27 cantilever: error
  final 2.16e+13 tras 720 iteraciones, "did not converge"). El anillo
  con presión convergía mal (meseta ~1e-12, tensiones contaminadas).
  SÍ funcionan: (a) el barrido axial (elasti3d: carga en la cara x=1
  con velx fijo en x=0 → σ_xx = 0.02 EXACTO), (b) los solves
  Dirichlet (velocidades prescritas, iface_3d). CONSECUENCIA: los
  tests 3D del lote 3 usan DEFORMACIÓN PRESCRITA (Dirichlet puro, sin
  solve): campos lineales/cuadráticos que el hex8/hex27 integran
  EXACTO, con expectativas analíticas EXACTAS (flexión pura del sheet
  pile, expansión radial proporcional del túnel). El cantilever con
  carga en punta y el vaso de presión con cargas consistentes quedan
  PENDIENTES del arreglo del solve mixto 3D (lote aparte, fuera del
  sub-sprint MSF). Cada `.dat` documenta el campo prescrito y su
  expectativa.

## Pendiente

- **Arreglo del solve mixto 3D del GNU** (matriz v-v degenerada del
  hex8/hex27 + BiCG): sin él, los tests de carga (cantilever con
  fuerza en punta, vaso de presión con cargas consistentes) no pueden
  validarse — lote aparte, fuera del sub-sprint MSF. Los tests 3D del
  lote 3 usan deformación prescrita (ver GOTCHA MAYOR arriba).
- axisimétrico: implementado (l=2πr) pero sin test dedicado en este
  lote (pendiente; el MSF L4 previsto pulirá el caso).
- Los tests `.dat` viven en `validation-suite/test-2014/` (gitignored;
  solo el bucle y los checks de `scripts/build_safe.sh` se versionan).
