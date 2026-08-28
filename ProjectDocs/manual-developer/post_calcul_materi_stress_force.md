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
- **Lot 3 (commit fa7650b, 3D numerical integration)** in `calcul_force.cc`:
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
- **Lot 4 (commit 631258e, 2026-08-28)**: the section stress SOURCE
  changed from the recovered nodal stresses (`NODE_DOF`, the L2
  decision) to the ELEMENT integration-point stresses (`ELEMENT_DOF`;
  `msf_element_rule` replicates the element's own quadrature incl. the
  SRI Gauss switch and the axisymmetric 1-point rule; the section
  points reconstructed with the Lagrange weights of the IP grid). The
  recovered nodal stresses remained only as the documented FALLBACK
  for the all-zero ELEMENT_DOF blocks (the solved 3D `derivatives`
  models - gforce10/gforce13). `elem.cc`: the ELEMENT_DOF WRITE covers
  the FULL stress block (6 components × nder slots). Axisymmetric: the
  section integrand carries the physical circumference 2π·r.
  `validate()` requires `options_element_dof -yes`. No new enums.
- **Lot 5 (this lot, 2026-08-28)**: the section forces change from the
  stress-FIELD integration (lots 2-4: the quadrature of σ_nn/σ_nt over
  the section faces) to the EQUILIBRIUM-based method: the section
  resultant over an end face = the sum of the ELEMENT INTERNAL FORCES
  f_elem = ∫Bᵀσ dV (the consistent nodal forces of the IP stress
  field) of the face nodes. New helpers in `calcul_force.cc`:
  `msf_element_internal_forces_2d/3d()` (the internal forces: the B
  matrix of polynom.cc:549-599 contracted with the IP stresses, the
  physical derivatives dn = invJ·p and the IP volumes w·4·|detJ| (2D) /
  w·8·|detJ| (3D) with the 2π·r axisymmetric factor - the SAME
  kinematics of pol()/materi(), the "element forces ... setup in a
  timestep" of the manual 6.913) and `msf_face_resultants_2d/3d()`
  (nor = n̂·R_face, she = |t̂·R_face|, mom = Σ (n̂·f)·arm about the
  middle of the face). The SRI quad4 gets the ELEMENT-CONSISTENT
  internal force (the full-rule shear part is replaced by the reduced
  1-point shear internal force of the momentum feedback:
  4·detJ_c·b_shear·mean(IP σ_xy) - the shear modulus cancels because
  σ_xy = G·γ and the 2×2 Gauss mean of the bilinear γ = the centroid
  value exactly; without it the SRI section would read the K_full
  (locked) equilibrium instead of K_SRI·u = P). The face quadrature
  and the Lagrange stress reconstruction of LOT 4 are retired. No new
  enums.

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
    el elemento se omite (valores 0). (Los LOT 2-5 usan la MISMA
    selección de caras; lo que cambia es la FUENTE de los valores.)
 2. **Fuente de σ (LOT 4-5)**: la tensión de los PUNTOS DE INTEGRACIÓN
    del elemento (`element_dof[ipoint*nuknwn + stres_indx + comp*nder]`
    = las tensiones constitutivas que el elemento USÓ en su última
    pasada - "the element forces needed for this option are setup in a
    timestep", manual 6.913). `options_element_dof` debe ser -yes
    (default; error claro en la validación si no). LOT 4: la tensión en
    un punto de cuadratura de la sección se reconstruía con los pesos
    de LAGRANGE de la regla del propio elemento (para las reglas con
    nodos los puntos de la sección COINCIDEN con IPs de la cara —
    lectura directa). LOT 5: la tensión se lee SOLO en los IPs del
    elemento y alimenta la integral de las fuerzas internas (decisión
    9) — la reconstrucción de Lagrange y la cuadratura de cara se
    retiran. La fuente antigua (la tensión nodal recuperada, decisión
    L2) queda SOLO como fallback: medido (2026-08-28), para los
    modelos 3D RESUELTOS con `derivatives` (gforce10/gforce13) el
    bucle escalonado NO propaga la deformación de la velocidad
    convergida a los IPs del elemento (el bloque de tensión del
    ELEMENT_DOF sale TODO CEROS mientras el NODE_DOF recuperado tiene
    los valores correctos) → en ese caso las fuerzas internas se
    calculan con las tensiones nodales evaluadas en las posiciones de
    los IPs con las funciones de forma (warning único documentado).
    Verificado además (L4): para el quad9/quad4 la tensión nodal
    recuperada es la media EXACTA de los IPs de los elementos
    adyacentes (la recuperación Lobatto es exacta), por lo que el
    cambio de fuente NO limpia la polución del arness (la polución
    vive en el propio campo σ).
 3. **Cuadratura (histórica, LOT 2-4)**: `integration_gauss(2)` /
    `integration_lobatto(3)` a lo largo de la arista/cara. LOT 5: la
    cuadratura de la sección se retira; la regla del elemento
    (msf_element_rule) se usa SOLO para la integral de las fuerzas
    internas (decisión 9).
 4. **Momento**: LOT 5: mom = Σ_{nodos de la cara} (n̂·f_elem)·arm con
    arm = (x_nodo − x_medio_de_la_cara)·t̂ — el momento de las fuerzas
    internas consistentes respecto al PUNTO MEDIO DE LA CARA (la
    decisión de brazos del L3; el equivalente discreto EXACTO del
    ∫σ_nn·dt del manual: el momento de las fuerzas nodales consistentes
    sobre el punto medio de la cara = el momento de las tracciones de
    la cara — las componentes de cizalla contribuyen solo a la
    componente en el plano del vector momento, que no se usa). NOTA:
    el momento incluye el momento DÉBIL completo de las fuerzas
    internas sobre el punto de la sección — para un estado de CORTE
    PURO el moms = el momento del par de reacciones del bloque (p.ej.
    0.1923 en msf_shear), mientras que el ∫σ_nn·dt del manual es 0
    (documentado en los tests; el momento del método de equilibrio es
    el consistente con la definición de fuerzas internas).
 5. **Valores** (por unidad de longitud l; 2D plano l=1, axisimétrico
    l=2π·r con r = coordenada radial del centroide, convención de
    area.cc): LOT 5: nor = n̂·R_face/l (CON SIGNO, tracción +),
    she = |t̂·R_face|/l (solo tamaño, manual), mom = Σ (n̂·f)·arm/l (con
    signo), con R_face = Σ de las fuerzas internas de los nodos de la
    cara (decisión 9). Las componentes de plot: norx/nory = nor·t̂,
    nors = |nor|; shex/shey = |she|·t̂, shes = |she|; momx/momy =
    mom·t̂, moms = |mom|. El signo tracción/compresión de nor y mom
    viaja en la DIRECCIÓN del vector de plot (t̂ hacia/desde el
    reference_point); la componente s es el tamaño físico (manual
    6.913).
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
 9. **FUERZAS INTERNAS DEL ELEMENTO (LOT 5)**: f_elem[inod, idim] =
    Σ_ip volumen[ip]·(Bᵀσ)[ip, inod, idim] — las fuerzas nodales
    consistentes del campo de tensiones de los IPs, la MISMA cantidad
    que materi() acumula en element_rhside con el SIGNO OPUESTO
    (element_rhside -= volume·force, materi.cc:618-619: la fuerza
    interna es la que el elemento ejerce SOBRE sus nodos, K·u,
    opuesta a las cargas). La cinemática replica pol()/materi():
    funciones de forma y derivadas en los IPs (interpolation_polynomial
    con la regla de msf_element_rule — la MISMA regla con la que el
    elemento integró), las derivadas físicas dn = invJ·p, el volumen de
    IP w·4·|detJ| (2D) / w·8·|detJ| (3D) con el factor 2π·r
    axisimétrico, y la matriz B de polynom.cc:549-599. FÍSICA: las
    fuerzas internas están en equilibrio con las cargas POR
    CONSTRUCCIÓN para el solve convergido (Σ sobre los elementos =
    −P en los dofs libres; la suma sobre los nodos de UN elemento es la
    identidad de traslación rígida, idénticamente cero). El resultante
    de sección sobre una cara extrema = Σ de las fuerzas internas de
    los nodos de la cara = ∫σ·n̂ dA de la cara (la identidad de las
    fuerzas nodales consistentes) = el resultante del cuerpo libre de
    las cargas aplicadas: ESTÁTICA EXACTA — la misma familia que usa el
    Professional ("the element forces needed for this option are setup
    in a timestep"; verificado: el arness del túnel msf_tunnel3d da los
    valores del node_dof_calcul del Professional DÍGITO A DÍGITO —
    nor 0.0998068, shes 0.013727, mom1 −2.0e-4 — mientras la
    integración de campo del L4 daba el valor analítico 0.1 EXACTO).
    El momento de las fuerzas nodales consistentes sobre el punto medio
    de la cara = el equivalente discreto EXACTO del ∫σ_nn·dt dA (las
    componentes de cizalla contribuyen solo a la componente en el plano
    del vector momento, que no se usa). SRI quad4 (opt-in): la parte de
    cizalla de regla completa de f_elem se reemplaza por la fuerza
    interna de cizalla reducida de 1 punto (el feedback del momento del
    punto fijo K_SRI·u = P): 4·detJ_c·b_shear·media(σ_xy de los IPs) —
    el módulo de cizalla se cancela (σ_xy = G·γ y la media Gauss 2×2 de
    la γ bilineal = el valor del centroide EXACTAMENTE). Sin la
    corrección la sección SRI leería el equilibrio K_full (lockeado)
    en vez de K_SRI·u = P. VERIFICACIÓN DEL EQUILIBRIO: las fuerzas de
    sección medidas = la estática del cuerpo libre de las cargas: la
    ménsula msf_beam2d da moms = P·(8−x) y shes = P EXACTOS (6 dígitos)
    en TODAS las secciones; la biempotrada gffq4 cumple la identidad
    |M_e|+|M_c| = pL²/8 a 0.9996; la ménsula quad4 gforce7q4 da
    N/V EXACTOS y M al 99.84% (las fuerzas internas de la solución
    lockeada EQUILIBRADA también cumplen la estática de las cargas).
    LIMITACIÓN MEDIDA: el resultante débil solo es la estática del
    cuerpo libre cuando el campo σ está EN EQUILIBRIO (el solve
    convergido); para los runs gruesos multi-paso del GNU (quad9/hex8)
    el estado σ recuperado NO está en equilibrio (propiedad preexistente
    del esquema escalonado, documentada en el baseline L4 como la
    polución N 1.24×/V 2.7× del gforce7) y el resultante débil lo
    AMPLIFICA (el gforce7 mide 5×; la integración de campo del L4 era
    menos sensible). El método es el correcto; el estado σ del GNU para
    esos runs es el limitante (el Professional lo produce en equilibrio
    con su solve).

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
| msf_beam2d (8×quad9, L=8, P=1e-2) | moms en x=0..8 | P·(8−x) | 0.079999, 0.060000, 0.040000, 0.020000, 0.000000 — EXACTOS (6 dígitos, LOT 5) |
| msf_beam2d | shes | P | 0.0100 EXACTO en TODAS las secciones (LOT 5; la banda ±30% de la integración de campo desaparece) |
| msf_beam2d | nors | 0 | ~5e-7 EXACTO |
| msf_beam2d_pure (4 puntos) | moms tramo central | 1e-2 cte | 0.00999994 EXACTO (LOT 5; antes 1.5%) |
| msf_beam2d_pure | shes tramo central | 0 | 1e-6 (≈0 EXACTO) |
| msf_quad9 | -all vs -primary | 27 vs 15 | 27 vs 15 |
| msf_nor (8×quad4 axial, N=2e-2) | nors | 2e-2 | 0.0200 EXACTO (LOT 5) |
| msf_shear (corte simple) | shes | G·γ·h = 0.384615 | 0.384615 EXACTO |
| msf_shear | moms | 0 (manual) / par de reacciones | 0.1923 = 0.3846·0.5 (LOT 5: el momento débil del par de reacciones del bloque) |
| msf_sheet3d (hex27, flexión prescrita) | mom1s | E·κ/12 = 0.0833333 | 0.0833333 EXACTO (81/81 nodos) |
| msf_sheet3d | nors / shes | 0 | < 1e-8 |
| msf_sheet3d_hex8 (corte puro) | shes | G·γ·t = 0.5 | 0.5 EXACTO (20/20 nodos) |
| msf_sheet3d_hex8 | mom1s | 0 (manual) / par de reacciones | 0.25 = 0.5·0.5 (LOT 5) |
| msf_tunnel3d (anillo hex27) | nors | E·u0·t/R = p·R = 0.1 | 0.0998068 — el valor FE DISCRETIZADO, IDÉNTICO al Professional (0.09980686) (LOT 5) |
| msf_tunnel3d | shes / mom1s | 0 | 0.013727 / 2.03e-4 — IDÉNTICOS al Professional (0.013727 / −2.0e-4) (LOT 5) |
| msf_hex27_avg | -all vs -primary | 45 vs 27 | 45 vs 27 (18 promediados) |
| msf_cant3d_hex27 (L5, cantilever 3D hex27 CARGADO) | mom1s en z=0..4 | P·(L−z) = 0.04..0 | 0.0403, 0.0297, 0.0199, 0.00997, 0.00016 (dentro del 1%; antes 9%) |
| msf_cant3d_hex27 | shes | P | 0.0097..0.0103 = P dentro del 4% (la cizalla 0.08·P de los IPs desaparece — la resultante del equilibrio) |
| msf_axisym (L5, anillo axisimétrico prescrito) | nors | σ_zz·t = 0.5 EXACTO por unidad de circunferencia | 0.5 EXACTO (4/4 nodos) |
| msf_axisym | shes / moms | 0 | 0 / 0 (LOT 5; antes mom 0.001) |
| arness gforce7 (L5, quad9 2-el, 30 pasos) | N/V/M en x=50 | −12.34/+100/−5000 | 61.7/500/24955 (5×) — el estado σ del run NO está en equilibrio (preexistente: el baseline L4 media N 1.24×/V 2.7× por integración de campo; el resultante débil lo amplifica) |
| arness gforce7q4 (L5, quad4 2-el) | N/V/M en x=50 | −12.34/+100/−5000 | 12.34 EXACTO / 100 EXACTO / 4992 (99.84%) — la solución lockeada EQUILIBRADA cumple la estática de las cargas |
| arness gffq4 (L5, biempotrada 10 quad4) | M_end/M_center | −825/+425 (deep) | −824.99 (1.0000×) / +424.46 (0.9987×); identidad \|M_e\|+\|M_c\| = 1249.45 = pL²/8 a 0.9996 — la estática del cuerpo libre de la solución equilibrada |
| arness gforce10/13 (L5, hex8 3D) | nors en z=50 | −12.34 | 12.34 EXACTO (1.0000×) — vía el fallback a las tensiones nodales; V/mom de la sección NO son la estática del cuerpo libre (el estado σ recuperado no está en equilibrio — el fallback) |
| qsri_beam2d / _sri (L5) | mom en x=0 / she en x=4 | P·8 / P | 0.08 EXACTO / 0.01 EXACTO con y sin SRI — la estática del cuerpo libre de las soluciones equilibradas (el lock queda en la deflexión y el campo σ, no en las fuerzas de sección) |

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
- **GOTCHA (L4, promediado quad9/hex27)**: con la fuente IP las caras
  de los elementos adyacentes en una sección compartida DIEREN (el
  campo σ es discontinuo entre elementos; las tensiones nodales
  recuperadas eran su media exacta). Los nodos PRIMARIOS de la sección
  muestran la media de las dos caras (linealidad: idéntico a la fuente
  antigua), pero los nodos PROME-DIADOS del plano medio usan las caras
  del PROPIO elemento → ya no coinciden exactamente con la media de
  los valores de salida (banda FE medida: desviación < 2e-3 absoluta;
  antes era exacta). El comportamiento es el CORRECTO del manual 6.906
  ("the averaged values from the forces and moments on the two
  opposing end faces" DEL ELEMENTO).
- **GOTCHA (L4, orden del promedio de la cizalla)**: la cizalla se
  promedia como |∫σ_nt| POR ELEMENTO y luego se promedian las
  magnitudes (la integral firmada se pierde): los checks de
  msf_beam2d_pure/qsri se actualizaron a los nuevos valores medidos
  (la cancelación de signos del campo nodal promediado no se
  conserva).
- **GOTCHA (L4, SRI quad4)**: la cizalla de la sección lee ahora el
  campo σ_xy CRUDO de los puntos de Gauss (la recuperación nodal
  h-weighted del fix D-b diluía la cizalla; el DIAG 12.4 documenta que
  el σ_xy del Q4 NO es superconvergente) → el shes del SRI sube a
  8.4·P en la sección interior (medido en qsri_beam2d_sri); el
  MOMENTO SRI no cambia (0.9375·P·L, la referencia de Hughes — el
  campo σ_nn reconstruido de los IPs Gauss = el mismo que la
  extrapolación nodal del D-b).
- El flag averaged se pre-aloca en `calculate()` (db_allocate) como
  NODE_DOF_CALCUL: el PUT dentro del bucle paralelo no puede alocar.
  version_all=1 → db_version_copy + renumbering lo llevan a
  VERSION_PRINT para el print (como NODE_DOF_CALCUL).
- `validate()` añade en L2: `materi_stress` obligatorio en initia
  (stres_indx ≥ 0) y, en 2D, los grupos objetivo solo pueden contener
  quad4/quad9 (error claro para otros tipos). En L3: en 3D solo
  hex8/hex27 (el manual 6.913). En L4: `options_element_dof` debe ser
  -yes (la fuente IP del ELEMENT_DOF; error claro si no). El aviso
  "not yet implemented" se elimina en L3.
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
  expectativa. **Diagnóstico completo (2026-08-28) en
  [DIAG-SOLVE-MIXTO.md](../DIAG-SOLVE-MIXTO.md)**: el sistema que Bi-CG
  resuelve es SOLO de velocidades y SPD (los dofs σ nunca entran en la
  matriz global; σ avanza por una actualización diagonal lumped —
  esquema escalonado, no saddle-point); los criterios de parada del
  Bi-CG declaran éxito en breakdown/estancamiento con error ≫
  tolerancia, y el punto fijo del esquema converge al estado locked
  (0.2315× para quad4 plain Y SRI; SuperLU ≡ Bi-CG byte-idénticos).
  **ARREGLADO (2026-08-28, lote A+B del diagnóstico — ver
  [solve_iterative_bicg.md](solve_iterative_bicg.md))**: criterios de
  parada honestos + CG para el sistema SPD. El cantilever 3D con
  carga REAL ya resuelve (hex8 1-e: converge a 4.2e-37; hex8×8/hex27×8:
  fuerzas de sección — hex27 mom ≈ 1.1·P·(L−x)) y la validación 3D con
  carga queda desbloqueada (el test msf_cant3d_hex27 del L4 la
  verifica: mom1 = P·(L−z) dentro del 9%; la cizalla de sección del
  esquema mixto 3D queda contaminada — 0.08·P medido — la "she ≈ P"
  del diagnóstico usaba una dirección de carga distinta y el campo
  nodal, ver el gotcha L4). **Fix C/D (2026-08-28, DIAG §12)**: el
  punto fijo del esquema = la solución del elemento (SRI quad4
  0.9375×; plain byte-idéntico; quad9/hex8 a su valor) — la estática
  exacta del Professional requiere además la integración por fuerzas
  internas (ver Pendiente).

## Pendiente

- **GOTCHA (malla 1-en-espesor, pre-L4 — la banda de cizalla)**:
  la cizalla integrada de la formulación mixta tiene polución en las
  superficies libres (condición de superficie libre débil): el
  quad9 1-en-espesor da |∫σ_xy dy| en la banda ±30% de P según la
  sección (medido en msf_beam2d: 0.0109..0.0133 con P=1e-2; la banda
  NO se aprieta con la fuente IP del L4 — la polución vive en el campo
  σ, los IPs = la misma polución que la nodal recuperada; el valor
  EXACTO de she se verifica con msf_shear, campo uniforme). El quad4
  1-en-espesor sufre shear locking (mom ≈ 0.23× y she ≈ 0.74× del
  valor estático) → los tests de flexión usan quad9; el SRI quad4
  (opt-in, familia qsri) mejora el momento a 0.9375× (fix C/D) pero la
  cizalla leída de los IPs Gauss queda contaminada (8.4·P medido — el
  σ_xy del Q4 no es superconvergente, DIAG 12.4). El hex8 3D sufre el
  mismo fenómeno (u_x ~ y² no representable); el SRI hex8 es trabajo
  futuro.
- **Fix C/D del esquema escalonado**: DONE (2026-08-28, ver
  [DIAG-SOLVE-MIXTO.md](../DIAG-SOLVE-MIXTO.md) §12): el punto fijo del
  esquema = la solución del elemento (SRI quad4 0.9375×, plain
  byte-idéntico, quad9/hex8 a su valor). La estática de sección del
  Professional es EXACTA en TODOS los casos (1e-10) porque NO integra
  el campo σ crudo (nodal NI de IP — la polución vive en el campo): su
  `node_dof_calcul` es consistente con un cálculo por FUERZAS
  INTERNAS/equilibrio del elemento ("the element forces needed for
  this option are setup in a timestep", manual 6.913). **IMPLEMENTADO
  en el LOT 5 (2026-08-28)**: las fuerzas de sección se calculan desde
  las fuerzas internas del elemento (f_elem = ∫Bᵀσ dV) — la estática
  del cuerpo libre de las cargas, EXACTA para los estados σ en
  equilibrio (verificado: msf_beam2d P·(8−x)/P a 6 dígitos, gforce7q4
  N/V exactos, la identidad gffq4 pL²/8 a 0.9996, y el túnel
  msf_tunnel3d = los valores del node_dof_calcul del Professional
  DÍGITO A DÍGITO). LIMITACIÓN MEDIDA (frente abierto): el resultante
  débil de las fuerzas internas solo es la estática del cuerpo libre
  cuando el campo σ está EN EQUILIBRIO con las cargas (el solve
  convergido); los runs gruesos multi-paso del GNU con quad9/hex8
  (gforce7, gforce7_ref, gforce10/13) tienen el estado σ recuperado
  NO en equilibrio — propiedad PREEXISTENTE del esquema escalonado
  (el baseline L4 la documentó como la polución N 1.24×/V 2.7× por
  integración de campo; el resultante débil la amplifica a 5×/1.25×).
  El Professional produce el estado σ en equilibrio (su solve
  monolítico/iterado) y por eso su estática es exacta en las mismas
  mallas; la vía de cierre para esos runs es del SOLVER (hacer que el
  σ satisfaga Bᵀσ = P), no del post-proceso de sección.
- axisimétrico: DONE (L4-5, msf_axisym): la convención l=2πr (fuerzas
  por unidad de circunferencia), verificado EXACTO (nor = 0.5, she =
  0, mom = 0 con las fuerzas internas).
- Los tests `.dat` viven en `validation-suite/test-2014/` (gitignored;
  solo el bucle y los checks de `scripts/build_safe.sh` se versionan).
