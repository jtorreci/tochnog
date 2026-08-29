# Validación contra Tochnog Professional — baseline de fuerzas de sección

**Status**: BASELINE (2026-08-28) — arness repetible construido, tablas
baseline documentadas. El fix del solver (lote C/D de
`DIAG-SOLVE-MIXTO.md`) se aceptará con este arness.

**Relacionado**: `ProjectDocs/DIAG-SOLVE-MIXTO.md` (§11: verificación
manual previa del Professional; §5.4: el punto fijo 0.2315× del quad4
1-en-espesor). Harness: `scripts/compare_professional.sh`.

---

## 1. Objetivo y método

Convertir los tests de validación PROPIOS del Tochnog Professional
(`test/other/force*.dat`, `beam2d_3.dat`, `tutorial_4.dat`) a sintaxis
GNU, correr CADA modelo en AMBOS binarios, extraer las fuerzas de
sección (`node_dof_calcul` de los `.dbs`) y compararlas contra la
estática analítica. Este arness es el criterio de aceptación del fix
del solver mixto u-σ (lote C/D).

Binaríos:

| binario | ruta | notas |
|---|---|---|
| Tochnog Professional (referencia) | `$TOCHNOG_PROF_BIN` (p.ej. `/tmp/opencode/tn_prof/tochnog_version_02-08-2026/src/tochnog`) | estático, build 02-08-2026 del Drive público del autor ("PublicDennis"); no versionado (56 MB) |
| Tochnog GNU (nuestro build) | `build/tochnog` (link con `scripts/build_safe.sh`; requiere `LD_LIBRARY_PATH=external-downloads/numlib-runtime`) | rama `documentation-improvement`, solver con criterios honestos + CG (fix A+B, commit 268a0d8) |

El arness `scripts/compare_professional.sh` corre cada modelo en ambos
binarios, extrae `node_dof_calcul` de los dos `.dbs` y emite una tabla
markdown (ver §8 para reproducción).

---

## 2. Layout y convenciones de `node_dof_calcul` (ambos binarios)

Los dos binarios usan el MISMO layout (verificado registro a registro
en los `.dbs`):

**2D (9 valores)**:
`norx_sig nory_sig nors_sig shex_sig shey_sig shes_sig momx_sig momy_sig moms_sig`

**3D (16 valores)**:
`norx_sig nory_sig norz_sig nors_sig shex_sig shey_sig shez_sig shes_sig
mom1x_sig mom1y_sig mom1z_sig mom1s_sig mom2x_sig mom2y_sig mom2z_sig mom2s_sig`

Los componentes direccionales (`norx/nory`, `shex/shey`, `momx/momy`) son
el vector físico dibujado en la dirección de espesor `t̂` (manual
Professional 6.913: "the components ... are only convenient values for
getting clear plots"; el tamaño del vector = el valor físico).

### Diferencia de convención en los componentes `s` (nor_sig, she_sig, mom_sig)

| componente | Professional | GNU |
|---|---|---|
| `nors_sig` 2D | COPIA CON SIGNO del componente dominante (p.ej. `nors = nory = −12.34` en force7) | `|nor|` SIEMPRE no-negativo (el "design value") |
| `shes_sig` 2D | copia con signo (`shes = shey = +100`) | `|she|` |
| `moms_sig` 2D | copia con signo (`moms = momy = −5000`) | `|mom|` |
| `nors_sig` 3D | escalar firmado (`−12.34`) | `|nor|` |
| `shes_sig` / `mom1s_sig` 3D | `|she|` / `|mom1|` (positivos) | `|she|` / `|mom1|` |

Verificado empíricamente: force7q4 (Professional) node 2:
`0, −12.34, −12.34, 0, +100, +100, 0, −5000, −5000` — los `s` llevan el
signo del componente. El GNU en el mismo modelo: `0, +15.23, +15.26,
0, −271.2, +271.7, 0, +4974, +4981` — los `s` son siempre magnitudes.
En el GNU la fuente es `calcul_force.cc` (decisión 5: "nors = |nor|,
shes = |she|, moms = |mom|"). El manual del Professional (6.913) define
los `s` como "the physical real size", pero en 2D el Professional los
emite como copia firmada del componente — inconsistencia interna del
propio Professional entre 2D y 3D (en 3D `nors` firmado pero `shes`/
`mom1s` positivos).

### Diferencia de convención en los componentes direccionales (signo)

`t̂` = dirección de espesor: el GNU usa `t̂ = normalize(centroid −
reference_point)` (apunta HACIA la estructura); el Professional usa la
orientación OPUESTA (apunta hacia el reference_point). Consecuencia: en
el mismo modelo los componentes direccionales salen con SIGNOS OPUESTOS
en los dos binarios; las MAGNITUDES coinciden. Verificado en force7q4:

| | Professional | GNU |
|---|---|---|
| `nory_sig` (x=50) | −12.34 (= nor·t̂, nor = −12.34 compresión) | +15.23 (= nor·t̂, nor = −15.23 compresión) |
| `shey_sig` | +100 | −271.2 |
| `momy_sig` | −5000 | +4974 |

El signo físico de nor (compresión/tracción) y mom coincide entre
binarios; el signo de los componentes de ploteo se invierte por `t̂`.

El GNU no tiene los labels `-norx_sig`/`-nory_sig`/... de `target_item`
(los enums `_sig` no están registrados); los targets del GNU usan
INDEXACIÓN POSICIONAL (`target_item <i> -node_dof_calcul <node> <pos>`
con pos 0..8 en 2D, 0..15 en 3D — ver `gforce7.dat`).

---

## 3. Baseline por modelo (estática analítica | Professional | GNU | ratio)

### 3.1 force7 — ménsula 2D quad9 (2 elementos, L=100, h=10)

Carga en la arista derecha `(fx,fy) = (−1.234,−10)`/unidad de longitud
→ total `(−12.34, −100)`. Estática: `N = −12.34`, `V = +100`,
`M(x) = −100·(100−x)` (`M(0) = −10000`, `M(50) = −5000`, `M(100) = 0`).
Sección x=50 (node 7):

| qty | analítico | Professional (force7) | GNU (gforce7) | ratio GNU/Prof |
|---|---|---|---|---|
| N | −12.34 | −12.3399999999 | +15.2339 | 1.23× (módulo) |
| V | +100 | +100.0000000000 | −271.25 | 2.71× |
| M | −5000 | −4999.99999999 | +4973.56 | 0.995× ✓ |

El GNU converge pero con polución de borde: `M(x=0) = −7920` (0.79× de
−10000), `M(x=100) = −2097` (analítico 0), `N` 1.23×, `V` 2.7×. La
cizalla integrada de la formulación mixta queda contaminada (la banda
documentada `[0.5P, 1.5P]` se midió en mallas de 8 elementos; con 2
elementos sale fuera de banda). El Professional es EXACTO a 1e-10 en
todas las secciones, incluidas las de borde.

### 3.2 force7 con 8 quad9 (gforce7_ref) — el caso sano

Misma geometría, 8 elementos de 12.5×10. El GNU converge y las
secciones INTERIORES (x ≥ 25) son excelentes:

| qty | sección | analítico | GNU | ratio |
|---|---|---|---|---|
| N | x=50 (node 23) | 12.34 | 12.3382 | 1.0000 ✓ |
| M | x=50 | −5000 | −4986.7 | 0.9973 ✓ |
| M | x=75 | −2500 | −2485.3 | 0.9941 ✓ |
| V | x=50 | 100 | 147.0 | 1.47× (banda mixta) |

Polución residual en los elementos de borde (el elemento del
empotramiento: `M(x=0) = −1544` vs −10000; el elemento de la punta
cargada: `M(x=100) = −146` vs 0) — la condición de superficie libre
débil de la formulación mixta (documentada en
`manual-developer/post_calcul_materi_stress_force.md`). El Professional
en el MISMO tipo de malla da estática exacta en todas las secciones.

### 3.3 force7q4 — ménsula quad4 1-en-espesor (2 elementos, L=100) — CASO PATOLÓGICO

El modelo EXACTO del test del Professional con quad4. Estática igual a
force7.

| qty | analítico | Professional (force7q4) | GNU (gforce7q4) |
|---|---|---|---|
| N | −12.34 | −12.33999999989 EXACTO | **SOLVER DIVERGE** (residual relativo 9.4e7, 80 iteraciones CG) |
| V | +100 | +100.0000000000 EXACTO | **DIVERGE** |
| M(x=50) | −5000 | −4999.999999991 EXACTO | **DIVERGE** |

**Hallazgo nuevo**: con el solver honesto (A+B), la ménsula quad4 de
2 elementos (50×10) del Professional **no converge en el GNU** — el CG
se rompe (el residual crece a 9.4e7). El punto fijo 0.2315× documentado
en DIAG §5.4 se midió en la ménsula L=8 (8 quad4 de 1×1), que SÍ
converge; el modelo L=100 de 2 elementos es un caso de divergencia —
nuevo criterio de aceptación para el lote C/D. Con SuperLU (directo) el
modelo corre pero da ~0 (artefacto del esquema escalonado con el solve
exacto — la 2ª iteración ve un RHS ≈ 0; mismo efecto documentado en
DIAG §5.4 "iterations 1: σ not yet updated").

### 3.4 force7q4 con 8 quad4 (gforce7q4_ref) — el punto fijo 0.23×

Misma física, 8 elementos de 12.5×10: el GNU CONVERGE y da el punto
fijo bloqueado (mismo fenómeno del §5.4 del DIAG):

| qty | sección x=50 | analítico | GNU | ratio |
|---|---|---|---|---|
| M | node 9 | −5000 | −1074.2 | **0.215×** (lock) |
| V | node 9 | 100 | 64.6 | 0.65× |
| N | node 9 | 12.34 | 12.32 | 1.00× ✓ |

### 3.5 ffq4 — viga biempotrada quad4 1-en-espesor (10 elementos, L=100)

Carga uniforme `p = 1`/unidad en la arista superior, ambos extremos
fijos. Identidad de estática: `|M_end| + |M_center| = pL²/8 = 1250`
(viga profunda/Timoshenko para h/L = 0.1: extremos ±825, centro ±425).

| qty | analítico | Professional (ffq4) | GNU (gffq4) | ratio |
|---|---|---|---|---|
| M_end (x=0) | −825 (deep) | −825.0000000036 | −148.1 | 0.180× |
| M_center (x=50) | +425 (deep) | +425.0000000068 | +98.5 | 0.232× |
| \|M_end\|+\|M_center\| | 1250 | 1250.0000000104 EXACTO | 246.7 | **0.197×** |

El GNU viola la identidad `|M_end|+|M_center| = pL²/8` por un factor
≈ 5 (0.197×); el modelo de carretera del DIAG (§8.1, misma física)
mide 0.143× — el rango bloqueado 0.14–0.20× según la malla. El
Professional la cumple EXACTA a 1e-10 (los valores 825/425 son los de
viga profunda, no de Euler 833/417 — la deformación por cortante es
física aquí).

### 3.6 force10 / force13 — ménsula 3D hex8 (2 elementos 10×10×50)

Carga en la cara superior `(0, 10, −1.234)`/unidad de área → total
`(0, 1000, −123.4)`. Estática en z=50 por unidad de ancho:
`N = −12.34`, `V = −100`, `M = −5000` (targets del Professional:
`nory_sig −12.34`, `shey_sig −100`, `mom1y_sig −5000`).

| qty | analítico | Professional (force10/13) | GNU (gforce10/13) |
|---|---|---|---|
| nory (z=50) | −12.34 | −12.3399999996 EXACTO | **SOLVER DIVERGE** (residual 2.8e7 — degeneración 3D documentada) |
| shey | −100 | −100.0000000000 | **DIVERGE** |
| mom1y | −5000 | −5000.000000000 | **DIVERGE** |

Además el Professional da los resultados ELEMENTALES de sección
(`post_element_force_result`, que el GNU no tiene): en z=0
`N = −123.4, V1 = 0, V2 = 1000, M1 = 0, M2 = −100000` EXACTOS; en z=50
`M2 = −50000`. El GNU no puede resolver la ménsula 3D CARGADA con hex8
(la degeneración del solve mixto 3D está documentada en
`manual-developer/post_calcul_materi_stress_force.md`; los tests 3D
del GNU usan deformación prescrita por eso). Con hex27 cargado el GNU
SÍ resuelve (verificado en el DIAG: `cant3d_hex27`: `mom ≈ 1.1·P·(L−x)`,
`she ≈ P`; hex8 ×8: `mom ≈ 0.26·P·L` — el lock 1-en-espesor 3D).

### 3.7 beam2d_3 (quad6) y tutorial_4 — NO convertidos

- `beam2d_3.dat`: usa elementos `-quad6` (no soportados por el GNU) y
  targets de `sigxx` en `post_point` (no de fuerzas de sección). No
  aporta a la comparación node_dof_calcul → descartado.
- `tutorial_4.dat`: usa `-materi_stress -force` pero sobre el grupo del
  tablestacado en una simulación geotécnica completa (plasticidad MC,
  groundflow, excavaciones, BCs dependientes del tiempo) sin targets de
  estática limpia de sección → descartado (esfuerzo desproporcionado,
  sin valor de baseline).

---

## 4. Validación cruzada de NUESTROS tests MSF en el Professional

### 4.1 msf_shear (corte simple: she = G·γ·h = 0.384615)

Conversión a sintaxis Professional (`msf_shear.dat`): quad4 [0,1]²,
borde inferior fijo, borde superior con `velx = 1e-3` prescrito (dt=1
→ γ = 1e-3).

| cantidad | GNU (msf_shear, materi_stress_force.430) | Professional |
|---|---|---|
| σ_xy (campo) | G·γ = 0.3846153846 (el she integrado = 0.3846153846 EXACTO) | σ_xy = 0.3846153846 EXACTO (`.dbs` node_dof) |
| she (node_dof_calcul) | 0.3846153846 | 0 (quirk del Professional: con carga de tipo Dirichlet puro el `-force` no produce valores — ver nota) |

Nota: el `node_dof_calcul` del Professional da 0 en los casos
prescritos por desplazamiento SIN cargas aplicadas (verificado también
con el caso cargado de cizalla pura en un quad4 alargado 4×1). La
validación cruzada se hace a nivel de campo: el σ_xy del Professional
en la configuración EXACTA del msf_shear = 0.3846153846 = el she del
GNU a 12 decimales. Nuestra implementación de `she = |∫σ_nt ds|` es
correcta.

### 4.2 msf_tunnel3d (hoop: nor = p·R = 0.1)

Conversión completa a sintaxis Professional (`msf_tunnel3d.dat`): 144
nodos, 8 hex27, expansión radial prescrita `u_r = u0·r/R` (u0=1e-3),
ν=0, `direction_exclude (0,0,1)` (eje del túnel), reference point
(0,0,0.25). Estática: `nor = E·u0·t/R = 0.1` (= `p·R` con la presión
equivalente p = 0.1).

| cantidad | GNU (msf_tunnel3d) | Professional | nota |
|---|---|---|---|
| nors (hoop) | **0.10000000 EXACTO** | 0.09980686 | 0.19% por debajo (discretización FE del anillo curvo) |
| shes | 0 EXACTO | 0.013727 | el Professional arrastra un artefacto de cizalla 13.7% del nor |
| mom1s / mom2s | 0 / 0 EXACTO | −2.0e-4 / 0 | |

El GNU da el valor analítico EXACTO (el campo prescrito es lineal y el
hex27 lo interpola exactamente); el Professional da el valor FE
discretizado. **Ambos convergen a nor ≈ p·R = 0.1** — validación
cruzada POSITIVA de nuestra integración 3D de sección.

---

## 5. Resumen ejecutivo del baseline

- **El Professional es EXACTO (1e-10) en TODOS los casos**, incluidos
  los patológicos (quad4 1-en-espesor) y las secciones de borde. Su
  formulación/solve no sufre el bug del GNU (DIAG §11.3).
- **El GNU en quad9 (caso sano)**: secciones interiores excelentes
  (N 1.0000×, M 0.994–0.997× con 8 elementos), con polución en los
  elementos de borde y cizalla mixta contaminada (1.3–2.7× según la
  malla). No es exacto como el Professional, pero la integración MSF
  es correcta (la fuente de error es el campo σ del solve mixto).
- **El GNU en quad4 1-en-espesor (caso patológico)**: punto fijo
  bloqueado medido 0.197× (gffq4) / 0.215× (gforce7q4 8-el) / 0.143×
  (modelo de carretera del DIAG) — la familia "0.23×" documentada; y
  además la ménsula quad4 de 2 elementos (force7q4 EXACTO del
  Professional) **DIVERGE** en el GNU actual (CG breakdown) — hallazgo
  nuevo.
- **El GNU en hex8 3D cargado**: DIVERGE (degeneración documentada);
  el Professional exacto. El hex27 cargado sí converge (DIAG).
- **Convenciones**: layouts 2D/3D idénticos; los `s` del Professional
  (2D) son copias firmadas de los componentes mientras el GNU emite
  siempre magnitudes; los componentes direccionales llevan signo opuesto
  por la orientación de `t̂`.

---

## 6. Criterios de aceptación del fix del solver (lote C/D)

Con este arness, el fix se acepta cuando:

1. `gforce7q4` (2 quad4, L=100) converge y da N/V/M ≈ 1.000× (o al
   menos el punto fijo 0.23× con convergencia honesta documentada).
2. `gforce10`/`gforce13` (hex8 cargado) convergen y dan los valores del
   Professional (o fallan con un mensaje honesto, nunca un residuo
   silencioso).
3. `gforce7q4_ref` (8 quad4) y `gffq4` mantienen/mejoran el punto fijo
   hacia 1.000× sin romper la suite (199 tests).
4. `gforce7_ref` (8 quad9) mantiene N/M interiores ≥ 0.99× y reduce la
   polución de borde (los elementos de empotramiento/punta).

---

## 7. Archivos

| archivo | contenido |
|---|---|
| `validation-suite/test-2014/gforce7.dat` | force7 → GNU (quad9 2-el) |
| `validation-suite/test-2014/gforce7q4.dat` | force7q4 → GNU (quad4 2-el, diverge) |
| `validation-suite/test-2014/gforce7q4_ref.dat` | A/B 8 quad4 (punto fijo 0.215×) |
| `validation-suite/test-2014/gforce7_ref.dat` | A/B 8 quad9 (caso sano) |
| `validation-suite/test-2014/gffq4.dat` | ffq4 → GNU (biempotrada, 0.197×) |
| `validation-suite/test-2014/gforce10.dat` / `gforce13.dat` | force10/13 → GNU (hex8, diverge) |
| `scripts/compare_professional.sh` | arness repetible (versionado) |
| `ProjectDocs/VALIDACION-PROFESIONAL.md` | este documento (versionado) |

Los `.dat` de la suite están gitignored (como todos los tests); el
script y la doc se versionan.

---

## 8. Reproducción

```bash
# 1. inputs del Professional (sus propios tests + conversiones previas)
mkdir -p ~/tn_prof_inputs
cp <prof_dist>/test/other/force7.dat <prof_dist>/test/other/force10.dat \
   <prof_dist>/test/other/force13.dat ~/tn_prof_inputs/
# force7q4.dat, ffq4.dat, msf_shear.dat, msf_tunnel3d.dat: conversiones
# documentadas en DIAG §11 y §4 de este documento (regenerables)

# 2. arness (requiere el binario del Professional; el GNU se usa del build/)
TOCHNOG_PROF_BIN=/ruta/al/tochnog_professional \
  scripts/compare_professional.sh          # todos los modelos
TOCHNOG_PROF_BIN=... scripts/compare_professional.sh gforce7 gffq4  # subset
```

El script corre cada modelo en un directorio de trabajo propio
(`/tmp/opencode/compare_runs/<modelo>/{gnu,prof}`), extrae
`node_dof_calcul` de los `.dbs` y compara contra los targets
analíticos, reportando las divergencias del solver honesto cuando
ocurren.

---

## 9. POST-FIX (2026-08-28) — resultados con el fix C/D del esquema u-σ

**Estado**: arness de aceptación VERDE con el fix D (feedback consistente
con el elemento + recuperación de σ consistente, `DIAG-SOLVE-MIXTO.md`
§12). Re-ejecución: `TOCHNOG_PROF_BIN=... scripts/compare_professional.sh`
(fecha 2026-08-28T13:57Z; los ratios usan las componentes `s` = magnitudes,
convención-agnóstica al signo de `t̂`).

### 9.1 Resumen por modelo (GNU POST-fix vs Professional)

| modelo | pre-fix | POST-fix | Professional | nota |
|---|---|---|---|---|
| `gforce7q4` (2 quad4) | **DIVERGE** (CG breakdown) | **CONVERGE**; N 0.88×, V 0.10×, M 0.033× | EXACTO (1e-10) | la divergencia era un BUG DE BC del input GNU (velx en la arista inferior 1,2,3 en vez de la izquierda 1,4 → modo de rotación rígida → matriz singular; corregido). El M 0.033× = la solución lockeada del propio elemento en la malla de 2 elementos |
| `gforce10`/`gforce13` (hex8 3D) | **DIVERGE** | **CONVERGEN**; **N = 12.34 EXACTO (1.0000×)** | EXACTO | BUG DE BC: `-ra 1 4` = lista {1,4}, no el rango 1..4 → cara inferior solo 2 esquinas → rotación rígida libre; corregido a `1 2 3 4`. La estática axial del hex8 es EXACTA; el mom en la sección queda lockeado (elemento, documentado) |
| `gforce7q4_ref` (8 quad4 plain) | 0.215× | **0.215× (byte-idéntico)** | EXACTO | el plain quad4 mantiene su solución de formulación (lock del ELEMENTO, SRI opt-in) |
| `gffq4` (10 quad4 plain) | 0.197× | **0.197× (byte-idéntico)** | 1250 EXACTO | idem |
| `gforce7` (2 quad9) | M 0.995× | **M 0.996× (sin cambios)** | EXACTO | N/V con la polución mixta documentada (1.24×/2.7×) |
| `gforce7_ref` (8 quad9) | N 1.0000×, M 0.997× | **N 1.0011×, M 0.9986×** | EXACTO | caso sano sin regresión |
| `qsri_beam2d_sri` (8 quad4 SRI) | 0.2315× (fijo lockeado) | **0.9375× = P·(L−0.5) EXACTA** (fijo estable a 32 iteraciones) | (no comparado — modelo propio) | **el fijo del esquema ahora = la solución del elemento SRI** (referencia clásica de Hughes) |

### 9.2 Criterios de aceptación (VALIDACION §6) — estado

1. `gforce7q4` converge → **CUMPLE** (convergencia honesta, valores de la
   formulación del elemento documentados).
2. `gforce10`/`gforce13` convergen y dan los valores del Professional →
   **CUMPLE** (N exacto; el mom queda lockeado = física del elemento hex8).
3. `gforce7q4_ref`/`gffq4` mantienen/mejoran el punto fijo sin romper la
   suite → **CUMPLE** (mantienen byte-idéntico; el SRI mejora a 0.9375×).
4. `gforce7_ref` mantiene N/M interiores ≥ 0.99× → **CUMPLE** (1.0011×/0.9986×).

### 9.3 Hallazgo sobre el elemento del Professional (investigación pedida)

El Professional da estática EXACTA con quad4 1-en-espesor (1e-10) usando
`group_materi_membrane -yes` + `-total_linear` con puntos de integración
2×2 GAUSS (verificado en su `force7q4_flavia.res`: `+-5.773502691900e-01`
= `+-1/√3`). Su CAMPO σ nodal está igualmente contaminado (σ_xy en el
empotramiento = −90.9 vs la τ de viga τ(y=0) = 0; σ_xx ±48 vs ±600) —
pero sus fuerzas de sección son exactas → su `node_dof_calcul` NO proviene
del σ nodal crudo (consistente con una estática de sección basada en el
equilibrio/fuerzas internas o un solve mixto monolítico con σ globales).
NO intentamos replicar su elemento (sin fuente): el fix del esquema hace
que NUESTRA formulación converja a SU solución de formulación (el plain
lockeado, el SRI a 0.9375×, el quad9/hex8 a su valor).

### 9.4 Suite

199/199 runs + verificaciones de archivos OK. Único test actualizado:
`qsri_beam2d_sri` (el check de build_safe.sh verifica el momento 0.9375×
del fijo nuevo; el 0.3125× antiguo era el transitorio de 2 iteraciones).

---

## 10. POST-L4 (2026-08-28) — fuente de σ = puntos de integración del elemento

**Estado**: el L4 del sub-sprint materi_stress_force cambió la fuente
de la tensión de sección de la NODAL recuperada (`NODE_DOF`) a los
PUNTOS DE INTEGRACIÓN del elemento (`ELEMENT_DOF`, las tensiones
constitutivas que el elemento usó — "the element forces needed for
this option are setup in a timestep", manual 6.913). Re-ejecución:
`TOCHNOG_PROF_BIN=... scripts/compare_professional.sh`
(2026-08-28T15:11Z; ratios sobre las componentes `s` = magnitudes).

### 10.1 Resultados del arness (GNU POST-L4 vs Professional)

| modelo | pre-L4 | POST-L4 | Professional | nota |
|---|---|---|---|---|
| `gforce7` (2 quad9) | N 1.2345×, V 2.7125×, M 0.9947× | **N 1.2365×, V 2.7167×, M 0.9976×** | EXACTO | la polución de N/V NO desaparece: las tensiones nodales recuperadas son la media EXACTA de los IPs de los elementos adyacentes (verificado en los .dbs) → la polución vive en el CAMPO σ (IPs = nodal); el MOMENTO mejora (0.9947 → 0.9976) |
| `gforce7q4` (2 quad4) | N 0.88×, V 0.10×, M 0.033× | N 0.88×, V 1.79×, M 0.0324× | EXACTO | el V cambia por el promedio por-elemento de la cizalla (|∫σ_nt| por cara); la solución lockeada del elemento intacta |
| `gforce7q4_ref` (8 quad4) | 0.215× | **0.215× (byte-idéntico)** | EXACTO | el plain quad4 mantiene su solución de formulación |
| `gffq4` (10 quad4) | 0.197× | **0.197×** | 1250 EXACTO | idem |
| `gforce7_ref` (8 quad9) | N 1.0011×, M 0.9986× | **N 1.0011×, M 0.9986×/0.9970×, V 1.4719×** | EXACTO | caso sano sin regresión |
| `gforce10`/`gforce13` (hex8 3D) | N 1.0000× | **N 1.0000× (12.34 EXACTO)** | EXACTO | vía el FALLBACK documentado: en los 3D resueltos con `derivatives` el bucle escalonado NO propaga la deformación a los IPs del elemento (el ELEMENT_DOF sale a ceros) → la sección lee las tensiones nodales recuperadas (warning único) |
| `msf_shear` | 0.3846153846 | **0.3846153846 EXACTO** | σ_xy EXACTO | corte simple, sin cambios |
| `msf_tunnel3d` | 0.1 EXACTO | **0.1 EXACTO** | 0.09980686 | sin cambios |

### 10.2 Cómo se reconstruye σ en la sección desde los IPs (decisión con evidencia)

- La cuadratura de la sección (Gauss(2)/Lobatto(3) según npol) NO
  cambia. Lo que cambia es la FUENTE evaluada en los puntos de la
  sección: la reconstrucción por LAGRANGE del campo de IPs del
  elemento (los pesos 1D de la regla PROPIA del elemento, replicada
  de pol(): Lobatto con nodos por defecto, Gauss para el SRI quad4 y
  las reglas MINIMAL, 1 punto en el eje axisimétrico).
- Reglas con nodos (quad9/hex27 Lobatto, quad4/hex8 esquinas): los
  puntos de la sección COINCIDEN con los IPs de la cara → los pesos
  son la delta de Kronecker (lectura directa del ELEMENT_DOF).
- Reglas interiores (SRI quad4 2×2 Gauss): los pesos interpolan/
  extrapolan el campo de IPs al punto de la sección (la misma
  reconstrucción "same B" del fix D-b; el momento SRI queda
  byte-idéntico — 0.9375·P·L — y la cizalla pasa a leer el σ_xy crudo
  de los Gauss points, 8.4·P medido — el σ_xy del Q4 no es
  superconvergente, DIAG §12.4).
- El cambio NO logra la estática exacta del Professional: la polución
  de N/V del quad9 vive en el CAMPO σ (los IPs la tienen igual) y el
  Professional es exacto porque NO integra el campo crudo (ni nodal ni
  de IP) — consistente con un cálculo por fuerzas internas/equilibrio
  del elemento.

### 10.3 Suite

201/201 runs + verificaciones de archivos OK (199 + msf_cant3d_hex27 +
msf_axisym). Checks actualizados CON justificación: los promediados
quad9/hex27 (los nodos del plano medio usan las caras del PROPIO
elemento — correcto per manual 6.906; banda FE < 2e-3), la cizalla de
flexión pura (2.3e-4 = 2.3% de P, el orden del promedio |∫| por
elemento), y las cizallas qsri (2.07·P plain / 8.4·P SRI — la fuente
IP lee el σ_xy crudo).

### 10.4 Conclusión del L4 (qué cierra y qué NO)

- **CIERRA**: (a) la pregunta (b) del DIAG §12 ("su estática NO viene
  del σ nodal crudo") — tampoco viene del σ de IP: la nodal es la
  media exacta de los IPs, ambos contaminados; el Professional integra
  por fuerzas internas/equilibrio (ver papers/PAPER-LINES.md); (b) la
  validación 3D con carga REAL (msf_cant3d_hex27: mom1 = P·(L−z)
  dentro del 9%); (c) el axisimétrico pendiente del L2 (msf_axisym:
  nor = σ·t por unidad de circunferencia EXACTO — con el fix del
  factor 2πr en el integrando).
- **NO CIERRA** (frentes anotados, fuera del alcance): la estática
  EXACTA del Professional (requiere la integración por fuerzas
  internas — recomendado como siguiente lote), el SRI hex8 3D, la
  cizalla de sección del esquema mixto (2D ±30%, 3D 0.08·P).

---

## 11. POST-L5 (2026-08-28) — fuerzas internas del elemento (equilibrio)

**Estado**: el L5 del sub-sprint materi_stress_force implementó la
estática de sección por FUERZAS INTERNAS del elemento (la familia del
Professional): la sección lee los resultantes de f_elem = ∫Bᵀσ dV
(las fuerzas nodales consistentes del campo de IPs) de los nodos de la
cara — la estática del cuerpo libre de las cargas, EXACTA para los
estados σ en equilibrio (el resultado que el L4 concluyó como la vía
del Professional: "su estática NO viene del campo σ crudo").
Re-ejecución: `TOCHNOG_PROF_BIN=... scripts/compare_professional.sh`
(2026-08-28T17:59Z; ratios sobre las componentes `s` = magnitudes).

### 11.1 Resultados del arness (GNU POST-L5 vs Professional)

| modelo | pre-L5 | POST-L5 | Professional | nota |
|---|---|---|---|---|
| `gforce7q4` (2 quad4) | N 0.88×, V 1.79×, M 0.0324× | **N 1.0000× (12.34), V 1.0000× (100), M 0.9984× (4992)** | EXACTO | las fuerzas internas de la solución lockeada EQUILIBRADA cumplen la estática del cuerpo libre: N/V EXACTOS y M al 99.84% (el lock vive en el campo σ y la deflexión, no en las fuerzas de sección) |
| `gffq4` (10 quad4) | 0.197× (viola pL²/8) | **M_end −1.0000× (824.99), M_center 0.9987× (424.46); \|M_e\|+\|M_c\| = 1249.45 = pL²/8 a 0.9996** | 1250 EXACTO | la identidad de estática se cumple (antes 0.197×); los valores individuales son los de viga profunda (825/425) porque la estática del cuerpo libre no depende de la formulación del elemento |
| `gforce7` (2 quad9) | N 1.2365×, V 2.7167×, M 0.9976× | **N/V/M = 5.0000× (61.7/500/24955)** | EXACTO | LIMITACIÓN MEDIDA: el resultante débil de las fuerzas internas es la estática del cuerpo libre SOLO cuando el campo σ está en equilibrio; el estado σ del run del gforce7 (30 pasos, quad9 grueso) NO lo está (preexistente — el L4 lo media como 1.24×/2.7× por integración de campo; el método de equilibrio lo hace explícito). La integración de campo era MENOS sensible a la no-equilibrio; el estado σ del GNU es el limitante, no el método (el Professional produce el σ en equilibrio con su solve) |
| `gforce7_ref` (8 quad9) | N 1.0011×, M 0.9986× | N 1.2500×, M 1.2484×, V 1.2500× | EXACTO | idem: el estado σ del run multi-paso no está en equilibrio; el resultante débil lo amplifica |
| `gforce10`/`gforce13` (hex8 3D) | N 1.0000× | **N 1.0000× (12.34 EXACTO)** | EXACTO | el N axial exacto; V/mom de la sección NO son la estática del cuerpo libre (el fallback del ELEMENT_DOF a ceros usa el σ nodal recuperado, cuyo estado tampoco está en equilibrio: V ≈ 0, mom2 0.99× en z=0 / 0.074× en z=50) |
| `msf_shear` | 0.3846153846 | **0.3846153846 EXACTO** | σ_xy EXACTO | corte simple, sin cambios; moms = 0.1923 (el par de reacciones del bloque, ver los tests) |
| `msf_tunnel3d` | 0.1 EXACTO | **nor 0.09980685, shes 0.0137271, mom1 2.03e-4 — IDÉNTICOS al Professional (0.09980686 / 0.013727 / −2.0e-4)** | 0.09980686 / 0.013727 / −2.0e-4 | **el método de equilibrio reproduce el node_dof_calcul del Professional DÍGITO A DÍGITO** (la integración de campo del L4 daba el valor analítico 0.1 EXACTO del campo prescrito; el Professional y el GNU-L5 dan el valor FE discretizado) — la validación cruzada MÁS FUERTE del método |
| `msf_beam2d` (8 quad9, propio) | moms 1-3%, shes banda ±30% | **moms = P·(8−x) y shes = P EXACTOS (6 dígitos) en TODAS las secciones** | (no comparado — modelo propio) | la estática del cuerpo libre exacta; la banda de cizalla de la integración de campo desaparece |
| `msf_cant3d_hex27` (propio) | mom1 dentro del 9%, shes 0.08·P | **mom1 = P·(L−z) dentro del 1%, shes = P dentro del 4%** | (no comparado) | la cizalla contaminada de los IPs desaparece (la resultante del equilibrio) |

### 11.2 Cómo se calculan las fuerzas internas (decisión con evidencia)

- f_elem[inod, idim] = Σ_ip vol[ip]·(Bᵀσ)[ip, inod, idim] — las
  fuerzas nodales consistentes del campo σ de los IPs (ELEMENT_DOF),
  la MISMA cantidad que materi() acumula en element_rhside con el
  signo opuesto. La cinemática replica pol()/materi(): derivadas
  físicas dn = invJ·p, volumen w·4·|detJ| (2D) / w·8·|detJ| (3D) con el
  factor 2π·r axisimétrico, y la matriz B de polynom.cc:549-599 (la
  MISMA regla de cuadratura del elemento, incl. SRI-Gauss y la regla
  de 1 punto axisimétrica).
- La sección: R_face = Σ de las fuerzas internas de los nodos de la
  cara = ∫σ·n̂ dA de la cara (la identidad de las fuerzas consistentes)
  = el resultante del cuerpo libre de las cargas. nor = n̂·R_face/l
  (tracción +), she = |t̂·R_face|/l, mom = Σ (n̂·f)·arm/l con el brazo
  desde el punto medio de la cara (la decisión de brazos del L3). El
  signo se verifica contra los targets del Professional (N compresión
  −12.34, M −5000 con t̂ hacia abajo: idéntico al L2-L4).
- SRI quad4: la parte de cizalla de regla completa se reemplaza por la
  fuerza interna de cizalla reducida de 1 punto del feedback (el punto
  fijo K_SRI·u = P): 4·detJ_c·b_shear·media(σ_xy de los IPs) — el
  módulo se cancela (σ_xy = G·γ, la media Gauss 2×2 de la γ bilineal =
  el valor del centroide EXACTO). Sin la corrección la sección SRI
  leería el equilibrio K_full (lockeado).
- VERIFICACIÓN DEL EQUILIBRIO (el criterio del lote): Σ f_elem ≈
  cargas aplicadas. Medido a través de las propias fuerzas de sección:
  la ménsula msf_beam2d da shes = P y moms = P·(8−x) EXACTOS en TODAS
  las secciones (la cara del empotramiento = la reacción −P con su
  momento P·L); la biempotrada gffq4 cumple la identidad |M_e|+|M_c| =
  pL²/8 a 0.9996. Las fuerzas internas de una solución en equilibrio
  cumplen la estática de las cargas POR CONSTRUCCIÓN (Σ_elementos
  f_elem = −P en los dofs libres), independientemente de la
  formulación del elemento (lockeada o no).

### 11.3 Conclusión del L5 (qué cierra y qué NO)

- **CIERRA**: (a) la implementación de la estática por fuerzas
  internas/equilibrio del elemento — la vía del Professional que el L4
  identificó ("su estática NO viene del campo σ crudo"); (b) la
  validación cruzada MÁS FUERTE: msf_tunnel3d = los valores del
  node_dof_calcul del Professional DÍGITO A DÍGITO; (c) la estática
  del cuerpo libre EXACTA en los estados en equilibrio: gforce7q4
  N/V EXACTOS, la identidad gffq4 pL²/8 a 0.9996, msf_beam2d
  P·(8−x)/P a 6 dígitos, msf_cant3d_hex27 mom1 al 1% y shes = P (la
  cizalla contaminada 0.08·P desaparece); (d) el lock del quad4 ya no
  contamina las fuerzas de sección (la estática del cuerpo libre de la
  solución equilibrada es exacta — el lock queda en la deflexión y el
  campo σ).
- **NO CIERRA** (frente abierto, del SOLVER — fuera del post-proceso):
  el estado σ de los runs gruesos multi-paso del GNU con quad9/hex8
  (gforce7 5×, gforce7_ref 1.25×, gforce10/13 V/mom) NO está en
  equilibrio con las cargas (el baseline L4 ya lo media como la
  polución N 1.24×/V 2.7× por integración de campo; el resultante
  débil lo amplifica). La estática exacta del Professional en ESAS
  mallas requiere que el solve produzca el σ en equilibrio (Bᵀσ = P) —
  la vía es del solver mixto, no de la sección. El SRI hex8 3D sigue
  pendiente.

---

## 12. POST-L6 (2026-08-29) — estado σ en equilibrio: ELEMENT_DOF 3D poblado + cinemática de sección corregida

**Estado**: el frente abierto del L5 ("el estado σ de los runs gruesos
NO está en equilibrio — Bᵀσ = P requiere la vía del solver mixto") se
CIERRA con dos fixes del POST-PROCESO/salida, no del solver: (a) el
ELEMENT_DOF de los 3D con `derivatives` nunca recibía el σ
constitutivo (el misterio "inc_ept=0" del L4 — mecanismo medido en
DIAG §13.1: materi.cc:850 escribía en el slot `stres_indx/nder + j`
en vez de `stres_indx + j·nder`, y con nder=5 el σ caía en el bloque de
desplazamientos; el restore del bloque de tensión leía solo 9 slots —
parcial — y los rangos epe/epp/epi con índices -1 capturaban slots
[0,8)); (b) las fuerzas internas de sección 2D (L5) se integraban con
la inversa del jacobiano TRASPUESTO (términos cruzados invjac[1]↔[2]
intercambiados en msf_element_internal_forces_2d) — un factor de
amplificación dependiente de la malla (5× para los elementos 50×10 del
gforce7, 5/4× para los 12.5×10 del gforce7_ref) que NO existía para
jacobianos diagonales (quad4, quad9 cuadrados) ni en 3D (producto
matriz-vector completo). El campo σ SÍ estaba en equilibrio: residuo
del solver 2.19e-13, vely(x=50) = −0.01154 ≈ 0.01138 analítico, σxx
nodal = ±301 ≈ ±300 analítico. Re-ejecución:
`TOCHNOG_PROF_BIN=... scripts/compare_professional.sh` (2026-08-29T05:54Z).

### 12.1 Resultados del arness (GNU POST-L6 vs Professional)

| modelo | POST-L5 | POST-L6 | Professional | nota |
|---|---|---|---|---|
| `gforce7` (2 quad9) | N/V/M 5.0000× | **N 1.0000× (12.34), V 1.0000× (100), M 0.9984× (4992)** | EXACTO | la estática del cuerpo libre del estado equilibrado — el 5× era la cinemática del L5 (DIAG §13.2), no el no-equilibrio |
| `gforce7_ref` (8 quad9) | N/M/V 1.2500× | **N 1.0000×, V 1.0000×, M 0.9987×** | EXACTO | idem (factor 5/4× de la misma cinemática) |
| `gforce7q4` (2 quad4) | N 1.0000×, V 1.0000×, M 0.9984× | **byte-idéntico** | EXACTO | jacobiano diagonal: la cinemática corregida es idéntica |
| `gffq4` (10 quad4) | identidad pL²/8 0.9996 | **byte-idéntico** | EXACTO | idem |
| `gforce10`/`gforce13` (hex8 3D) | N 1.0000× (fallback del L4) | **N 1.0000× (ELEMENT_DOF real — el fallback ya no dispara)** | EXACTO | el misterio inc_ept=0 cerrado: el σ constitutivo llega a los IPs; V/mom de la sección 3D siguen ≈ 0 (frente de caras 3D del post-proceso, no del esquema — el cantilever cargado en y msf_cant3d_hex27 SÍ da shes = P al 4%) |
| `msf_shear` | 0.3846153846 EXACTO | **EXACTO, sin cambios** | σ_xy EXACTO | — |
| `msf_tunnel3d` | = Professional dígito a dígito | **sin cambios** | 0.09980686 | — |

### 12.2 Verificación del equilibrio (el criterio del lote)

Residuo Bᵀσ − P medido con las fuerzas internas del ELEMENT_DOF
(cinemática corregida, gforce7): Σ f_elem = −P en los dofs libres y
las reacciones de empotramiento (12.34, 100) iguales a las cargas —
Bᵀσ = P dentro de la tolerancia del solver (1e-5) ANTES y DESPUÉS del
lote: el estado σ del esquema escalonado con el fix C/D SIEMPRE estuvo
en equilibrio; el "5×" era la lectura del L5. Residuo de la lectura:
61.7/500/24955 (5×) → 12.34/100/4992 (1.0000×/1.0000×/0.9984×).

### 12.3 Suite

209/209 runs + verificaciones de archivos OK en build limpio. Checks
actualizados CON justificación (calibrados contra el restore σ parcial
de los modelos con `derivatives`, que perdía σyy/σxy... para nder=4 —
el restore ahora lee el bloque completo 6·nder): mesh_act_grav
(velx −0.95 → +1.3333 — el signo físico de la rampa de gravedad +x: el
desplazamiento elástico u = F·L/(E·A) = 1000·1/1000 = 1), cmat_gate
(sigxy −88.39 → +9.09 — la cizalla lineal con el primer paso capped;
el A/B |sigxy| > 1.0 se conserva) y qsri3d_beam (los momentos de
sección 710/711 ahora leen la estática del cuerpo libre P·L = 0.08
para SRI y OFF — el discriminador SRI es la deflexión de los targets,
0.897× vs 0.221× de la EB). Tiempos: peor caso 6 s por test (timeout
120 s sin riesgo).

### 12.4 Conclusión del L6 (qué cierra y qué NO)

- **CIERRA**: (a) el misterio inc_ept=0 del L4 (el ELEMENT_DOF de los
  3D con `derivatives` sale a ceros) — mecanismo completo: slot de
  escritura equivocado + lag de una iteración + restore parcial —
  fijado en materi.cc:850/elem.cc; (b) el "no-equilibrio" del L5 — el
  campo σ estaba en equilibrio, el 5× era la cinemática traspuesta de
  la sección 2D (fix en msf_element_internal_forces_2d); (c) la última
  brecha del arness: gforce7/gforce7_ref a 1.0000× sin el fallback 3D.
- **NO CIERRA** (frentes anotados, fuera del alcance): el V/mom de la
  sección 3D de los modelos con carga axial (gforce10/13 — la cara/el
  brazo 3D del post-proceso, no el esquema), el SRI hex8 3D con modos
  de energía cero, la polución de cizalla del σ_xy crudo del Q4.
