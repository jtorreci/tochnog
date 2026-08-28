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
