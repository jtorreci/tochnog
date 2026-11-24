# warnings_cleanup: `make audit` (-Wall -Wextra -Wpedantic) 300 → 0

Branch `dev/warnings` (worktree aislado). Objetivo: reducir los warnings del
target de auditoría al mínimo **sin ningún cambio de comportamiento**.
Verificación de neutralidad: suite interna 16/16 + spot-check de .dbs
byte-idénticos contra el binario del árbol principal.

> Actualización 2026-09-07 (dev/find): los 4 hallazgos residuales eran bugs
> reales o latentes y se corrigieron con la lógica del flujo de cada modelo
> (ver "Resolución de hallazgos" al final). El audit vuelve a **0 warnings**.

## Contadores

| Categoría | Antes | Después | Tratamiento |
|---|---|---|---|
| `-Wunused-parameter` | 129 | 0 | nombres comentados (`/*name*/`) en .cc; `(void)` en .c (C pre-C23 no admite omitir nombres) |
| `-Wunused-variable` | 66 | 0 | eliminados declaradores muertos |
| `-Wregister` (C++17) | 46 | 0 | keyword eliminada (deprecada, sin efecto) |
| `-Wunused-but-set-variable` | 30 | 0 | eliminados restos; `(void)` donde la asignación convivía con código vivo |
| `-Wformat-overflow=` | 6 | 0 | buffers ampliados (`print_fr.cc`) |
| `-Wwrite-strings` | 6 | 0 | `db_number()` pasa a `const char[]` |
| `-Wmaybe-uninitialized` | 5 | 0 | 3 falsos positivos inicializados; 2 restos reales corregidos (hallazgos 1-2) |
| `-Wunused-function` | 4 | 0 | `inv_eps`×2, `pp_kk_set`, `pert_DM` (estáticos muertos) eliminados |
| `-Wstringop-truncation` | 2 | 0 | `strncpy`→`memcpy` sobre buffer pre-rellenado |
| `-Wextra` (address of register var) | 2 | 0 | resuelto con la eliminación de `register` |
| `-Wparentheses` | 1 | 0 | paréntesis explícitos |
| `-Wvla` | 1 | 0 | VLA de `hypoplas.cc` → heap `get_new_dbl` (hallazgo 3) |
| `-Wint-in-bool-context` | 1 | 0 | `task==` añadido a test9 (hallazgo 4) |
| `-Wempty-body` | 1 | 0 | `if (db(...));` → `if (db(...)) { }` |
| **Total** | **300** | **0** | |

Los 4 restantes eran hallazgos deliberadamente no tocados por la limpieza
original (ver abajo); corregirlos exigía decidir semántica de código de
modelos constitutivos o cambiar comportamiento en un camino erróneo. Eso se
hizo en `dev/find` con la lógica del flujo real de cada modelo y verificación
contra el corpus del Professional (ver "Resolución de hallazgos").

## Qué se hizo por categoría

- **register (46)**: eliminada la keyword en `math.cc` (34) y `so_bicg.cc`
  (12+2). Es un no-op del lenguaje desde C++17; el código compilado es
  idéntico.
- **unused-variable (66)**: declaradores muertos eliminados. Verificación
  previa por variable: ningún uso residual dentro de su función, ni usos
  detrás de `#if`/`#ifdef` (se comprobó cada "extra-hit" de token contra el
  alcance de función real). En `top.cc` el inicializador `db_dbl(...)` se
  conservó como llamada suelta para no alterar el flujo de error del lookup.
- **unused-but-set (30)**: mayoría restos de refactors (p.ej. `zero_one` en
  `data.cc`, `nuknwn` en `groundda.cc`, `old_name` en `print_dx.cc`,
  `max_node_old` en `interface.cc`, `phi_flow` en la variante directa de
  stress, `smooth_size` en `print_hi.cc`, `te` en `materi.cc`). Donde la
  asignación convive con código vivo o hay riesgo de configs alternas
  (`so.cc` `ksptype`/`pctype` se usan bajo `#if PETSC_USE`) se añadió
  `(void)var;` — no-op garantizado. En modelos constitutivos
  (masin/masin_visco/sanisand) mismo criterio conservador: `(void)` en vez
  de borrar matemática muerta.
- **unused-parameter (129)**: en .cc los nombres se comentaron en la firma
  (`tipo /*name*/`); en .c (C < C23 no permite omitir nombres en
  definiciones) se restauró el nombre y se añadió `(void)` en el cuerpo.
  `umat.c` (plantilla Abaqus K&R) recibe un bloque `(void)` de los 38
  argumentos al inicio del cuerpo.
- **write-strings (6)**: `db_number()` solo lee el string (verificado:
  únicamente `strcmp`); firma `char name[]` → `const char name[]` en
  `database.cc` + **ambos** `tochnog.h` y `tochnog-mod.h` (sincronizados).
- **format-overflow (6) / stringop-truncation (2)** (`print_fr.cc`):
  ampliados `tmp[20]→[32]`, `field[16]→[24]`, `pstep[80]→[96]`,
  `r2[96]→[160]`, `compnames[6][9]→[6][12]` (los nombres FRD se siguen
  truncando a 8 con `strncpy` antes del `%sX`); `strncpy(cl,"  100CL",7)` y
  `strncpy(&cl[36],"STATIC",6)` → `memcpy` (cl está pre-rellenado de
  espacios y terminado en `cl[75]='\0'`: comportamiento idéntico).
- **unused-function (4)**: eliminadas `inv_eps` (duplicada en masin.c y
  masin_visco.c), `pp_kk_set` y `pert_DM` (sanisand.c) — estáticas y sin
  referencias (no hay `#if` en esos archivos). `pert_DM` además leía
  `y_star` sin inicializar, por lo que era código inservible.
- **parentheses (1)** (`new_mesh.cc`): paréntesis explícitos alrededor de
  los dos pares `&&` dentro del `||`.
- **empty-body (1)** (`so.cc`): `if ( db(...MG...) );` → `if (... ) { }`
  (el `;` huérfano hacía que la llamada `solve_iterative_petsc` bajo
  `#if PETSC_USE` fuera incondicional; se preserva exactamente esa semántica
  en ambas configuraciones).
- **maybe-uninitialized (3 de 5 resueltos)**:
  - `hypo.c` `c1`/`c2`: falso positivo — se asignan dentro de
    `if (*ihypotype==1)` y se usan bajo la misma condición; el compilador no
    puede probar que el puntero no cambie entre llamadas intermedias.
    Inicializados a `0.` en la declaración (camino de uso siempre pasa por
    la asignación previa → neutro).
  - `sanisand.c` `y_2`/`y_3` en `trial_state`: se escriben con `for (i<n)`
    y `n==NYDIM==20` siempre (verificado: `nyact = 6+nasvy = 20`); el
    compilador no acota `n`. Inicializados a `{0.0}` (neutro: se sobreescriben
    por completo antes de leerse).

## Hallazgos (warnings que delataban bugs o riesgo real) — CORREGIDOS

> Estado: los 4 hallazgos fueron corregidos en `dev/find` (2026-09-07) con la
> lógica del flujo real de cada modelo. El audit vuelve a **0 warnings** y la
> suite interna sigue **16/16**. Ver "Resolución de hallazgos" al final para
> la verificación por hallazgo (resultados, Pro, blast-radius).

1. **`masin_visco.c:520` — lectura de variable sin inicializar**
   `norm_tr_hypo_Dsom_rot` solo se asigna si `norm_Dsom > 1e-8` y se lee
   SIEMPRE en `wyfact = sqrt(1/3)*norm_tr_hypo_Dsom_rot + 1`. Un incremento
   de deformación desviadora ≈ 0 (compresión isotrópica, punto en reposo)
   produce un valor indeterminado que contamina `wy`/`acorrwy` y con ello el
   tangente `LL_unl`. Inicializar a un valor concreto elegiría una semántica
   física que el código no documenta → se deja el warning como señal.
2. **`sanisand.c:937` — `Kpm1` sin inicializar en camino de error**
   `plast_mod_DM` retorna sin escribir `*Kpm1` cuando `Kp < 0` en la rama
   `mario_DT_test != 0` (fija `*error = 3`). El llamador solo comprueba
   `*switch2` antes de leer `one/Kpm1`; con basura > 0 proseguiría la
   iteración de `drift_corr_DM` con un factor corrupto en un estado que ya
   señaló error. No se toca: el fix correcto exige decidir el flujo de error.
3. **`hypoplas.cc:168` — VLA con sospecha de off-by-two en el buffer**
   `double nonloc_info[length_nei-2]` se rellena/escribe con longitud
   `length_nei` (`array_set`, `db(...GET/PUT)`) y se indexa
   `[1+npointmax*ndim+npointmax] == length_nei-2`, i.e. 1-2 doubles fuera del
   buffer local cuando el modo no-local está activo (`find_local_softvar`).
   Otros sitios (`elem.cc`, `generate.cc`) usan buffers de `length_nei`. El
   arreglo neutro requiere saber si el record debe llevar `+2` y el buffer
   local `-2` (probablemente un resto de un cambio de layout); convertir el
   VLA a tamaño fijo sin resolver el off-by-two sería cosmético.
4. **`plasti.cc:1106` — typo sospechoso en la gate del Mohr-Coulomb apex**
   `test9 = GET_FLOW_RULE && plasti_type==..._APEX` (sin `task==`, a
   diferencia de test1..test8). Efecto: la gate exterior de mohrcoul se abre
   para CUALQUIER task (incl. `GET_FLOW_RULE_GRAD`) cuando
   `plasti_type==APEX`. Análisis del flujo: las gates interiores filtran por
   `task`, así que hoy no hay cambio observable en los tasks usados por
   stress.cc — pero es un bug latente si alguien añade un task nuevo.
   Corregirlo (`test9 = task==GET_FLOW_RULE && ...`) es neutro por análisis,
   pero cambia semántica aparente; se documenta sin tocar.
5. *(contexto)* `hypoplas.cc` está junto al hallazgo 3: mismo bloque.

## Resolución de hallazgos (dev/find)

Cada fix se decidió por el flujo real del código (no inicialización a ciegas
para callar el warning) y se verificó antes/después con el corpus del
Professional y la suite interna.

### 1. masin_visco.c — dirección degenerada: else con `norm_tr = 0`

**Causa exacta**: en `get_tan()` (formulación LD de Jerman-Masin 2020),
`norm_tr_hypo_Dsom_rot = tr(hypo_Dsom_rot)/norm_Dsom` solo se asigna si
`fabs(norm_Dsom) > 1e-8`. `hypo_Dsom_rot` es la dirección de carga
normalizada (línea 454); cuando degenera (norma ≈ 0, p.ej. compresión
isotrópica), el cociente no está definido y la lectura posterior en
`wyfact` usaba un valor indeterminado.

**Fix**: `else norm_tr_hypo_Dsom_rot = 0.0;` con comentario. Valor decidido
por el límite del flujo: `norm_tr = 0 → wyfact = 1 → wy = 1 → acorrwy = 0`,
i.e. la corrección direccional de la formulación LD se anula — el único
límite consistente cuando no hay dirección definida (equivalente a un tensor
de dirección puramente desviador).

**Verificación**: hypo12 (masin visco, Dr/Iv) produce `.dbs` **byte-idénticos**
antes/después (el camino `norm_Dsom ≤ 1e-8` no se dispara en ese test);
hypo13 falla idéntico antes/después (breakdown del solver iterativo,
pre-existente y documentado — calibración visco pendiente, no relacionado
con este fix). No hay test del corpus que aísle la rama degenerada; el fix
elimina el UB y fija la semántica al límite físicamente neutro.

### 2. sanisand.c — Kpm1 definido en la rama de error + guard del llamador

**Causa exacta**: `plast_mod_DM()` retorna con `*error = 3` sin escribir
`*Kpm1` cuando `Kp < 0` en la rama `mario_DT_test != 0` (último recurso del
substepping adaptativo). `drift_corr_DM()` solo comprueba `*switch2` antes
de `one/Kpm1`; con `Kpm1` basura positiva proseguía la iteración con un
factor corrupto en un estado ya señalado como error.

**Fix** (dos capas, por el flujo real):
- `plast_mod_DM()`: en la rama `Kp < 0` escribe `*Kpm1 = one/Kp` (negativo)
  antes de fijar `*error = 3`. Es el mismo valor que la rama normal habría
  escrito; con `Kp < 0` el guard del llamador (`one/Kpm1 <= 0` ⇔ `Kp <= 0`)
  lo detecta y aborta limpio con `error = 3` (el umat hace fallback a `y_n`).
- `drift_corr_DM()`: comprueba `*error == 3` justo después de
  `plast_mod_DM()` y retorna — el llamador no debe seguir iterando con un
  módulo plástico que la rutina de módulos ya declaró inválido.
- `Kpm1 = 0.0` en la declaración (solo para satisfacer el análisis estático:
  todo path que alcanza el guard ya lo escribió).

**Verificación**: no hay test sanisand en el corpus Pro (búsqueda
`GROUP_MATERI_PLASTI_SANISAND` en PROF/validation-suite: 0 hits) — la
verificación es de flujo (arriba) + audit 0 warnings + suite 16/16 + corpus
sin cambios de estado.

### 3. hypoplas.cc — buffer del record no-local a `length_nei` (heap)

**Causa exacta**: el record `NONLOCAL_ELEMENT_INFO` tiene longitud
`length_nei = 1+npointmax*ndim+npointmax+2` (database.cc `data_length` y
todos sus usuarios: elem.cc, nonloc.cc, top.cc, generate.cc). hypoplas.cc
declaraba `double nonloc_info[length_nei-2]` y lo usaba con longitud
`length_nei` (`array_set` y `db()` GET/PUT), escribiendo 1-2 doubles fuera
del buffer local. El slot `length_nei-2` (índice
`1+npointmax*ndim+npointmax`) es el flag "elemento usado para el cálculo
no-local" que top.cc resetea al cerrar el paso y nonloc.cc lee como
`jelem_nonloc` — el buffer local DEBE cubrir el record completo.

**Fix**: buffer en heap `get_new_dbl(length_nei)` + `delete[]` (mismo patrón
que elem.cc/generate.cc/top.cc para este record), eliminando el VLA
(-Wvla/-Wpedantic) y el overflow. Semántica intacta: GET → flag=1 → PUT.

**Verificación**: el bug estaba **ACTIVO** — probe no-local sintético
(hypo1 wolfersdorff + `materi_plasti_softvar_local/_nonlocal` +
`options_nonlocal_softvar`): el binario sin el fix aborta con "Error detected
for data item: nonlocal_element_info / record 0" tras "Nonlocal weights
ready"; con el fix corre rc=0 hasta el final. Wolfersdorff sin no-local
(hypo1/2/4): `.dbs` byte-idénticos antes/después (neutro). No existe test
no-local+hypo en el corpus Pro ni en la suite; el camino se ejercitó con el
probe (el Professional no acepta las keywords `materi_plasti_softvar_*` en
initia, así que no hay referencia Pro para este camino).

### 4. plasti.cc — gate del apex Mohr-Coulomb con `task==`

**Causa exacta**: `test9 = GET_FLOW_RULE && plasti_type==..._APEX` — sin
`task==`, la constante de enum `GET_FLOW_RULE` (no nula) hacía que la gate
exterior de mohrcoul se abriera para cualquier task cuando
`plasti_type == APEX`.

**Fix**: `test9 = task==GET_FLOW_RULE && plasti_type==..._APEX` (consistente
con test1..test8). Neutro por análisis: las gates interiores del bloque
mohrcoul filtran por `task` y los únicos tasks usados por stress.cc son
GET_YIELD_RULE / GET_FLOW_RULE / GET_FLOW_RULE_GRAD; para GRAD el bloque
interior no produce cambios observables.

**Verificación**: tests mohrcoul del corpus (mohrcou1/2/4, mohr_coul_direct*)
sin cambio de estado; audit 0 warnings (-Wint-in-bool-context eliminado);
suite interna 16/16.

### Blast-radius

Familia hypo (hypo1/2/3/4/6/7/8/9/10/11/12/13), slope_* (no-local clásico),
mohrcou*/mohr_coul_direct*, incremental_driver (hypo y mohr-coulomb): sin
cambios de estado (PASS/RUNFAIL/PARSE) entre baseline (commit 9609431,
binario del árbol principal) y el árbol con los 4 fixes. Corpus completo:
**204 PASS / 151 RUNFAIL / 8 PARSE** en ambos (ver corpus_find_baseline.txt
vs corpus_find_fixed.txt). Suite interna 16/16 OK. Audit `-Wall -Wextra
-Wpedantic`: **0 warnings** (antes 4).
## Verificación de neutralidad

- Build audit completo (`make audit`, con `--clean`): 0 errores, binario
  `build/tochnog-audit`.
- Suite interna del repo: **16/16 OK** (modo reducido post-incidente, tests
  hypo1-4, gforce7, tsup_* + file-checks).
- Spot-check de corpus (binario del árbol principal como referencia,
  SOLO lectura): `elasti6 interface1 conspr1 ground1 hypo1 mpc1 force1` —
  `.dbs` byte-idénticos (diff=0 tras filtrar `time_calculation`/`version`).

## Archivos

45 ficheros: 43 fuentes + `tochnog.h`/`tochnog-mod.h` (sincronizados, un
cambio puntual: firma `const` de `db_number`). Sin cambios en scripts,
`.github`, ProjectDocs (excepto este documento), ni en el árbol principal.
