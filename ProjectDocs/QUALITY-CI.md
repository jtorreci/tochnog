# Calidad / CI — objetivos make audit/asan/ubsan + GitHub Actions

Estado: implementado (PLAN.md punto 2, DECISIONES.md "Build/Calidad").

## Objetivos make (makefile, sección linux-gcc)

Cada objetivo delega en `scripts/build_safe.sh` (la MISMA resolución de
dependencias y el MISMO runner de la suite del build canónico) y SOLO
añade flags de compilación y de link vía `TN_EXTRA_FLAGS`:

| Objetivo  | Flags añadidos                                                        | Binario resultante      |
|-----------|-----------------------------------------------------------------------|-------------------------|
| `make audit` | `-Wall -Wextra -Wpedantic` (sin sanitizers)                        | `build/tochnog-audit`   |
| `make asan`  | `-fsanitize=address,undefined -fno-omit-frame-pointer -g -O1` (ASan+UBSan juntos, decisión PLAN.md) | `build/tochnog-asan` |
| `make ubsan` | `-fsanitize=undefined -fno-omit-frame-pointer -g -O1`              | `build/tochnog-ubsan`  |

Los tres recompilan TODO desde cero (build limpio) y dejan el árbol SIN
`*.o`: tras un `make audit/asan/ubsan`, para volver al build por defecto
corre `./scripts/build_safe.sh --clean` (solo las builds limpias son
fiables en este árbol, ver AGENTS.md). El binario por defecto
`build/tochnog` nunca se pisa: los targets de calidad escriben en
`build/tochnog-<modo>`.

El código es legacy (C-style C++, `.cc`): `-Wextra -Wpedantic` producen
MUCHOS warnings esperables. No se corrigen en el target: el objetivo es
que el build COMPILE y que los warnings queden visibles. El makefile
canónico ya compila con `-Wall` (base de `build_safe.sh`); `audit` añade
`-Wextra` y `-Wpedantic` sobre esa base.

### Verificación

```bash
# 1) Build de calidad + suite interna (16 tests recuperados):
make audit        # resumen de warnings al final; log en /tmp/tn_build_safe.log
make asan         # rc=0 cuando la suite pasa; HOY rc=1 esperado: hypo1-4
                  # abortan por el OOB de hypo.c (ver "Hallazgo conocido")
make ubsan
# Warnings completos:
grep -iE "warning:" /tmp/tn_build_safe.log

# 2) Smoke commitado (también es lo que corre CI):
.github/ci/run_smoke.sh build/tochnog-audit
.github/ci/run_smoke.sh build/tochnog-asan     # con ASAN_OPTIONS/UBSAN_OPTIONS
.github/ci/run_smoke.sh build/tochnog-ubsan

# 3) No-regresión (build por defecto intacto):
./scripts/build_safe.sh --clean                # 16/16, rc=0
```

Los tres targets corren la suite interna completa (16 tests recuperados
del modo reducido: hypo1-4, gforce7, tsup_*) con su binario. `make asan`
pasa a los tests `TN_MEMLIMIT_KB=unlimited` y `ASAN_OPTIONS=detect_leaks=0`.

## Hallazgo conocido (ASan): OOB de lectura en el path hypo — CORREGIDO (2026-09-06)

`make asan` compilaba y corría la suite, pero la suite NO pasaba completa:
**hypo1-4 abortaban con `stack-buffer-overflow`** (el resto, 12/16, corría
limpio bajo ASan+UBSan; sin reportes UBSan en la suite). El gate honesto
del modo reducido hacía fallar el build con "runs fallidos 12/16".

Causa raíz (bug de ABI pre-existente del fork, corregido en
`<FIX_COMMIT>`):

- `hypoplas.cc` declaraba `int find_local_sv[1], options_nonlocal[1]`
  (4 bytes) y los pasaba a `hypo_`; el resto de argumentos INTEGER del
  f2c-port (nhis, ndata, use_pres, use_epi, hypo_type) ya eran
  `long int *` (la convención f2c del árbol: `typedef long int integer`
  en tochnog.h; hypo.c es un port f2c→C puro).
- `hypo.c` (port puro C) declara los parámetros como
  `long int *find_local_sv, *options_nonlocal` (8 bytes en LP64) y lee
  `options_nonlocal[0]` como load de 8 bytes (hypo.c:165) → leía 4 bytes
  más allá del objeto en el stack (READ OOB determinista, expuesto por
  ASan; también sigma_ lee ambos flags en hypo.c:620-621).

En el build normal la lectura OOB caía en el layout del stack y
"funcionaba de casualidad" (podía incluso leer basura en el flag,
cambiando la rama `if(options_nonlocal[0] && !find_local_sv[0])`); ASan
la exponía. Fix aplicado: declarar los dos flags como `long int` en
hypoplas.cc (declaración extern y locales), consistente con los demás
argumentos INTEGER del port y con hypo.c — hypo.c NO se tocó.

**Verificación post-fix**: `make asan` → suite **16/16** (hypo1-4 ya no
abortan); corrida directa del hypo1 del corpus y del suite bajo
`build/tochnog-asan` sin NINGÚN reporte AddressSanitizer; la suite del
build normal 16/16 y los .dbs de hypo1-4/7-9/12 byte-idénticos antes y
después del fix (la basura que leía el OOB era padding a cero en el build
-O1: el bug NO contaminaba resultados, pero era UB determinista y mataba
los builds sanitizados; además PODRÍA dispararse con otras flags/layouts).
Detalle del fix y re-medición hypo vs el Professional:
ProjectDocs/SEGUIMIENTO-CONVERGENCIA.md (registro 2026-09-06).

## Variables de hooks en scripts/build_safe.sh (aditivas)

Con los defaults el script es byte-idéntico al build canónico:

- `TN_EXTRA_FLAGS` — flags extra añadidos a CADA compilación gcc/g++ y al
  link final (los sanitizers necesitan sus flags en el link).
- `TN_BIN` — nombre del binario bajo `build/` (default `tochnog`); la
  suite interna corre contra `build/$TN_BIN`.
- `TN_MEMLIMIT_KB` — tope de memoria virtual por proceso (default
  `4000000`). Los binarios sanitizados necesitan `unlimited` (ver
  Gotchas).
- `TN_SKIP_SUITE=1` — compila + check f2c y omite la suite interna (lo usa
  CI: la suite no está versionada).

## Dependencias de terceros y sanitizers (decisión documentada)

Las dependencias de terceros (SuperLU local `external-downloads/superlu-6.0.1/`,
numlib runtime `external-downloads/numlib-runtime/`, LAPACK/BLAS/gfortran
del sistema) entran al build como LIBRERÍAS PRE-COMPILADAS: nunca se
recompilan con los flags del target, así que no hace falta excluirlas del
sanitize. El link del binario sanitizado con `-fsanitize=address,undefined`
intercepta sus llamadas a malloc/free (ASan) sin instrumentar su interior.

Si en el futuro alguien compila SuperLU (u otra dependencia) DENTRO del
árbol con los flags del target, ahí sí se necesita excluirla: compilar la
dependencia con `gcc -c ...` SIN los flags de sanitize y linkar su
`libsuperlu.a` limpia. No hay objeto de terceros que recompilar hoy, por
eso no existe un interruptor `SANITIZE_DEPS` cableado; la variable natural
de escape es apuntar `SUPERLU_A` a una librería compilada sin sanitize.

## GitHub Actions (.github/workflows/ci.yml)

Disparos: push a `master`/`documentation-improvement`, pull_request y
workflow_dispatch.

- **Job `build-test`** (ubuntu-22.04): instala dependencias vía apt
  (`libsuperlu-dev liblapack-dev libblas-dev libsqlite3-dev libgfortran5`),
  compila con `TN_SKIP_SUITE=1 ./scripts/build_safe.sh --clean` y corre el
  smoke commitado `.github/ci/run_smoke.sh build/tochnog`.
- **Job `sanitize`** (opcional): compila con los flags de `make asan`
  (mismos valores) y corre el smoke bajo ASan/UBSan con
  `ASAN_OPTIONS=detect_leaks=0` y `UBSAN_OPTIONS=halt_on_error=1` (cualquier
  UB reportada hace fallar el job). Se activa en push a `master`, dispatch
  manual o PR etiquetada `sanitize` (el build sanitizado es caro; no corre
  en cada push de rama).
- Los jobs NO corren el corpus grande (~363 tests): demasiado lento. En CI
  solo corre el smoke (suite propia no versionada — ver abajo).

### Por qué CI corre un "smoke" y no la suite propia

`validation-suite/` y `external-downloads/` NO están versionados
(.gitignore: "Downloads y suites de terceros (no se versionan)"), así que
un checkout limpio de CI no tiene los `.dat` de la suite ni el SuperLU
local. La suite propia (16 tests recuperados del modo reducido) se corre
localmente con `./scripts/build_safe.sh`; en CI se corre el subset
commitado `.github/ci/smoke/` (2 inputs cuad4 elásticos: 1 y 4 elementos)
que ejercita el pipeline core: parse → malla → materi elástico → bounda →
solve → dump de base de datos. El runner `.github/ci/run_smoke.sh`
comprueba rc=0, que el `.dbs` se escribió y que contiene dofs nodales
resueltos + `end_data`.

`ubuntu-22.04` está fijado a propósito: su `libsuperlu-dev` es SuperLU
6.0.1, la misma versión que el local de desarrollo
(`external-downloads/superlu-6.0.1`); el SuperLU de Ubuntu 24.04 es otra
major y no está verificado contra este código.

## Gotchas

1. **ASan aborta bajo `ulimit -v`** (error medido:
   `AddressSanitizer failed to allocate ... ReserveShadowMemoryRange
   failed ... Perhaps you're using ulimit -v`): ASan reserva TB de shadow
   address space. El build canónico limita la memoria virtual a 4 GB por
   proceso; los targets `asan`/`ubsan` y el job sanitize usan
   `TN_MEMLIMIT_KB=unlimited`.
2. **`*.o` son compartidos por todos los builds en el mismo árbol**: un
   build de calidad recompila todo con sus flags y luego borra los `.o`.
   El siguiente build por defecto DEBE ser `--clean` para no mezclar
   objetos instrumentados con flags normales.
3. **El target `tochnog` del makefile no produce un archivo `tochnog` en
   la raíz** (`-o` lo fija el invocador): GNU make relinka cada
   invocación; no es un error.
4. **LeakSanitizer off**: `tochnog` es legacy con fugas conocidas (los
   fixes de db_close son trabajo en curso, DECISIONES.md). El job de
   sanitize caza errores de memoria (OOB, use-after-free) y UB; la caza
   de leaks es un ejercicio aparte.
5. **UBSan**: por defecto reporta y continúa (rc=0). La suite local puede
   correr con reportes; el job sanitize de CI eleva a fallo solo sobre el
   smoke (verificado limpio bajo `halt_on_error=1`).
6. **Tiempos**: cada build de calidad es una compilación completa
   secuencial (~6-10 min). Planificar los jobs de CI con timeout ≥ 30 min.
