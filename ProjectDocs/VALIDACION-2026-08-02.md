# Validación con suite oficial test-2014 (sourceforge)

Fecha: 2026-08-02
Binario: tochnog modernizado (rama documentation-improvement), compilado con g++ 13
Baseline comparado: fuente original Tochnog-Latest-jan-2014 (sin modernizar) compilada con los mismos flags

## Resumen

- Suite oficial descargada de tochnog.sourceforge.net: `validation-suite/test-2014/` (206 archivos, 192 tests listados en makefile target `small`).
- **169 PASS / 23 FAIL** con binario modernizado (inyectando `number_of_integration_points 27` en inputs que requieren elementos de 2º orden).
- Todos los FAIL se reproducen de forma **idéntica o peor** con el binario original 2014 sin modernizar → **no hay regresiones por la modernización**.

## Resultados por categoría de fallo

| Categoría | Tests | Causa |
|---|---|---|
| `.dat` inexistente en suite | crack1-3, examp13, pridbs1.dbs | El makefile los lista pero el archivo no está en el zip |
| `element_dof` no activo | examp7, examp15, refine4, ho_othr2 | Idéntico en binario original 2014 |
| `node_nonlocal` | examp19 | Idéntico en binario original 2014 |
| `Error in data part` | examp22 | Idéntico en binario original 2014 |
| fallo sin mensaje claro | condif5, ground3, incnav5, examp17, examp14 | Idéntico en binario original 2014 |
| requiere SUPERLU | examp9 | Feature desactivada (`SUPERLU_USE=0` en tnsuplu.h) |
| requiere hipoplasticidad/f2c | hypo1-4 | `HYPO_USE=0`; falta libf2c en el sistema |
| **Buffer overflow (mejora)** | mohrcou1 | **Original aborta (SIGABRT, glibc detecta overflow); modernizado pasa limpio** |

## Nota destacada

`mohrcou1` aborta con **buffer overflow detectado por glibc** en el binario original 2014 (rc=134), mientras que el binario modernizado lo ejecuta correctamente (rc=0). La modernización (RAII, strings seguros, vectores) eliminó un desbordamiento real.

`condif10` en el original 2014 también aborta (rc=134, SIGABRT); el modernizado sale limpio (rc=1).

## Comparación de salidas numéricas

Tests muestreados (truss1, elasti1, spring1, wave1, viscos1, absorp1): el archivo `.dvd` de salida es **byte-idéntico** entre binario original y modernizado (diff_lines=0, ignorando fecha/versión/CPU time).

## Cómo se compiló

```
make tochnog "SYS_FILE=sysposix" "OBJ=o" "BCPP=" "VCPP=" \
  "COMPILER_C=gcc" "COMPILER_CPP=g++" \
  "COMPILER_FLAGS=-c -O1 -Wall -D_REENTRANT" \
  "LINK_FLAGS_BEFORE=" "LINK_FLAGS_AFTER=-lm -lpthread -o tochnog"
# + link manual: g++ *.o -l:liblapack.so.3 -l:libblas.so.3 -lm -lpthread
# hipoplasticidad desactivada: tnhypo.h HYPO_USE 0 + hypo_dum.c (stub)
```

Notas de build:
- `tochnog.h` no tenía include guards → se añadieron (`#ifndef TOCHNOG_H`).
- `f2c.h` define macros `min`/`max` que rompen `<algorithm>`/`<random>` → `#undef` tras incluirlo en `tochnog.h`.
- `string_utils.h` usaba `std::min`/`std::unique_ptr` sin incluir `<algorithm>`/`<memory>` → corregido.
- LAPACK se enlaza contra librerías versionadas (`-l:liblapack.so.3 -l:libblas.so.3`); no hay symlinks dev.

## Pendiente

- Instalar `libf2c` para habilitar hipoplasticidad (tests hypo1-4) y correrlos.
- Compilar con SUPERLU para examp9.
- Tests `large` y `very_large` no ejecutados completos (solo muestreo).
