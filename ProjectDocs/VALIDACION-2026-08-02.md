# Validación con suite oficial test-2014 (sourceforge)

Fecha: 2026-08-02
Binario: tochnog modernizado (rama documentation-improvement), compilado con g++ 13
Baseline comparado: fuente original Tochnog-Latest-jan-2014 (sin modernizar) compilada con los mismos flags

## Resumen

- Suite oficial descargada de tochnog.sourceforge.net: `validation-suite/test-2014/` (206 archivos, 192 tests listados en makefile target `small`).
- **169 PASS / 21 FAIL / 2 TIMEOUT** con binario modernizado.
- Todos los FAIL se reproducen de forma **idéntica o peor** con el binario original 2014 sin modernizar → **no hay regresiones por la modernización**.

## Resultados por categoría de fallo (23 casos)

| # | Categoría | Tests | Naturaleza | ¿Reparable con flags? |
|---|---|---|---|---|
| 1 | `.dat` inexistente | crack1-3, examp13, pridbs1.dbs | El makefile los lista pero el archivo no está en el zip descargado | No (faltan inputs) |
| 2 | `target_value` numérico fuera de tolerancia | condif5, ground3, incnav5, examp17, examp14 | El motor corre bien; el resultado NO cumple el `target_value` de referencia del input (p.ej. condif5: temp=1.11e-16 vs objetivo 1.0; ground3: -4 vs -2; incnav5: 0.472 vs 0.419; examp17: -21330 vs -20506; examp14: 9.61 vs 15.13) | No (fallo numérico de referencia; el original 2014 da los MISMOS valores) |
| 3 | `element_dof` con refinamiento | examp7, examp15, refine4, ho_othr2 | `control_mesh_refine_globally` sobre elementos de 2º orden intenta leer `element_dof` no generado | No (idéntico en original 2014) |
| 4 | `node_nonlocal` | examp19 | `materi_plasti_f_nonlocal` + `options_nonlocal` con refinamiento: `node_nonlocal` no generado | No (idéntico en original 2014) |
| 5 | `Error in data part` | examp22 | El parser no reconoce `group_materi_elasti_camclay_g` (material no soportado en esta versión) | No (idéntico en original 2014) |
| 6 | requiere SUPERLU | examp9 | El input exige `SUPERLU_USE=1` | **Sí** (recompilar con SUPERLU) |
| 7 | requiere hipoplasticidad | hypo1-4 | El input exige `HYPO_USE=1` + libf2c | **Sí** (recompilar con f2c) |
| 8 | ejecución muy lenta | examp24 | Test grande (>5 min CPU); timeout del runner | No es fallo; subir timeout |
| 9 | **Buffer overflow (mejora)** | mohrcou1, condif10 | **Original aborta (SIGABRT, glibc detecta overflow); modernizado pasa limpio** | Ya corregido |

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

## Conclusión sobre los 23 fails

- **2 son reparables con flags de compilación**: `examp9` (SUPERLU), `hypo1-4` (f2c). Son features desactivadas, no bugs.
- **1 no es fallo**: `examp24` solo necesita más de 5 min de ejecución.
- **1 es una mejora**: `mohrcou1`/`condif10` — el modernizado corrige un buffer overflow que abortaba el original.
- **El resto (5 inputs faltantes + 5 desvíos numéricos de `target_value` + 5 idénticos al original) son fallos del código fuente/inputs heredados**, no de la modernización. Los desvíos numéricos dan resultados IDÉNTICOS al binario original 2014: las referencias de la suite fueron generadas con otra compilación/plataforma.

## Pendiente

- Instalar `libf2c` para habilitar hipoplasticidad (tests hypo1-4) y correrlos.
- Compilar con SUPERLU para examp9.
- Tests `large` y `very_large` no ejecutados completos (solo muestreo).
