#!/bin/bash
# Build seguro de tochnog para WSL2.
#
# Por que es seguro:
#   - Compilacion SECUENCIAL (make -j1): evita el OOM de la VM WSL2 que
#     tumba la sesion al lanzar muchos gcc/g++ en paralelo.
#   - ulimit -v por proceso: limita la memoria virtual de cada compilador
#     para que un archivo patologico no pueda matar la VM.
#   - timeout en el paso de compilacion: si algo se cuelga, se aborta.
#   - Link con LAPACK/BLAS real (sin libf2c, ya no es necesaria).
#   - Verificacion automatica: corre los 4 tests hypo y comprueba targets.
#
# Uso: ./scripts/build_safe.sh [--clean]
#   --clean  borra *.o y recompila todo desde cero (lento, ~6-8 min).

set -uo pipefail

REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_DIR"

# Seleccion de SuperLU: se prefiere la version mas moderna disponible.
# Orden: (1) SuperLU del sistema (apt, ej. libsuperlu-dev 6.0.1),
#        (2) SuperLU 6.0.1 local compilado, (3) SuperLU 4.3 local.
SUPERLU_A=""
SUPERLU_INC=""
if [ -f /usr/lib/x86_64-linux-gnu/libsuperlu.so ] || [ -f /usr/lib/x86_64-linux-gnu/libsuperlu.a ]; then
  SUPERLU_A="-lsuperlu"
  SUPERLU_INC="-I/usr/include/superlu"
  echo "==> Usando SuperLU del sistema (apt): /usr/lib/x86_64-linux-gnu/libsuperlu"
elif [ -f "$REPO_DIR/external-downloads/superlu-6.0.1/lib/libsuperlu.a" ]; then
  SUPERLU_A="$REPO_DIR/external-downloads/superlu-6.0.1/lib/libsuperlu.a"
  SUPERLU_INC="-I$REPO_DIR/external-downloads/superlu-6.0.1/SRC"
  echo "==> Usando SuperLU 6.0.1 local: external-downloads/superlu-6.0.1"
elif [ -f "$REPO_DIR/external-downloads/superlu-4.3/lib/libsuperlu_4.3.a" ]; then
  SUPERLU_A="$REPO_DIR/external-downloads/superlu-4.3/lib/libsuperlu_4.3.a"
  SUPERLU_INC="-I$REPO_DIR/external-downloads/superlu-4.3/SRC"
  echo "==> Usando SuperLU 4.3 local: external-downloads/superlu-4.3"
else
  echo "!! No se encontro SuperLU. Instala libsuperlu-dev (sudo ./scripts/install_solvers.sh)"
  echo "   o compila uno local en external-downloads/superlu-4.3/ o superlu-6.0.1/."
  exit 1
fi

# SQLite (optional, for tabular export). Empty by default; if libsqlite3-dev
# is installed, enable with -lsqlite3. The header is included via ALL_INCLUDE
# in the Makefile when SQLITE_USE=1 in tn_sqlite.h.
SQLITE_INC="-I/usr/include"
SQLITE_LIB="-lsqlite3"

MAKE_FLAGS=( "SYS_FILE=sysposix" "OBJ=o" "BCPP=" "VCPP="
  "COMPILER_C=gcc" "COMPILER_CPP=g++"
  "COMPILER_FLAGS=-c -O1 -Wall -D_REENTRANT $SUPERLU_INC $SQLITE_INC"
  "LINK_FLAGS_BEFORE=" )

LINK_FLAGS_AFTER="-l:liblapack.so.3 -l:libblas.so.3 $SUPERLU_A $SQLITE_LIB -lm -lpthread -o build/tochnog"

echo "==> Limite de memoria por proceso: 4 GB"
ulimit -v 4000000 2>/dev/null || echo "    (no se pudo aplicar ulimit, continuando)"

if [ "${1:-}" = "--clean" ]; then
  echo "==> Build limpio: borrando *.o"
  rm -f *.o
fi

echo "==> Compilando con make -j1 (secuencial)..."
timeout 900 make -j1 tochnog "${MAKE_FLAGS[@]}" \
  "LINK_FLAGS_AFTER=$LINK_FLAGS_AFTER" > /tmp/tn_build_safe.log 2>&1
RC=$?
if [ $RC -ne 0 ]; then
  echo "!! make fallo (rc=$RC). Ultimas 30 lineas del log:"
  tail -30 /tmp/tn_build_safe.log
  exit 1
fi
echo "==> Compilacion y link OK. Binario: $(ls -la build/tochnog | awk '{print $5}') bytes"

echo "==> Verificando ausencia de f2c en el binario..."
F2C_SYM=$(nm build/tochnog 2>/dev/null | grep -cE "s_wsle|pow_dd|s_stop|do_lio")
echo "    simbolos f2c: $F2C_SYM (debe ser 0)"
if [ "$F2C_SYM" != "0" ]; then
  echo "!! OJO: el binario aun referencia runtime f2c"
fi

echo "==> Ejecutando tests hypo con limites de memoria..."
HIPO_OK=0
for t in hypo1 hypo2 hypo3 hypo4 smooth1 dof1 mlx1 vtk_dof1 gen1 genbeam1; do
  ( cd validation-suite/test-2014 &&
    ulimit -v 4000000 &&
    timeout 120 "$REPO_DIR/build/tochnog" "$t.dat" > "/tmp/${t}_safe.out" 2>&1 )
  RC=$?
  if [ $RC -eq 0 ]; then
    HIPO_OK=$((HIPO_OK+1))
    echo "    $t: OK (exit 0)"
  else
    echo "    $t: FALLO (rc=$RC)"
  fi
done
echo "==> Resumen: $HIPO_OK/10 tests hypo en verde."
echo "==> Log de compilacion completo en /tmp/tn_build_safe.log"
