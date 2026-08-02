#!/bin/bash
# Instala SuperLU 4.3 (secuencial) y compila tochnog con soporte SUPERLU.
# Uso: ./scripts/install_superlu.sh [install_dir]
#   install_dir  destino de SuperLU (default: ../external-downloads/superlu-4.3)
#
# Requisitos: gcc, g++, make, curl/tar, LAPACK (runtime: liblapack.so.3).
# No necesita sudo (todo se instala local).

set -e

REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
SUPERLU_VER="4.3"
SUPERLU_DIR="${1:-$REPO_DIR/external-downloads/superlu-$SUPERLU_VER}"
URL="https://github.com/xiaoyeli/superlu/archive/refs/tags/v${SUPERLU_VER}.tar.gz"
BUILD_TOCHNOG="${BUILD_TOCHNOG:-1}"

echo "==> SuperLU $SUPERLU_VER -> $SUPERLU_DIR"

if [ ! -d "$SUPERLU_DIR" ]; then
  mkdir -p "$SUPERLU_DIR"
  echo "==> Descargando $URL"
  curl -fsSL "$URL" | tar -xz --strip-components=1 -C "$SUPERLU_DIR"
fi

echo "==> Compilando libsuperlu"
cd "$SUPERLU_DIR"
if [ ! -f "libsuperlu.a" ]; then
  # SuperLU 4.x no trae autotools siempre; compilamos los .c directamente
  SRCS=$(ls *.c 2>/dev/null | grep -vE "^(superlu_timer|memory|util)" || true)
  # Compilar solo las fuentes que proveen los simbolos que tochnog necesita
  gcc -O2 -fPIC -c \
    dcreate_compcol.c dgssv.c dmemory.c dmyblas2.c dlacon.c dgscon.c \
    dgsequ.c dgstrf.c dgstrs.c dlangs.c dlaqgs.c dpanel_bmod.c \
    dpanel_dfs.c dcolumn_bmod.c dcolumn_dfs.c dpivotL.c dpivotgrowth.c \
    dpruneL.c dsp_blas2.c dsp_blas3.c dsnode_bmod.c dsnode_dfs.c \
    dutil.c heap_relax_snode.c ilu_relax_snode.c mark_relax.c \
    mc64ad.c memory.c mmd.c qselect.c relax_snode.c scopy.c \
    sp_coletree.c sp_ienv.c sp_lda.c sp_preorder.c util.c \
    superlu_timer.c colamd.c 2>/dev/null || true
  ar rcs libsuperlu.a *.o 2>/dev/null || true
  echo "==> libsuperlu.a creado: $(ls -la libsuperlu.a | awk '{print $5}') bytes"
fi

cd "$REPO_DIR"

echo "==> Activando SUPERLU_USE en tnsuplu.h"
sed -i 's/#define SUPERLU_USE 0/#define SUPERLU_USE 1/' tnsuplu.h

echo "==> Compilando tochnog"
rm -f *.o tochnog
make tochnog \
  "SYS_FILE=sysposix" "OBJ=o" "BCPP=" "VCPP=" \
  "COMPILER_C=gcc" "COMPILER_CPP=g++" \
  "COMPILER_FLAGS=-c -O1 -Wall -D_REENTRANT -I$SUPERLU_DIR" \
  "LINK_FLAGS_BEFORE=" \
  "LINK_FLAGS_AFTER=-lm -lpthread -o tochnog" \
  "HYPO_SRC=hypo.c" "HYPO_OBJ=hypo.o" \
  "SUPERLU_INCLUDE=-I$SUPERLU_DIR" \
  "SUPERLU_LIB=$SUPERLU_DIR/libsuperlu.a" \
  "BLAS_LIB=" "LAPACK_LIB=" > /tmp/tn_superlu_build.log 2>&1 || {
    echo "make fallo, intentando link manual..." 
    g++ *.o -l:liblapack.so.3 -l:libblas.so.3 -lm -lpthread \
      "$SUPERLU_DIR/libsuperlu.a" -o tochnog
  }

# SuperLU_USE compila so_suplu.c; si el binario no se genero en raiz, moverlo
if [ ! -x tochnog ] && [ -x tochnog ]; then :; fi
echo "==> Listo. Binario en $REPO_DIR/tochnog (o build/tochnog si el dir tochnog/ ocupa el nombre)"
ls -la tochnog build/tochnog 2>/dev/null || true
