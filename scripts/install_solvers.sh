#!/bin/bash
# Instala librerias optimizadas de solucion de ecuaciones para tochnog.
#
# Lo que instala:
#   - SuperLU 6.0.1 (secuencial, directo) via apt (libsuperlu-dev)
#   - OpenBLAS (BLAS/LAPACK optimizados) -- normalmente ya esta
#
# Uso (como superusuario): sudo ./scripts/install_solvers.sh
#
# NOTA: si no quieres instalar nada en el sistema, tochnog puede usar el
# SuperLU 6.0.1 compilado localmente en external-downloads/superlu-6.0.1/
# (ver scripts/build_safe.sh, ya lo detecta automaticamente).

set -e

echo "==> Actualizando indice de paquetes"
apt-get update

echo "==> Instalando SuperLU 6 (libsuperlu-dev)"
apt-get install -y libsuperlu-dev

echo "==> Verificando OpenBLAS (libopenblas0)"
apt-get install -y libopenblas0 libopenblas0-pthread

echo "==> Librerias instaladas:"
dpkg -l 2>/dev/null | grep -iE "superlu|openblas" | awk '{print "   ", $2, $3}'

echo ""
echo "Listo. Ahora puedes recompilar tochnog con:"
echo "   ./scripts/build_safe.sh"
echo ""
echo "OJO: el build usara SuperLU 6 del sistema si existe (headers en /usr/include/superlu),"
echo "y si no, el SuperLU 6.0.1 local en external-downloads/superlu-6.0.1/."
