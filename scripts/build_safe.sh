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
#   - Verificacion automatica: corre la suite de tests (13 tests / 47 runs)
#     y comprueba targets.
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
HIPO_TOTAL=0
# 13 tests: 12 preexistentes + familia iface_mc (1 test logico = 6 runs:
# iface_mc a/a' invarianza, iface_mc_slip b/b' invarianza, tension c, gap d)
# + iface_mc_mem (memory), + iface_mc_dil/dil_1step/num (dilatancia RF-4 y
# clamp numerico del MC) + familia 3D (iface_3d friccion 3D, iface_3d_slip
# deslizamiento libre, iface_conv3d conversion tria3->prism6,
# iface_3d_stress print 3D, iface_3d_order orden -x)
# + iface_gen/iface_gen_geom/iface_gen_geom_off (control_mesh_generate_interface)
# + iface_gen_method/iface_gen_method_gen (control_mesh_generate_interface_method)
# + iface_gen_quad9 (subdivision de caras cuadraticas 2D)
# + materi_direct/materi_direct_mc/materi_direct_auto (plasti directo normal)
# + materi_direct_visco/materi_direct_wall (_visco/_wall del plasti directo)
# + mat_rel/mat_rel_reset (materi_displacement_relative)
# + slide_axi (slide_axisymmetric)
# + reset_value_linear (control_reset_value_linear)
# + mesh_act_grav/mesh_act_grav2 (mesh_activate_gravity_time, método 1 y 2)
# + strain_settle/strain_settle_diag (strain_settlement_parameters/_diagram)
# + contact (familia contact_*: apply, penalty, plasti_friction, target).
# + groundflow_consolidate_off (groundflow_consolidation_apply -no: omite el
#   termino de divergencia material; presion queda 0 frente a 1.67 con default).
# + groundflow_vangenuchten (ley van Genuchten: gsat nodo 5 = 0.7364 y pres
#   media -2.76 vs -5 lineal por la k reducida) y groundflow_nonsaturated_off
#   (groundflow_nonsaturated_apply -no: pres media -5.0 exacta, saturado).
# + groundflow_total_pressure_tension (con epp inicial 0.01, la presion estatica
#   de water_height reemplaza a la de la ecuacion: node_rhside -4.147 vs 0.0525)
#   y groundflow_interface (permeabilidad de interfaz llena el bloque derecho
#   aislado hasta pres 2.0; sin el registro queda 0).
# + groundflow_flux_edge (groundflow_flux_edge_normal sobre geometry_line 1:
#   flux 0.1 entrando por el borde inferior da pres 2.0 en el nodo inferior;
#   con flux 0 queda 0).
# + groundflow_phreatic_multiple (groundflow_phreatic_level_multiple + _static:
#   dos columnas con niveles 3 y 1 dan pres estatica -3 y -1 en el fondo).
# + groundflow_seepage (groundflow_seepage_geometry + bounda_dof -pres 0:
#   flujo saliente por el fondo -> perfil drenado, pres media 1.0).
for t in hypo1 hypo2 hypo3 hypo4 smooth1 dof1 mlx1 vtk_dof1 gen1 genbeam1 reset1 cda1 \
         iface_mc iface_mc_1step iface_mc_slip iface_mc_slip_1step iface_mc_tension iface_mc_gap \
         iface_mc_mem iface_mc_dil iface_mc_dil_1step iface_mc_num \
         iface_3d iface_3d_slip iface_conv3d iface_3d_stress iface_3d_order \
         iface_gen iface_gen_geom iface_gen_geom_off iface_gen_method iface_gen_method_gen \
         iface_gen_quad9 \
         materi_direct materi_direct_mc materi_direct_auto \
         materi_direct_visco materi_direct_wall \
         mat_rel mat_rel_reset \
         slide_axi reset_value_linear \
         mesh_act_grav mesh_act_grav2 strain_settle strain_settle_diag \
         contact \
         groundflow_consolidate_off groundflow_vangenuchten groundflow_nonsaturated_off \
         groundflow_total_pressure_tension groundflow_interface groundflow_flux_edge \
         groundflow_phreatic_multiple groundflow_seepage; do
  HIPO_TOTAL=$((HIPO_TOTAL+1))
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
echo "==> Resumen: $HIPO_OK/$HIPO_TOTAL runs OK (13 tests: 12 preexistentes + familia iface_mc en 10 runs + familia 3D en 5 runs + familia generate_interface en 6 runs + familia materi_direct en 5 runs + materi_displacement_relative en 2 runs + slide/reset_value en 2 runs + gravity/settlement en 4 runs + contact en 1 run + groundflow_consolidate_off en 1 run + groundflow_vangenuchten/groundflow_nonsaturated_off en 2 runs + groundflow_total_pressure_tension/groundflow_interface en 2 runs + groundflow_flux_edge en 1 run + groundflow_phreatic_multiple en 1 run + groundflow_seepage en 1 run)."
echo "==> Log de compilacion completo en /tmp/tn_build_safe.log"
