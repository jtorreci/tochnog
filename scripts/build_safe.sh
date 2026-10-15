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
#   - Verificacion automatica: corre la suite de tests (lista en el bucle)
#     y comprueba targets.
#
# Uso: ./scripts/build_safe.sh [--clean]   (TN_CHECK_ONLY=1: sin build)
#   --clean  borra *.o y recompila todo desde cero (lento, ~6-8 min).
#
# Builds de calidad (ver makefile): make audit / asan / ubsan los invocan
# con las variables TN_EXTRA_FLAGS / TN_BIN / TN_MEMLIMIT_KB. Ejemplos:
#   make audit   # -Wall -Wextra -Wpedantic -> build/tochnog-audit
#   make asan    # -fsanitize=address,undefined -> build/tochnog-asan
#   make ubsan   # -fsanitize=undefined -> build/tochnog-ubsan

set -uo pipefail

# Las comprobaciones numericas con awk/mawk dependen de la conversion
# string->numero de strtod, que respeta LC_NUMERIC: con una locale de
# ',' decimal (p.ej. es_ES.UTF-8) "0.5" se convierte en 0 y los checks
# fallan espuriamente. Fijar el locale C hace el parseo deterministico.
export LC_ALL=C

# ---------------------------------------------------------------------------
# Quality/CI hooks (ADITIVOS: con los defaults este script es byte-identico
# al build canonico). Los usan los targets `make audit/asan/ubsan` y el
# workflow de GitHub Actions (.github/workflows/ci.yml):
#   TN_EXTRA_FLAGS  flags EXTRA anadidos a cada compilacion gcc/g++ Y al
#                   link final (los sanitizers necesitan sus flags al link).
#   TN_BIN          nombre del binario de salida bajo build/ (default tochnog).
#                   Targets de calidad: build/tochnog-audit|asan|ubsan, asi el
#                   build/tochnog por defecto nunca se pisa.
#   TN_MEMLIMIT_KB  tope de memoria virtual por proceso (compiladores y runs
#                   de tests). Los binarios sanitizados NECESITAN "unlimited":
#                   ASan reserva TB de shadow address space y aborta bajo
#                   RLIMIT_AS ("ReserveShadowMemoryRange failed... ulimit -v").
#   TN_SKIP_SUITE   =1 compila + check f2c y SALTA la suite interna (lo usa
#                   CI: validation-suite/ no esta versionado en el repo; CI
#                   corre sus propios smoke inputs en .github/ci/smoke).
#   TN_CHECK_ONLY   =1 SALTA la compilacion y corre la suite sobre el
#                   build/$TN_BIN existente (falla si no existe).
# ---------------------------------------------------------------------------
TN_EXTRA_FLAGS="${TN_EXTRA_FLAGS:-}"
TN_BIN="${TN_BIN:-tochnog}"
TN_MEMLIMIT_KB="${TN_MEMLIMIT_KB:-4000000}"
TN_SKIP_SUITE="${TN_SKIP_SUITE:-0}"
TN_CHECK_ONLY="${TN_CHECK_ONLY:-0}"

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

# SQLite (optional, for tabular export). The header may come from the
# system (libsqlite3-dev) or from external-downloads/numlib-runtime/include
# (header extracted from a .deb for machines without the dev package);
# -l: uses the runtime soname directly so the dev symlink is not required.
SQLITE_INC="-I/usr/include"
if [ -f "$REPO_DIR/external-downloads/numlib-runtime/include/sqlite3.h" ] && \
   [ ! -f /usr/include/sqlite3.h ]; then
  SQLITE_INC="-I$REPO_DIR/external-downloads/numlib-runtime/include"
  echo "==> Header sqlite3.h del runtime numerico local"
fi
SQLITE_LIB="-l:libsqlite3.so.0"

# Runtime numerico local (opcional): maquinas sin liblapack3/libblas3 en el
# sistema y sin sudo. external-downloads/numlib-runtime/ contiene las libs
# extraidas de .deb de Debian (ver su README.txt). Si existe, se anade al
# link (-L + -rpath-link) y se exporta LD_LIBRARY_PATH para los tests.
NUMLIB_RUNTIME="$REPO_DIR/external-downloads/numlib-runtime"
NUMLIB_LINK=""
NUMLIB_FORTRAN=""
if [ -d "$NUMLIB_RUNTIME" ]; then
  NUMLIB_LINK="-L$NUMLIB_RUNTIME -Wl,-rpath-link,$NUMLIB_RUNTIME"
  NUMLIB_FORTRAN="-l:libgfortran.so.5"
  export LD_LIBRARY_PATH="$NUMLIB_RUNTIME${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
  echo "==> Runtime numerico local: $NUMLIB_RUNTIME"
fi

# Optimisation -O2 (F4 T13, 2026-10-10, BACKLOG C49): -O1 until the
# fixed-size scratch-buffer overflows were fixed (scripts/audit_db_buffers.py,
# ProjectDocs/manual-developer/db_buffer_contract.md); the -O2 corpus A/B
# gives byte-identical .dbs, 29-47 % less time on large1/force11 and 8-16 %
# on the corpora (BACKLOG C49). TN_EXTRA_FLAGS is
# appended, so a later -O1 / -O0 in it still wins (sanitizer targets).
MAKE_FLAGS=( "SYS_FILE=sysposix" "OBJ=o" "BCPP=" "VCPP="
  "COMPILER_C=gcc" "COMPILER_CPP=g++"
  "COMPILER_FLAGS=-c -O2 -Wall -D_REENTRANT $SUPERLU_INC $SQLITE_INC $TN_EXTRA_FLAGS"
  "LINK_FLAGS_BEFORE=" )

LINK_FLAGS_AFTER="$NUMLIB_LINK -Wl,--start-group $SUPERLU_A -l:liblapack.so.3 -l:libblas.so.3 $NUMLIB_FORTRAN -Wl,--end-group $SQLITE_LIB -lm -lpthread $TN_EXTRA_FLAGS -o build/$TN_BIN"

echo "==> Limite de memoria por proceso: $TN_MEMLIMIT_KB"
ulimit -v "$TN_MEMLIMIT_KB" 2>/dev/null || echo "    (no se pudo aplicar ulimit, continuando)"

# build/ puede no existir en un checkout fresco (CI, clon nuevo)
mkdir -p build

if [ "$TN_CHECK_ONLY" = "1" ]; then
  # Check-only mode: no build, run the suite on the existing binary.
  if [ ! -x "build/$TN_BIN" ]; then
    echo "!! TN_CHECK_ONLY=1 but build/$TN_BIN does not exist"
    exit 1
  fi
  echo "==> TN_CHECK_ONLY=1: build skipped, using the existing build/$TN_BIN"
else
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
  echo "==> Compilacion y link OK. Binario: $(ls -la build/$TN_BIN | awk '{print $5}') bytes"
  if [ -n "$TN_EXTRA_FLAGS" ]; then
    echo "==> Flags extra aplicados: $TN_EXTRA_FLAGS"
    NWARN=$(grep -ciE "warning:" /tmp/tn_build_safe.log || true)
    echo "==> Warnings en el log de compilacion: $NWARN (log: /tmp/tn_build_safe.log)"
  fi

  echo "==> Verificando ausencia de f2c en el binario..."
  F2C_SYM=$(nm build/$TN_BIN 2>/dev/null | grep -cE "s_wsle|pow_dd|s_stop|do_lio")
  echo "    simbolos f2c: $F2C_SYM (debe ser 0)"
  if [ "$F2C_SYM" != "0" ]; then
    echo "!! OJO: el binario aun referencia runtime f2c"
  fi
fi

if [ "$TN_SKIP_SUITE" = "1" ]; then
  echo "==> TN_SKIP_SUITE=1: suite interna omitida (solo compilacion)."
  echo "==> Log de compilacion completo en /tmp/tn_build_safe.log"
  exit 0
fi

echo "==> Ejecutando tests hypo con limites de memoria..."
HIPO_OK=0
HIPO_TOTAL=0
# Sprint 11 lote 1 (control_print_*): los archivos de salida que se
# APPENDEN entre ejecuciones se limpian antes del bucle para que las
# verificaciones de la seccion posterior partan de archivos frescos.
# Sprint 11 lote 2 (frecuencia de prints): dof.* (freq_timeint usa
# dof.0..dof.12 y freq_timestep dof.22/23; dof.20/21 se regeneran en el
# bucle), los .frd de print_frd -separate_sequential (uno por print) y
# sigxx4.his (historia de los 2 tests, se suma).
# Sprint 11 lote 3 (extensiones VTK): tn30..tn43.vtk de vtk_coord1/
# vtk_dofcalc1/vtk_empty1/vtk_nodmeth1/vtk_other1.
# Sprint 11 lote 4 (dof_line/dof_point): disy.*/velx.*/vely.*/disx.* de
# dpline1/dpline_n/dpline_geom/dpline_group/dpline_eps/dpline_method/
# dpline_move/dpline_time/dpoint_time/dpoint1 y adis.* (calcul de dpoint1).
# Sprint 11 lote 5 (node prints + dof smoothing): velx.10/11/20/30/40/52,
# vely.10/11/20/21/54/55/56/57/50/51 (cpn1/cpn_angular/cpn_geom/cpn_sort/
# cpn_zero; velx.*/vely.* ya cubiertos), node_dof_0.11/node_dof_1.11
# (cpn1 numerico), avel.12 (calcul de cpn1) y dof.60/61/70/71
# (dsmooth1/dsmooth_n; dof.* ya cubierto).
# Sprint 11 lote 6 (beam force/moment): beam_force_moment.5/7/8 (bmom1/
# bmom_truss/bmom_2d), beam_force_moment.0 (bmom_switch, sequential) y
# element_truss_force_0.7 (bmom_truss; bmom_noint NO debe crear archivo).
# Sub-sprint materi_stress_force lote 1 (infraestructura): los archivos
# materi_stress_force.* de msf_parse/msf_print/msf_errors_2dwarn (el
# print escribe UN bloque por step_close en modo append).
rm -f validation-suite/test-2014/dof.* \
      validation-suite/test-2014/freq_timeint*.frd \
      validation-suite/test-2014/freq_timestep*.frd \
      validation-suite/test-2014/sigxx4.his \
      validation-suite/test-2014/tn30.vtk validation-suite/test-2014/tn31.vtk \
      validation-suite/test-2014/tn32.vtk validation-suite/test-2014/tn33.vtk \
      validation-suite/test-2014/tn34.vtk validation-suite/test-2014/tn35.vtk \
      validation-suite/test-2014/tn36.vtk validation-suite/test-2014/tn40.vtk \
      validation-suite/test-2014/tn41.vtk validation-suite/test-2014/tn42.vtk \
      validation-suite/test-2014/tn43.vtk \
      validation-suite/test-2014/velx.* validation-suite/test-2014/vely.* \
      validation-suite/test-2014/disx.* validation-suite/test-2014/disy.* \
      validation-suite/test-2014/adis.* \
      validation-suite/test-2014/avel.* \
      validation-suite/test-2014/node_dof_0.11 validation-suite/test-2014/node_dof_1.11 \
      validation-suite/test-2014/beam_force_moment.* \
      validation-suite/test-2014/materi_stress_force.* \
      validation-suite/test-2014/element_truss_force_0.7
# 13 tests: 12 preexistentes + familia iface_mc (1 test logico = 6 runs:
# iface_mc a/a' invarianza, iface_mc_slip b/b' invarianza, tension c, gap d)
# + iface_mc_mem (memory), + iface_mc_dil/dil_1step/num (dilatancia RF-4 y
# clamp numerico del MC) + familia 3D (iface_3d friccion 3D, iface_3d_slip
# deslizamiento libre, iface_conv3d conversion tria3->prism6,
# iface_3d_stress print 3D, iface_3d_order orden -x)
# + iface_gen/iface_gen_geom/iface_gen_geom_off (control_mesh_generate_interface)
# + iface_gen_method/iface_gen_method_gen (control_mesh_generate_interface_method)
# + iface_gen_quad9 (subdivision de caras cuadraticas 2D)
# + materi_direct/materi_direct_mc/materi_direct_auto (plasti directo normal;
#   materi_direct_auto is an EXPECTED_ERROR test since 2026-10-01)
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
# + groundflow_pressure_atm/_def (groundflow_pressure_atmospheric, keyword
#   heredado del GNU: cap de la presion estatica en phreatic_coord. Con
#   pa=0.5 la succion +1 se mantiene en 0.5 y la compresion -1 pasa intacta;
#   con default 0 la succion se anula -> 0).
# + groundflow_total_pressure_limit/_dry (manual Professional 6.587; F4 T1
#   retarget to the Professional 02-08-2026: the pres dof is never cut and
#   the consolidation term is kept -> free-node pres 1.6667 with limit 0.5
#   and 0; gvely 0.16667 with 0.5, 0 with 0 = limit reached, dry).
# + contact_block/ctrl_apply/heatgen (remate P6 contacto: bloque en caida
#   detenido por el penalty disy~-0.003 vs -0.055 libre; control_contact_apply
#   0 -no = caida libre exacta; contact_heat_generation nombre Professional
#   del legacy contact_heatgeneration, valor fluye: friction_energy escala
#   0.5->1.312, 1.0->2.624 via TOCHNOG_DEBUG).
# + cd_method/cd_geom (change_dataitem_time_method -tangent: phi almacenado
#   = atan(tabla) verificado directo en el record del grupo con target_item
#   generico; change_dataitem_geometry: split del grupo — clon con los
#   records group_*, elementos dentro del geometry_brick movidos al clon,
#   cambio aplicado solo al clon: bloque izquierdo sigxy~0, derecho ~1 y el
#   record original conserva c=1.0).
# + cda_arith/copy/activate (control_data_arithmetic -multiply/-plus(-all)
#   con _double sobre young: 1000*2+500=2500 target exacto en el record;
#   control_data_copy_index+_factor (g0->g1 x2) y copy+_factor (todos los
#   indices x0.5): 500/1000 targets exactos; control_data_activate -no
#   borra los records bounda_force -> columna descargada relaja a disy~0
#   en cuasiestatico, vs ~0.01 con la fuerza activa).
# + cdist_normal/corr/clamp (control_distribute layout Professional
#   [-normal/-lognormal item index|-all number] + _parameters mean/std +
#   _seed reproducible: deltas por elemento calibrados exactos +35.41/
#   -21.37, distintos entre elementos; _correlation_length 2e12 = campo
#   constante: ambos deltas IGUALES 35.41 pese a 25 unidades de distancia;
#   _minimum_maximum con std 1e6: deltas exactamente +-5 = clamp activo.
#   Layout GNU (triplets + _values, ho_othr4) intacto).
# + condif_heat_edge/vol/vol2 (condif_heat_edge_normal: analogo termico
#   del flux_edge via area() MTYPES 7 — q=0.1 por el borde inferior con
#   T=0 arriba: T fondo 2.0 y medio 1.0 exactos (Fourier); heat_volume con
#   _element: solo elemento 2 calentado -> T(centro)=0.25 analitico (0.5
#   sin restriccion, A/B documentado); _factor polinomio x: S(x)=x ->
#   T(centro)=0.5 analitico).
# + condif_convec/rad/convec_el (condif_convection/radiation_edge_normal,
#   nombres Professional de los legacy condif_convection/radiation: mismos
#   masters tipo 7/8 en area() con fisica compartida h*(Tenv-T) y
#   alfa*(Tr^4-T^4) + variantes _element/_element_group/_element_side/
#   _node/_element_node. Conveccion analitica T0=0.5 (h=k=L=1); radiacion
#   Newton 10 it: T0+T0^4=1 -> 0.7245; _element que no toca la geometria:
#   T=0 exacto, A/B de la restriccion).
# + aeg_node/aeg_seq/bt_factor (area_element_group_node + _method -any +
#   _time -yes: brick solo columna izquierda -> group 1 (young 2000):
#   sigyy -20 vs -10 A/B fisico; sequence con alias Professional
#   _element_group (guiones) + _geometry_method -any: cambio de grupo en
#   t=0.1 duplica sigyy; bounda_time_factor 2.0 sobre tabla 5.0 -> carga
#   efectiva 10, sigyy -10 vs -5).
# + iface_condif/expansion/tangref (group_interface_condif_conductivity:
#   espejo termico del groundflow_interface, bloque aislado se llena a
#   T=2 solo por la interfaz; _materi_expansion_normal: eps=alfa*T con
#   pseudo-carga INCREMENTAL (patron stress.cc), bloques fijos -> sigxx
#   1.98=2*kn*alfa*T (0 sin el record); _tangential_reference_point +z:
#   t1 rota a z, la cortante en y pasa a f_t2=10.4 y f_t~0 (frame por
#   defecto daria f_t=la cortante)).
# + node_force_inertia/slide/pressure (node_force discreto en dof.cc +
#   node_inertia calculado m*a+m*g PUT en VERSION_NEW: -20.11 ~= m*g=-20
#   (target -19.83 until C28 2026-10-06, calibrated on the reversed sign)
#   con node_force -30 actuando; node_pressure overrides de calcul -static:
#   node_static_pressure 7.5 -> post_calcul -static da 7.5 exacto vs -2
#   calculado; node_slide smoke clone de slide_axi con la geometria
#   desplazada y membresia por record — HALLAZGO: slide_axi original
#   tambien da rhside 0, su target ±0.05 pasaba trivial, verificacion
#   fisica del slide pendiente de modelo dedicado).
# + creset_geom/iface (control_reset_geometry/_node/_element_group:
#   filtro de nodos del reset (elementos completamente dentro o todos sus
#   nodos listados) — hisv0 0.5->0 solo columna izquierda, la derecha
#   conserva 0.5, A/B sin filtro falla; control_reset_interface_strain:
#   history strain_normal de la interfaz a 0.0 EXACTO tras el reset).
# + fedge_alias/fedge_restrict/fvol_elem/cmat_gate (Sprint 9: traducción
#   de prefijo force_edge_*->force_element_edge_* y force_volume_*->
#   force_element_volume_* DENTRO de db_number — clave: el detector de
#   fin-de-valores también usa db_number, traducir solo en el punto del
#   keyword rompia el parseo; variantes _element/_element_group/
#   _element_node/_element_side/_node(+_factor) para las 3 familias edge
#   + _element/_element_group para volume; fedge_alias con sintaxis
#   Professional pura (sigxx 5.0), fedge_restrict (_element que no toca
#   la geometria: disx 0), fvol_elem (_element_group: sigxx 0.5 derecha
#   vs 0 izquierda), cmat_gate (plasti_tension_apply -no: cutoff
#   ignorado sigxy -88.4 lineal vs ~1 capped, A/B sin gate falla)).
# + fproj_tunnel (force_edge_projected Terzaghi, sintaxis Professional:
#   campo lineal ph/pv proyectado por pared — sigxx=10 EXACTO (ph en
#   pared vertical) y sigyy=20 (pv en horizontal): el ratio 2:1 prueba
#   la proyeccion por orientacion; una presion uniforme daria 1:1).
# + dsmall/dignore (Sprint 9 lote 3: data_activate+_time borra la carga
#   en t=0.1 -> columna relaja a disy~0; data_ignore -bounda_time ANTES
#   del record -> load default 0, disy 0 exacto; control_solver y
#   axisymmetric aliases via db_number; print_mesh_dof dump smoke).
# + mdirect_comp/gate (Sprint 10: compression_direct espectral — autoval
#   < sigy recortados via matrix_jacobi; columna comprimida sigyy -10
#   capped a -5.0 EXACTO; _visco relaja 1-exp(-dt/tm), sin record corte
#   total. pressure_limit 0.1: la plasticidad directa se desactiva cuando
#   p=-tr/3 excede el limite -> sigyy -10 elástico, A/B sin limit falla).
# + mdp_shear/mfactor (Sprint 10 lote 2: VALIDACION ANALITICA de
#   Drucker-Prager — el GNU implementa f=sqrt(J2)+3*alpha*sm-K (forma
#   estandar DP; NO sqrt(3J2)); cizalla pura con phi=0 (alpha=0, sin
#   presion/dilatancia): sigma_xy cap = K = 2c/sqrt(3); c=30, dt=0.05:
#   34.636 vs 34.641 analitico = 0.014% DE ERROR. El nombre Professional
#   druck_prag via alias db_number del legacy druckprag (misma fisica,
#   tests legacy druckpr1/examp23). group_materi_factor 0.1: escala
#   rigidez+tensiones -> con vely prescrita sigyy=0.1 vs 1.0, A/B).
# + mmc_tension (Sprint 10 lote 3: Mohr-Coulomb CLASICO implementado de
#   cero — bloques plasti.cc con matrix_eigenvalues, f=0.5(s1-s3)+
#   0.5(s1+s3)sin(phi)-c*cos(phi) exacta del manual; VALIDACION via la
#   resistencia a traccion uniaxial clasica sig_t=2c*cos(phi)/(1+sin(phi)):
#   c=40, phi=30 -> 46.19 EXACTO a la primera. El gradiente del flujo lo
#   provee el driver por diferencias finitas. Sanity phi=0 (Tresca):
#   sig_t=2c=80, converge lento por el vertice singular (79.5, doc). El
#   rig de cizalla pura acaba en traccion uniaxial a 45 grados — el mismo
#   sig_t — hallazgo documentado).
# + mmchs_soft (Sprint 10 lote 4: hardening-SOFTENING analitico — c baja
#   linealmente 80->20 con kappa/kappa_crit en [0,1]; phi=0 -> sig_t=2c.
#   CALIBRACION del test (el 1-paso brutal original divergia a 7.55): (1)
#   pasos pequenos dt=0.02 (20 pasos, total 0.4) para acumular kappa
#   gradual; (2) kappa_crit=0.5 mantiene kappa/kappa_crit en zona lineal
#   (0.318/0.5=0.636 -> c=41.8 -> sig_t=83.6); (3) control_timestep_iterations
#   8: el retorno plastico debe perseguir la superficie que BAJA — con el
#   default de 1 iteracion queda por encima (acoplamiento debil). Medido
#   sigxx=81.9 (98% del analitico), target 83.6±5. A/B sin softening
#   (c_1=c_0=80) -> 159.2=2c_0: discrimina. GOTCHA: el target del dof usa
#   el basename -kap, NO -materi_plasti_kappa (ese nombre da lectura fuera
#   de limites/DBL_MAX).
# + mcap2/mcap_legacy (Sprint 10 lote 5: alias Professional
#   group_materi_plasti_cap2 del GNU group_materi_plasti_cap via db_number
#   — misma fisica c phi alpha R + tabla epsilonp_v pb. Oedometro elastico
#   con c=1e6 (f_yield = p - pb < 0 en el camino isotropico): sigxx
#   -0.5769 EXACTO. El gemelo con la keyword legacy da el MISMO valor:
#   cap2 == legacy).
# + mcrunch/mcrunch_low (group_materi_failure_crunching, ortografia
#   Professional con la 'n'; el typo GNU 'cruching' queda como alias
#   db_number. tmp = autovalor principal mas compresivo (negativo); falla
#   si tmp > threshold. Threshold alto 0.5: intacto sigxx -0.5769 EXACTO.
#   Threshold -0.5: elemento marcado para borrado en el PRIMER timestep
#   (element_delete_times[0]=0.05; GOTCHA: threshold negativo dispara
#   incluso a deformacion cero — semantica GNU fiel a jan-2014).
# + mvoid/mvoid_low (group_materi_failure_void_fraction, ortografia
#   Professional con underscore; el GNU 'voidfraction' queda como alias.
#   void inicial 0.3 via node_dof -from/-to; |void| > threshold. 0.9:
#   intacto; 0.2: elemento marcado en t=0.05).
# + mpower (group_materi_elasti_poisson_power, manual 6.653: nu =
#   nu0+nu1*(p/p1)^alpha, nu<=nu2, p=-sig_mean. Oedometro con nu0=0.2
#   nu1=0.1 p1=1 alpha=1 y E*eps=1.2: punto fijo analitico nu=0.4, p=2,
#   sigxx=-1.7143/sigyy=-2.5714. K06 (2026-10-06): the block REPLACES the
#   poisson ratio (before it added C(1000,0.3) + C(1000,nu(p)) ->
#   -1.4389/-3.2924); the incremental update (C from the previous step's
#   stress) gives sigxx -0.4389 / sigyy -1.4101 (Pro 02-08-2026
#   -0.4417/-1.4123); window 0.02, which excludes the constant nu=0.3
#   base (-0.692/-1.615). control_timestep_iterations 8).
# + mshf/mshf_nof (group_materi_elasti_shear_factor, manual 6.654:
#   multiplica la rigidez cortante de young+poisson. Cizalla pura 1 paso
#   con factor 2: sigxy 50 vs 33.33 sin factor; factor 0 -> sigxy 0
#   EXACTO (toda la tension cortante fluye por las entradas escaladas).
#   GOTCHA: el dof de tension post_point es un unknown del sistema
#   acoplado (formulacion de tensiones) — su respuesta a un cambio de
#   tangente NO es lineal: ratio 1.5, no 2; eptxy tambien se desplaza
#   0.0322 -> 0.025 con factor 2).
# + mk0/mk0_off (group_materi_elasti_k0 + control_materi_elasti_k0 -yes,
#   manual 6.650: nu = K0/(1+K0), K0>0.95 truncado. Con K0=0.5 -> nu=1/3:
#   oedometro confinado sigma_xx/sigma_yy = 0.5 EXACTO (sigxx -0.75,
#   sigyy -1.5). A/B con el control -no: nu=0.3 -> -0.5769/-1.3462).
# + myoung6/myoung6_e2/myoung6_e3/myoung6_apply (Sprint 10 lote 6:
#   group_materi_elasti_young_power UPGRADE 3->6 params, manual 6.662 /
#   theory 2.2.2: E = E0 + E1*(p/p1)^alpha con E>=E2 y E<=E3, p=-sig_mean
#   (compresion positiva). Oedometro E0=1000 E1=500 E2=800 E3=3000 p1=1
#   alpha=1: punto fijo analitico E=1714.29 (sigzz -2.3077/sigxx -0.9890);
#   el codigo incremental aterriza en E_eff~1237 = 72% del fijo (medido
#   sigxx -0.713886/sigyy -1.66573; GOTCHA: el dof de tension es unknown
#   del sistema acoplado — la respuesta NO es el fijo secante; el stress
#   acumulado usa el PROMEDIO de los E pasados). CAMBIO SEMANTICO: la
#   forma GNU 3-params (young0*|p/p0|^alpha) ya no se acepta. GOTCHA
#   mayor: C_matrix ACUMULA (array_add), el GNU 2014 con young+young_power
#   sumaba C(1000)+C(E_power) -> rigidez doblada; el bloque young_power
#   ahora LIMPIA C/Cmem antes de construir (semantica Professional: la
#   ley ES el modulo). Clamps A/B analiticos EXACTOS: myoung6_e2 (E2=3000
#   fuerza E=3000: sigxx -1.7308/sigyy -4.0385) y myoung6_e3 (E3=800
#   fuerza E=800: sigxx -0.4615/sigyy -1.0769). myoung6_apply:
#   materi_elasti_young_power_apply -no (manual 6.801) -> E=E0=1000
#   constante EXACTO (sigxx -0.5769/sigyy -1.3462, la base lineal).
# + msph/msph_flat (Sprint 10 lote 6: group_materi_elasti_stress_pressure_history_factor
#   manual 6.655 + initia materi_stress_pressure_history 4.50: el maximo
#   |p| historico se guarda en el dof sph; si la presion actual es MENOR
#   (descarga/recarga) la rigidez se multiplica por factor, si es el nuevo
#   maximo NO. Oedometro 2 fases: carga 4 pasos vely=-0.002 (pico sigyy
#   -0.5385, sph=0.3333) + descarga 2 pasos vely=+0.002. Con factor 3 el
#   tangente de descarga es 3E y la tension SOBREPASA a +0.2692 vs -0.2692
#   con factor 1 (A/B msph_flat). La decision usa la presion ESTIMADA del
#   paso actual (p_old + dp con dp=-mean(C:inc_ept)) contra el maximo
#   historico a INICIO de paso (old_unknowns[sph]): captura el PRIMER paso
#   de descarga y evita la degeneracion p==sph en el pico (redondeo ->
#   factor espurio -> historia runaway). GOTCHA IEEE: scalar_dabs(-0.0)
#   devuelve -0.0 y -0.0 < sph es TRUE -> el factor se aplicaria desde el
#   primer paso de carga; normalizado con p==0. -> p=0. GOTCHA matrix_a4b:
#   NO es in-place seguro (corrompe el buffer de salida).
# + mcap1/mcap1_elast/mcap1_comb (Sprint 10 lote 7: group_materi_plasti_cap1
#   implementado DE CERO — el cap elastoplastico del manual 6.691/teoria
#   cap1: f = q^2/M^2 + p*(p* - p*c), p* = p + c*cot(phi), p*c = pc +
#   c*cot(phi), con pc = dof nodal de la initia materi_plasti_cap1_history
#   (4.17). El ENDURECIMIENTO vive en stress.cc (patron de kappa): pc crece
#   con la deformacion volumetrica plastica del cap, pc_dot =
#   deps_p_cv*K_ref/(lambda*/kappa* - 1)*((pc + c*cot(phi))/p_ref)^m.
#   VALIDACION ANALITICA (hex8, compresion isotropica v=-0.04, K=833.3,
#   pc_0=100, m=0, lambda*/kappa*=10): el cap activa en p=pc=100 (paso 20)
#   y el retorno al punto fijo discreto da dpc = 5.0/9 = 0.5/paso ->
#   pc = 110.0 EXACTO tras 20 pasos plasticos; descarga 5 pasos elasticos
#   -> p = 85.0 -> sigxx -84.9996 (error 0.0005%). El gemelo elastico
#   (mcap1_elast) termina en sigxx -175.0000 EXACTO: discrimina superficie
#   Y endurecimiento. mcap1_comb anade druck_prag (phi 30 c 10): en el
#   camino isotropico f_dp < 0 siempre, el cap1 domina (mecanismo de
#   maxima f de plasti_rule) y la respuesta es IDENTICA a mcap1.
#   GOTCHA bounda_time_increment: la variable local persiste entre records
#   bounda (GET_IF_EXISTS no la resetea): records pairs tras uno con
#   increment se leen como load-only (sus tiempos se vuelven cargas ->
#   velx=100 en un paso). Solucion: formato increment en TODOS los records.
#   GOTCHA indice: k=floor(t/increment) -> 41 cargas + 5 descargas dan
#   40 pasos de carga + 5 de descarga.
# + mhardsoil_elast/mhardsoil_elast2/mhardsoil_unload/mhardsoil_unload_flat
#   /mhardsoil_plast/mhardsoil_plast_elast/mhardsoil_gp0/mhardsoil_gp0_off
#   (Sprint 10 lote 8: HARDENING-SOIL implementado de cero, manual 6.649/
#   6.703/4.22/4.40/6.146. Elastico: E = Eref*((sig3+c*cot(phi))/
#   (sigmaref+c*cot(phi)))^m con sig3_manual = la menor compresiva =
#   mayor autovalor algebraico en codigo (traccion positiva; el menor
#   autovalor daria el axial, NO el confinamiento); primer paso de carga
#   E50/nu50 vs descarga/recarga Eur/nuur con el switch del maximo |p|
#   historico (dof sph compartido con materi_stress_pressure_history);
#   base <= 0 -> E = Eref. Plastico: f = q/(E50*(1-q/qa)) - 2*q/Eur -
#   gamma_p con qa = qf/Rf y qf = 2*sin(phi)*(sig3+c*cot(phi))/
#   (1-sin(phi)) derivado de Mohr-Coulomb en falla (Schanz); gamma_p =
#   el dof kappa (materi_plasti_kappa, int sqrt(0.5*deps_p:deps_p)) +
#   el extra inicial de control_materi_plasti_hardsoil_gammap_initial
#   (6.146: -yes crea gamma_p_extra = f(estado inicial) en el primer
#   paso, guardado en element_intpnt_materi_plasti_hardsoil_gammap_initial
#   y SUMADO a gamma_p en la ley -> f = 0 al arrancar con tensiones
#   desviadoras). Validacion ANALITICA: mhardsoil_elast (uniaxial
#   plano con c=10 phi=30 m=0.5: E50 = 1000*sqrt(17.32/117.32) = 384.23,
#   sigyy = -1.098901*384.23*0.0004 = -0.1689 EXACTO; A/B sigmaref=200
#   -> -0.1241), mhardsoil_unload (oedometro confinado 4+2 pasos, m=0:
#   Eur=3000 -> +0.2692 vs -0.2692 con Eur=1000, patron msph),
#   mhardsoil_plast_elast (gemelo elastico -4.0385 EXACTO),
#   mhardsoil_gp0 (tension inicial sigxx=-2 + control: f=0 -> el estado
#   NO se relaja, sigxx -2.0 EXACTO, record 0.0007763 = f(q=2) analitico,
#   kappa 0) y mhardsoil_gp0_off (sin control: el retorno relaja la
#   desviadora y kappa crece; el estado final cumple f ~ 0).
#   mhardsoil_plast (sigxx=-4). F4 T5 (2026-10-10): tras el fix del
#   kappa del hardsoil (gamma_p = kappa del inicio del paso + incremento
#   del propio retorno) los targets de mhardsoil_plast (kap 0.00104991,
#   sigxx -3.11898) y mhardsoil_gp0_off (kap 0.000494563, sigxx
#   -1.56781) estan calibrados GNU (sin oraculo Pro, decision del owner);
#   los valores viejos 0.0006/-2.44 y 0.000293/-1.225 eran la rama par
#   de una oscilacion de periodo 2 en el numero de iteraciones. GOTCHAS: el
#   retorno plastico del HS (f = O(q/E) ~ 1e-4, minisculo) converge en
#   el punto de integracion pero el dof de tension del sistema acoplado
#   NO lo propaga (el estado medido queda sobre la superficie f ~ 0 con
#   el kappa crecido; el gemelo con carga pura da la respuesta elastica:
#   quirk numerico del formulacion de tensiones con leyes de f pequena,
#   ver manual-developer). node_dof: nder = ndim+2 con derivatives
#   (2D: 4), la lista de valores es UNA por nodo aplicada al rango
#   (nuknwn valores, no nuknwn*nodos). La pre-alocacion del record
#   element_intpnt_* (necesaria: el PUT del primer paso corre en el
#   bucle paralelo de elementos donde no se puede alocar) usa el
#   centinela -1 (top.cc) y el primer paso lo sobrescribe con
#   gamma_p_extra >= 0.
# Sprint 10 lote 9: initias por-modelo materi_strain_plasti_cap/
# _compression/_diprisco/_druckprag (manual 4.35-4.40; dof dedicado
# por modelo, llenado con el MISMO inc_epp que materi_strain_plasti:
# mstrain_<model> plastico con target != 0 en el basename del dof
# (eppcapzz -0.036 exacto en cap1, eppdrpxy identico al epp generico
# 0.910768, eppdipyy -0.0105714) + gemelo elastico con target 0;
# mdiprisco_hist: alias materi_plasti_diprisco_history 11 (4.18) =
# mismo hisv dof, target hisv10 -152.795 identico a diprisc1; camclay:
# group_materi_elasti_camclay_pressure_min (6.647) clamp de la
# presion en K = (1+e)*p/kappa (mc_pressure_min: sigxx -0.0093333 y
# sigyy -0.2093333 analiticos EXACTOS con pressure_min 10 vs
# mc_pressure_min_off sin record: +0.06566, signo invertido =
# K degenerado). diprisco_density: PENDIENTE (ley de interpolacion
# loose/dense no documentada en el manual, ver SEGUIMIENTO).
# Sprint 10 lote 10: mrepeat_save — control_repeat_save (6.348) +
# control_repeat_save_calculate (6.349): columna elastica con vely
# prescrita -0.01, control_repeat 3 saltos a control_timestep dt=0.05:
# cada salto guarda disy del post_point (paso -0.0005, secuencia
# aritmetica exacta) en REPEAT_SAVE_RESULT[indice del repeat]; al
# agotarse el contador se calculan media (-0.001) y varianza
# POBLACIONAL (d^2*(N^2-1)/12 = 1.6666667e-7) en
# REPEAT_CALCULATE_RESULT. Targets exactos sobre ambos records;
# el target de varianza discrimina el NUMERO de saves (un 4o save
# daria 3.125e-7 y media -0.00125). A/B ad-hoc: sin
# control_repeat_save_calculate el record repeat_calculate_result no
# existe (target -> exit 1, verificado por separado).
# Sprint 11 lote 1 (8 keywords control_print_*):
# + dbmeth (control_print_database_method 6.266: -all -> dbmeth20.dbs con
#   TODOS los records base y sin "Size of"; -size_tot -> dbmeth21.dbs con
#   "Size of <record> is <bytes>" por record + "Total size" + "Size of the
#   system matrix is <ecuaciones>"; -size_tot_large -> dbmeth22.dbs con
#   SOLO la linea de la matriz del sistema: todos los records del modelo
#   son < 1 Mb — el A/B entre 21 y 22 discrimina el filtro de tamano).
# + partialname (control_print_partialname 6.337: el stdout contiene los
#   records element* (element, element_group, element_mass, ...) y NINGUN
#   record node*; grep sobre /tmp/partialname_safe.out).
# + elmethod (control_print_element_method 6.286: -middle default -> 2
#   lineas con la coordenada MEDIA del elemento; -node -> 4 lineas con las
#   coordenadas nodales; 2*mid1 == x_node2 EXACTO = la media es el
#   promedio de las coords nodales. La fuerza inicial se relaja a 0 con la
#   deformacion (self-stress), el discriminador es estructural).
# + hreltime (control_print_history_relative_time 6.322: con tr=0.2 la
#   ultima linea de disx1.his es 0.1 (0.3-0.2); sin tr la ultima de
#   disx2.his es 0.4; el A/B de 0.3 discrimina el desplazamiento).
# + numit (control_print_number_iterations 6.336: monitor de consola,
#   control_timestep_iterations 8 -> 16 lineas "control_print_number_iterations:"
#   en el stdout (2 pasos x 8 iteraciones); NO es
#   -inverse_iteration_number).
# + meshdoff (alias control_print_mesh_dof -> print_mesh_dof via db_number:
#   smoke del dump print_mesh_dof.dat, primera linea "1 0 0").
# + dofrhside (alias control_print_dof_rhside -> CONTROL_PRINT_UNKNOWNSRHSIDE
#   via db_number: smoke de velx_rhside.20 con 3 columnas x y rhs).
# + dofid_no (control_print_dof_id 6.270 -no: dof.21 con 3 columnas x y dof;
#   dof1 con el default -yes produce dof.20 con 4 columnas x y dof node —
#   el contraste 4 vs 3 columnas discrimina el default).
# Sprint 11 lote 2 (frecuencia de prints, manual Professional 6.291/6.292):
# + freq_timeint (control_print_frequency_timeinterval 10 0.15 sobre
#   control_timestep 10 0.04 0.41): control_print_dof 10 -separate_sequential
#   escribe SOLO 3 archivos dof.0..dof.2 (t=0.16, 0.32, 0.41 — los tiempos
#   EXACTOS del ejemplo del manual) y control_print_frd 10
#   -separate_sequential uno por print (freq_timeint0-2.frd con la linea
#   100CL = tiempo del print). El bloque 11 sin frecuencia escribe
#   dof.3..dof.12 (10 archivos): el contraste 10 vs 3 discrimina el gate.
# + freq_timestep (control_print_frequency_timestep 22 5): dof.22 se
#   escribe SOLO en t=0.20, 0.40, 0.41 (freq_timestep0-2.frd) y dof.23
#   (sin frecuencia) 11 veces: 3*L23 == 11*L22. control_timestep 0.04 0.40
#   0.01 0.01 (2 incrementos): el GNU CLAMPEA el ultimo paso parcial de un
#   incremento (0.36+0.04 -> 0.41), asi que un solo incremento 0.41 no
#   tendria paso en 0.40; con 2 incrementos la secuencia es 0.04..0.40,0.41.
#   control_print_history NO se gatea (excepcion): sigxx4.his tiene una
#   linea por paso (10 + 11 = 21), no 6.
# Sprint 11 lote 3 (extensiones VTK, manual Professional 6.340-6.346):
# + vtk_coord1 (control_print_vtk_coord: -yes/default escribe el bloque
#   POINTS en tn30.vtk; -no lo omite en tn31.vtk pero CELLS queda).
# + vtk_dofcalc1 (control_print_vtk_dof_calcul: sin filtro -> tn32.vtk con
#   los 2 campos post (materi_strain_total_average + materi_stress_mises);
#   -none -> tn33.vtk sin campos post pero con los primarios; filtro
#   -materi_strain_total -> tn34.vtk solo con el campo ept).
# + vtk_empty1 (control_print_vtk_empty: elemento 2 vacio por densidad 0
#   (ELEMENT_EMPTY=-YES); default -> tn35.vtk con 2 celdas; -no ->
#   tn36.vtk con 1 celda en CELLS y CELL_TYPES).
# + vtk_nodmeth1 (control_print_vtk_node_method: total lagrange con
#   vely=-0.01 -> disy=-0.001; -node -> tn40.vtk con la coordenada
#   almacenada y=1.0; -node_deformed_mesh -> tn41.vtk con y=0.999).
# + vtk_other1 (control_print_vtk_other: default -> tn42.vtk con
#   SCALARS boundary_condition (nodo 1=1.0, nodo 4=0.0: bounda_force no
#   es condicion de contorno) + VECTORS mesh_deformation; -no ->
#   tn43.vtk sin ninguno de los dos campos).
# Sprint 11 lote 4 (dof_line + dof_point, manual Professional
# 6.273-6.283): interpolacion de node_dof Y node_dof_calcul a lo largo
# de una polilinea / en un punto. Modelo base: quad4 elastico 2D con
# vely=-0.01 arriba / 0 abajo -> campo lineal exacto disy=-0.01*t*y (o
# uniforme -0.01*t con vely en TODOS los nodos). node_start_refined se
# da como INPUT (registro NODE-class, version_all=1) -> el metodo por
# defecto -node_start_refined busca en el frame de referencia.
# + dpline1 (interpolacion EXACTA: linea (0.5,0)->(0.5,1) con n=3 ->
#   disy.30 = 3 lineas con 0/-0.005/-0.01; las funciones de forma
#   lineales reproducen un campo lineal).
# + dpline_n (control_print_dof_line_n 6.279): A/B n=3 vs n=5 -> disy.31
#   3 lineas (t=1) vs disy.32 5 lineas (t=2: 0/-0.005/-0.01/-0.015/-0.02;
#   los bloques control_timestep CORREN EN SERIE, el 32 arranca en t=1).
# + dpline_geom (polilinea en V (0,0)->(1,1)->(2,0), n=5): disy.33 con
#   (0.5,0.5), (1,1), (1.5,0.5) -> los puntos se reparten sobre los 2
#   SEGMENTOS (longitud total 2*sqrt(2)), no sobre la recta (0,0)->(2,0).
# + dpline_group (control_print_dof_line_element_group 6.275): elemento
#   1 en group 0, elemento 2 en group 1; con el filtro {1} el punto
#   (0,0.5) (elemento 1) NO se imprime -> disy.34 2 lineas (x=1, x=2) vs
#   disy.35 3 lineas (x=0,1,2) sin filtro.
# + dpline_eps (control_print_dof_line_eps_iso 6.276): punto (1.02,0.5)
#   0.02 fuera de la malla: con el default 1.e-3 NO se acepta (distancia
#   > element_largest_size*eps) -> disy.36 1 linea; con eps_iso=1.0 se
#   acepta (extrapolacion) -> disy.37 2 lineas (x=0.5 y x=1.02).
# + dpline_method (control_print_dof_line_method 6.277): analisis
#   follow-material (materi_velocity SIN materi_displacement) con
#   node_start_refined de input; tras 1 paso dt=1 la malla deformada
#   tiene el borde superior en y=0.99. Linea (0.5,0)->(0.5,0.995), n=2:
#   -node_start_refined -> vely.38 2 lineas (y=0 y y=0.995); -node ->
#   vely.39 1 linea ((0.5,0.995) queda 0.005 POR ENCIMA del borde
#   deformado y no se acepta). GOTCHA: con group_materi_memory -total el
#   check de matrix_inverse(old_rot) estalla (distorsion espuria en el
#   follow-material 1-elemento) -> -updated_without_rotation.
# + dpline_move (control_print_dof_line_move 6.278): vely=-0.01 uniforme
#   (translacion rigida), linea (0.5,0.5)->(0.5,0.6), n=2, 2 pasos
#   dt=0.1: con -yes las coordenadas siguen la particula (y=0.499 en la
#   linea 3 de disy.40) vs fijas (y=0.5 en disy.41).
# + dpline_time (control_print_dof_line_time 6.280) / dpoint_time
#   (6.283): primera linea de cada archivo = "# time 0.1" (comentario
#   gnuplot) en disy.42 / disy.43.
# + dpoint1 (control_print_dof_point 6.281): serie temporal en el punto
#   (0.5,0.5) -> disy.44 2 lineas -0.001/-0.002 (t=0.1/0.2) + calcul
#   node_dof_calcul: post_calcul -materi_displacement -average -> adis.44
#   con la media del vector desplazamiento (disx+disy)/2 = -0.0005/-0.001.
#   GOTCHA del modelo: con fixed_in_space el strain total es 0 (F=I), por
#   eso el calcul usa el desplazamiento, no el strain.
# F0 T4 (2026-10-05): the run list lives in scripts/build_safe_tests.txt
# (one test-2014 input per line, '#' comments). Restored from the full
# pre-incident loop (2488c4a^) minus reset1 (input lost) and
# materi_direct_auto (expected-error test, checked below), plus the
# tsup_* tests and gforce7.
TEST_LIST="$REPO_DIR/scripts/build_safe_tests.txt"
RUN_FAILED=""
while read -r t; do
  HIPO_TOTAL=$((HIPO_TOTAL+1))
  ( cd validation-suite/test-2014 &&
    ulimit -v "$TN_MEMLIMIT_KB" &&
    timeout 120 "$REPO_DIR/build/$TN_BIN" "$t.dat" > "/tmp/${t}_safe.out" 2>&1 )
  RC=$?
  if [ $RC -eq 0 ]; then
    HIPO_OK=$((HIPO_OK+1))
    echo "    $t: OK (exit 0)"
  else
    echo "    $t: FALLO (rc=$RC)"
    RUN_FAILED="$RUN_FAILED $t"
  fi
done < <(sed -e 's/#.*//' -e 's/[[:space:]]//g' "$TEST_LIST" | grep -v '^$')
echo "==> Resumen: $HIPO_OK/$HIPO_TOTAL runs OK (lista: scripts/build_safe_tests.txt)"

# ---------------------------------------------------------------------
# Sprint 11 lote 1: verificacion de ARCHIVOS y STDOUT de los 8 keywords
# control_print_* (los targets de los .dat verifican el modelo; aqui se
# verifican los archivos generados y el stdout capturado en /tmp).
# ---------------------------------------------------------------------
T2014="$REPO_DIR/validation-suite/test-2014"
CHECK_FAIL=0
N_CHECK_OK=0
N_CHECK_FAIL=0
check_ok()   { echo "    $1: OK"; N_CHECK_OK=$((N_CHECK_OK+1)); }
check_fail() { echo "    $1: FALLO ($2)"; CHECK_FAIL=1; N_CHECK_FAIL=$((N_CHECK_FAIL+1)); }

# dof1: control_print_dof_id DEFAULT -yes -> dof.20 con 4 columnas (x y dof node)
if awk 'NF!=4{exit 1}' "$T2014/dof.20" 2>/dev/null; then
  check_ok "dof.20 (control_print_dof_id default -yes: lineas x y dof node)"
else
  check_fail "dof.20" "esperaba 4 columnas por linea (x y dof node)"
fi

# dofid_no: control_print_dof_id -no -> dof.21 con 3 columnas (x y dof).
# C23: 48 lineas = 12 componentes (vel 2, velint 2, dis 2, sig 6) x 4
# nodos, una por componente, y sin basura del heap (denormales e-3xx).
if awk 'NF!=3{exit 1}' "$T2014/dof.21" 2>/dev/null && \
   [ "$(wc -l < "$T2014/dof.21")" = "48" ] && \
   ! rg -q 'e-3[0-9][0-9]' "$T2014/dof.21"; then
  check_ok "dof.21 (control_print_dof_id -no: 48 lineas x y dof, sin basura e-3xx)"
else
  check_fail "dof.21" "esperaba 48 lineas de 3 columnas (x y dof) sin denormales"
fi

# dbmeth: -all sin tamanos; -size_tot con tamanos + matriz del sistema;
# -size_tot_large solo la matriz (todos los records < 1 Mb en el modelo)
if [ -f "$T2014/dbmeth20.dbs" ] && [ "$(grep -c 'Size of' "$T2014/dbmeth20.dbs")" = "0" ] \
   && grep -q "end_data" "$T2014/dbmeth20.dbs" \
   && grep -q "^element  " "$T2014/dbmeth20.dbs"; then
  check_ok "dbmeth20.dbs (method -all: todos los records, sin tamanos)"
else
  check_fail "dbmeth20.dbs" "dump -all incompleto"
fi
if [ -f "$T2014/dbmeth21.dbs" ] && [ "$(grep -c 'Size of' "$T2014/dbmeth21.dbs")" -gt 10 ] \
   && grep -q "Total size is" "$T2014/dbmeth21.dbs" \
   && grep -q "Size of the system matrix is" "$T2014/dbmeth21.dbs"; then
  check_ok "dbmeth21.dbs (method -size_tot: tamanos + total + matriz del sistema)"
else
  check_fail "dbmeth21.dbs" "dump -size_tot incompleto"
fi
if [ -f "$T2014/dbmeth22.dbs" ] && [ "$(grep -c 'Size of' "$T2014/dbmeth22.dbs")" = "1" ] \
   && grep -q "Size of the system matrix is" "$T2014/dbmeth22.dbs"; then
  check_ok "dbmeth22.dbs (method -size_tot_large: solo la matriz del sistema)"
else
  check_fail "dbmeth22.dbs" "dump -size_tot_large incorrecto"
fi

# partialname: el stdout tiene records element* y NINGUN record node*
if [ "$(grep -c '^element' /tmp/partialname_safe.out)" -ge 5 ] \
   && ! grep -q '^node ' /tmp/partialname_safe.out; then
  check_ok "partialname (prefijo -element: solo records element*)"
else
  check_fail "partialname" "el stdout no discrimina el prefijo"
fi

# elmethod: -middle -> 2 lineas (coordenada media); -node -> 4 lineas
# (coordenadas nodales); 2*mid1 == x_node2 (la media es el promedio)
NMID=$(wc -l < "$T2014/element_truss_force_0.20")
NNOD=$(wc -l < "$T2014/element_truss_force_0.21")
MID1=$(awk 'NR==1{print $1}' "$T2014/element_truss_force_0.20")
XN2=$(awk 'NR==2{print $1}' "$T2014/element_truss_force_0.21")
if [ "$NMID" = "2" ] && [ "$NNOD" = "4" ] && \
   awk -v a="$MID1" -v b="$XN2" 'BEGIN{ d=2*a-b; exit !(d<1.e-3 && d>-1.e-3) }'; then
  check_ok "elmethod (-middle 2 lineas vs -node 4 lineas; 2*mid1=$MID1 == x_node2=$XN2)"
else
  check_fail "elmethod" "lineas o relacion media/nodal incorrecta"
fi

# hreltime: con tr=0.2 la ultima linea de disx1.his es 0.1 (0.3-0.2);
# sin tr la ultima de disx2.his es 0.4
T1=$(awk 'END{print $1}' "$T2014/disx1.his")
T2=$(awk 'END{print $1}' "$T2014/disx2.his")
if awk -v a="$T1" -v b="$T2" \
   'BEGIN{ d1=a-0.1; d2=b-0.4; exit !(d1<1.e-3 && d1>-1.e-3 && d2<1.e-3 && d2>-1.e-3) }'; then
  check_ok "hreltime (tr=0.2: ultima linea 0.1 vs 0.4 sin tr)"
else
  check_fail "hreltime" "tiempos relativos incorrectos (t1=$T1 t2=$T2)"
fi

# numit: 16 lineas de monitoreo en stdout (2 pasos x control_timestep_iterations 8)
NIT=$(grep -c "control_print_number_iterations:" /tmp/numit_safe.out)
if [ "$NIT" = "16" ]; then
  check_ok "numit (16 monitoreos de iteracion: 2 pasos x 8 iteraciones)"
else
  check_fail "numit" "esperaba 16 monitoreos, hay $NIT"
fi

# meshdoff: alias control_print_mesh_dof -> dump print_mesh_dof.dat (smoke)
if [ -f "$T2014/print_mesh_dof.dat" ] && [ "$(head -1 "$T2014/print_mesh_dof.dat")" = "1 0 0" ]; then
  check_ok "meshdoff (alias control_print_mesh_dof: print_mesh_dof.dat generado)"
else
  check_fail "meshdoff" "print_mesh_dof.dat no generado o primera linea incorrecta"
fi

# dofrhside: alias control_print_dof_rhside -> velx_rhside.20 (x y rhs)
if [ -f "$T2014/velx_rhside.20" ] && awk 'NF!=3{exit 1}' "$T2014/velx_rhside.20"; then
  check_ok "dofrhside (alias control_print_dof_rhside: velx_rhside.20 x y rhs)"
else
  check_fail "dofrhside" "velx_rhside.20 no generado o formato incorrecto"
fi

# ---------------------------------------------------------------------
# Sprint 11 lote 2: control_print_frequency_timeinterval (6.291) y
# control_print_frequency_timestep (6.292). Los .frd de
# -separate_sequential son UNO POR PRINT con la linea 100CL = tiempo
# exacto del print (campo 3); los conteos de dof.* discriminan el gate.
# ---------------------------------------------------------------------

# freq_timeint: 3 prints SOLO en t=0.16, 0.32, 0.41 (ejemplo del manual
# 6.291: control_timestep 10 0.04 0.41 + frequency_timeinterval 10 0.15)
FTI0=$(grep "100CL" "$T2014/freq_timeint0.frd" 2>/dev/null | awk '{print $3}')
FTI1=$(grep "100CL" "$T2014/freq_timeint1.frd" 2>/dev/null | awk '{print $3}')
FTI2=$(grep "100CL" "$T2014/freq_timeint2.frd" 2>/dev/null | awk '{print $3}')
NFRD_TI=$(ls "$T2014"/freq_timeint?.frd 2>/dev/null | wc -l)
if [ "$NFRD_TI" = "3" ] && \
   awk -v a="$FTI0" -v b="$FTI1" -v c="$FTI2" \
   'BEGIN{ d1=a-0.16; d2=b-0.32; d3=c-0.41;
          exit !(d1<1.e-6 && d1>-1.e-6 && d2<1.e-6 && d2>-1.e-6 && d3<1.e-6 && d3>-1.e-6) }'; then
  check_ok "freq_timeint (3 prints SOLO en t=0.16, 0.32, 0.41: 100CL $FTI0/$FTI1/$FTI2)"
else
  check_fail "freq_timeint" "esperaba 3 prints en 0.16/0.32/0.41, hay $NFRD_TI frd con tiempos $FTI0/$FTI1/$FTI2"
fi

# freq_timeint: dof.0..dof.2 = 3 archivos de 40 lineas (3 prints gated;
# 10 componentes velx vely velintx velinty + 6 sig x 4 nodos - eran 176
# hasta el fix C23 del bucle de print_dof, que repetia bloques y leia
# node_dof fuera del record);
# el bloque 11 sin frecuencia escribe dof.3..dof.12 (10 archivos) ->
# 13 archivos en total. Sin el gate habria 13 archivos solo de dof.0..12
# con los MISMO conteos... el discriminator real es 3 vs 10 por bloque.
NDF_TI=$(ls "$T2014"/dof.[0-9] "$T2014"/dof.1[0-2] 2>/dev/null | wc -l)
L0_TI=$(wc -l < "$T2014/dof.0" 2>/dev/null || echo 0)
L12_TI=$(wc -l < "$T2014/dof.12" 2>/dev/null || echo 0)
if [ "$NDF_TI" = "13" ] && [ "$L0_TI" = "40" ] && [ "$L12_TI" = "40" ]; then
  check_ok "freq_timeint (dof.0..dof.12 = $NDF_TI archivos: 3 gated + 10 sin frecuencia, 40 lineas/print)"
else
  check_fail "freq_timeint dof.*" "esperaba 13 archivos de 40 lineas (3 gated + 10), hay $NDF_TI archivos (dof.0=$L0_TI dof.12=$L12_TI lineas)"
fi

# freq_timestep: 3 prints SOLO en t=0.20, 0.40, 0.41 (ejemplo del manual
# 6.292: control_print_frequency_timestep 22 5)
FTS0=$(grep "100CL" "$T2014/freq_timestep0.frd" 2>/dev/null | awk '{print $3}')
FTS1=$(grep "100CL" "$T2014/freq_timestep1.frd" 2>/dev/null | awk '{print $3}')
FTS2=$(grep "100CL" "$T2014/freq_timestep2.frd" 2>/dev/null | awk '{print $3}')
NFRD_TS=$(ls "$T2014"/freq_timestep?.frd 2>/dev/null | wc -l)
if [ "$NFRD_TS" = "3" ] && \
   awk -v a="$FTS0" -v b="$FTS1" -v c="$FTS2" \
   'BEGIN{ d1=a-0.20; d2=b-0.40; d3=c-0.41;
          exit !(d1<1.e-6 && d1>-1.e-6 && d2<1.e-6 && d2>-1.e-6 && d3<1.e-6 && d3>-1.e-6) }'; then
  check_ok "freq_timestep (3 prints SOLO en t=0.20, 0.40, 0.41: 100CL $FTS0/$FTS1/$FTS2)"
else
  check_fail "freq_timestep" "esperaba 3 prints en 0.20/0.40/0.41, hay $NFRD_TS frd con tiempos $FTS0/$FTS1/$FTS2"
fi

# freq_timestep: dof.22 (3 prints gated) vs dof.23 (11 prints sin
# frecuencia): mismo numero de lineas por print -> 3*L23 == 11*L22
L22_TS=$(wc -l < "$T2014/dof.22" 2>/dev/null || echo 0)
L23_TS=$(wc -l < "$T2014/dof.23" 2>/dev/null || echo 0)
if [ "$L22_TS" != "0" ] && [ "$((3*L23_TS))" = "$((11*L22_TS))" ]; then
  check_ok "freq_timestep (ratio 3*L23==11*L22: dof.22=$L22_TS lineas, dof.23=$L23_TS)"
else
  check_fail "freq_timestep ratio" "3*$L23_TS != 11*$L22_TS"
fi

# excepciones: control_print_history NO se gatea -> sigxx4.his tiene una
# linea por paso (10 pasos de freq_timeint + 11 de freq_timestep = 21),
# no 6 (2 bloques x 3 prints gated)
NHIS=$(wc -l < "$T2014/sigxx4.his" 2>/dev/null || echo 0)
if [ "$NHIS" = "21" ]; then
  check_ok "excepcion history (sigxx4.his con $NHIS lineas = 10+11 pasos: control_print_history NO gateado)"
else
  check_fail "excepcion history" "esperaba 21 lineas en sigxx4.his (10+11 pasos), hay $NHIS"
fi

# ---------------------------------------------------------------------
# Sprint 11 lote 3: verificacion de ARCHIVOS de las 5 extensiones VTK
# (control_print_vtk_coord 6.340, _dof_calcul 6.342, _empty 6.343,
# _node_method 6.345, _other 6.346). Los .vtk se regeneran por paso y
# quedan los de la ultima escritura.
# ---------------------------------------------------------------------

# vtk_coord1: default (o -yes) escribe el bloque POINTS; -no lo omite
# (tn31.vtk conserva CELLS/CELL_TYPES: solo desaparecen las coordenadas)
if [ "$(grep -c '^POINTS' "$T2014/tn30.vtk" 2>/dev/null)" -ge 1 ] && \
   [ "$(grep -c '^POINTS' "$T2014/tn31.vtk" 2>/dev/null)" = "0" ] && \
   grep -q '^CELLS' "$T2014/tn31.vtk"; then
  check_ok "vtk_coord (tn30 con POINTS por defecto; tn31 -no sin POINTS pero con CELLS)"
else
  check_fail "vtk_coord" "el bloque POINTS no discrimina -yes/-no"
fi

# vtk_dofcalc1: sin filtro -> los 2 campos post; -none -> ninguno (los
# primarios intactos: VECTORS materi_velocity); filtro -materi_strain_total
# -> solo el campo ept (materi_stress_mises ausente)
if grep -q "materi_strain_total_average" "$T2014/tn32.vtk" 2>/dev/null && \
   grep -q "materi_stress_mises" "$T2014/tn32.vtk" && \
   ! grep -q "materi_strain_total_average" "$T2014/tn33.vtk" && \
   ! grep -q "materi_stress_mises" "$T2014/tn33.vtk" && \
   grep -q "materi_velocity" "$T2014/tn33.vtk" && \
   grep -q "materi_strain_total_average" "$T2014/tn34.vtk" && \
   ! grep -q "materi_stress_mises" "$T2014/tn34.vtk"; then
  check_ok "vtk_dof_calcul (sin filtro ambos campos post; -none ninguno; filtro solo ept; primarios intactos)"
else
  check_fail "vtk_dof_calcul" "el filtro de campos post no discrimina"
fi

# vtk_empty1: default -> 2 celdas (elemento vacio incluido); -no -> 1
# celda (el elemento vacio por densidad 0 se excluye de CELLS/CELL_TYPES)
C35=$(awk '/^CELLS /{print $2; exit}' "$T2014/tn35.vtk" 2>/dev/null)
C36=$(awk '/^CELLS /{print $2; exit}' "$T2014/tn36.vtk" 2>/dev/null)
T36=$(awk '/^CELL_TYPES /{print $2; exit}' "$T2014/tn36.vtk" 2>/dev/null)
if [ "$C35" = "2" ] && [ "$C36" = "1" ] && [ "$T36" = "1" ]; then
  check_ok "vtk_empty (tn35 $C35 celdas incluye el vacio; tn36 -no $C36 celda)"
else
  check_fail "vtk_empty" "el conteo de celdas no discrimina -yes/-no (C35=$C35 C36=$C36 T36=$T36)"
fi

# vtk_nodmeth1: -node -> coordenadas almacenadas (y del nodo 4 = 1.0);
# -node_deformed_mesh -> deformadas (y ~ 1 + disy = 0.998, disy=-0.002
# con vely=-0.01 durante 0.2)
Y40=$(awk '/^POINTS /{f=1;next} f&&/^CELLS /{exit} f&&NF{n++; if(n==4){print $2; exit}}' "$T2014/tn40.vtk" 2>/dev/null)
Y41=$(awk '/^POINTS /{f=1;next} f&&/^CELLS /{exit} f&&NF{n++; if(n==4){print $2; exit}}' "$T2014/tn41.vtk" 2>/dev/null)
if awk -v a="$Y40" -v b="$Y41" \
   'BEGIN{ d1=a-1.0; d2=(b-a)+0.002; exit !(d1<1.e-3 && d1>-1.e-3 && d2<1.e-3 && d2>-1.e-3) }'; then
  check_ok "vtk_node_method (tn40 y=$Y40 almacenada vs tn41 y=$Y41 deformada)"
else
  check_fail "vtk_node_method" "coordenadas no discriminan -node/-node_deformed_mesh (y40=$Y40 y41=$Y41)"
fi

# vtk_other1: default -> boundary_condition (nodo 1 = 1.0 por bounda_unknown,
# nodo 4 = 0.0 por bounda_force que NO es condicion de contorno) +
# mesh_deformation; -no -> ninguno de los dos campos
BC=$(awk '/^SCALARS boundary_condition/{f=1;next} f&&/^$/{exit} f&&/^LOOKUP_TABLE/{next} f' "$T2014/tn42.vtk" 2>/dev/null)
B1=$(echo "$BC" | awk 'NR==1{print $1}')
B4=$(echo "$BC" | awk 'NR==4{print $1}')
if grep -q "boundary_condition" "$T2014/tn42.vtk" 2>/dev/null && \
   grep -q "mesh_deformation" "$T2014/tn42.vtk" && \
   [ "$B1" = "1.0" ] && [ "$B4" = "0.0" ] && \
   ! grep -q "boundary_condition" "$T2014/tn43.vtk" && \
   ! grep -q "mesh_deformation" "$T2014/tn43.vtk"; then
  check_ok "vtk_other (tn42 boundary_condition $B1/$B4 + mesh_deformation; tn43 -no sin campos)"
else
  check_fail "vtk_other" "los campos other no discriminan -yes/-no (B1=$B1 B4=$B4)"
fi

# ---------------------------------------------------------------------
# Sprint 11 lote 4: verificacion de ARCHIVOS de dof_line/dof_point
# (manual Professional 6.273-6.283). Lineas "x y <dof>", un archivo por
# label de dof (disy.30, vely.38, ...) y por item de node_dof_calcul
# (adis.44 = "a"+dis, convencion de print_unknowns).
# ---------------------------------------------------------------------

# dpline1: interpolacion EXACTA del campo lineal disy=-0.01*y en la
# linea (0.5,0)->(0.5,1) con n=3: 3 lineas con 0/-0.005/-0.01
N30=$(wc -l < "$T2014/disy.30" 2>/dev/null || echo 0)
L30_2=$(awk 'NR==2{print $2, $3}' "$T2014/disy.30" 2>/dev/null)
L30_3=$(awk 'NR==3{print $2, $3}' "$T2014/disy.30" 2>/dev/null)
Y2_30=$(echo "$L30_2" | awk '{print $1}'); V2_30=$(echo "$L30_2" | awk '{print $2}')
Y3_30=$(echo "$L30_3" | awk '{print $1}'); V3_30=$(echo "$L30_3" | awk '{print $2}')
if [ "$N30" = "3" ] && awk -v y="$Y2_30" -v v="$V2_30" \
   'BEGIN{d1=y-0.5; d2=v+0.005; exit !(d1<1e-3 && d1>-1e-3 && d2<1e-3 && d2>-1e-3)}' && \
   awk -v y="$Y3_30" -v v="$V3_30" \
   'BEGIN{d1=y-1.0; d2=v+0.01; exit !(d1<1e-3 && d1>-1e-3 && d2<1e-3 && d2>-1e-3)}'; then
  check_ok "dpline1 (interpolacion exacta: (0.5,0.5,-0.005) y (0.5,1,-0.01))"
else
  check_fail "dpline1" "disy.30 no interpola el campo lineal (N=$N30, l2=$L30_2, l3=$L30_3)"
fi

# dpline_n: A/B del numero de puntos -> 3 vs 5 lineas; el punto 2 de n=5
# cae en y=0.25 con -0.005 (t=2) y el ultimo en (0.5,1,-0.02)
N31=$(wc -l < "$T2014/disy.31" 2>/dev/null || echo 0)
N32=$(wc -l < "$T2014/disy.32" 2>/dev/null || echo 0)
L32_5=$(awk 'NR==5{print $2, $3}' "$T2014/disy.32" 2>/dev/null)
Y5_32=$(echo "$L32_5" | awk '{print $1}'); V5_32=$(echo "$L32_5" | awk '{print $2}')
if [ "$N31" = "3" ] && [ "$N32" = "5" ] && \
   awk -v y="$Y5_32" -v v="$V5_32" \
   'BEGIN{d1=y-1.0; d2=v+0.02; exit !(d1<1e-3 && d1>-1e-3 && d2<1e-3 && d2>-1e-3)}'; then
  check_ok "dpline_n (n=3 -> $N31 lineas vs n=5 -> $N32 lineas; ultimo punto (0.5,1,-0.02))"
else
  check_fail "dpline_n" "el numero de puntos no discrimina (N31=$N31 N32=$N32, l5=$L32_5)"
fi

# dpline_geom: polilinea en V -> los puntos 2,3,4 caen en (0.5,0.5),
# (1,1), (1.5,0.5) (sobre los segmentos, no sobre la recta (0,0)-(2,0))
G33_2=$(awk 'NR==2{print $1, $2, $3}' "$T2014/disy.33" 2>/dev/null)
G33_3=$(awk 'NR==3{print $1, $2, $3}' "$T2014/disy.33" 2>/dev/null)
G33_4=$(awk 'NR==4{print $1, $2, $3}' "$T2014/disy.33" 2>/dev/null)
if awk -v a="$G33_2" -v b="$G33_3" -v c="$G33_4" \
   'BEGIN{ split(a,A); split(b,B); split(c,C);
     exit !( (A[1]-0.5<1e-3&&A[1]-0.5>-1e-3&&A[2]-0.5<1e-3&&A[2]-0.5>-1e-3&&A[3]+0.005<1e-3&&A[3]+0.005>-1e-3) &&
             (B[1]-1.0<1e-3&&B[1]-1.0>-1e-3&&B[2]-1.0<1e-3&&B[2]-1.0>-1e-3&&B[3]+0.01<1e-3&&B[3]+0.01>-1e-3) &&
             (C[1]-1.5<1e-3&&C[1]-1.5>-1e-3&&C[2]-0.5<1e-3&&C[2]-0.5>-1e-3&&C[3]+0.005<1e-3&&C[3]+0.005>-1e-3) ) }'; then
  check_ok "dpline_geom (polilinea en V: (0.5,0.5), (1,1), (1.5,0.5))"
else
  check_fail "dpline_geom" "la distribucion sobre la polilinea es incorrecta (l2=$G33_2 l3=$G33_3 l4=$G33_4)"
fi

# dpline_group: con el filtro {1} el punto (0,0.5) (elemento 1, group 0)
# NO se imprime -> 2 lineas empezando en x=1; sin filtro -> 3 lineas
N34=$(wc -l < "$T2014/disy.34" 2>/dev/null || echo 0)
N35=$(wc -l < "$T2014/disy.35" 2>/dev/null || echo 0)
X1_34=$(awk 'NR==1{print $1}' "$T2014/disy.34" 2>/dev/null)
X1_35=$(awk 'NR==1{print $1}' "$T2014/disy.35" 2>/dev/null)
if [ "$N34" = "2" ] && [ "$N35" = "3" ] && [ "$X1_34" = "1" ] && [ "$X1_35" = "0" ]; then
  check_ok "dpline_group (filtro {1}: $N34 lineas desde x=$X1_34 vs sin filtro $N35 desde x=$X1_35)"
else
  check_fail "dpline_group" "el filtro de grupos no discrimina (N34=$N34 N35=$N35 x34=$X1_34 x35=$X1_35)"
fi

# dpline_eps: con el default el punto (1.02,0.5) fuera de la malla NO se
# acepta -> 1 linea; con eps_iso=1.0 se acepta -> 2 lineas (x=1.02)
N36=$(wc -l < "$T2014/disy.36" 2>/dev/null || echo 0)
N37=$(wc -l < "$T2014/disy.37" 2>/dev/null || echo 0)
X2_37=$(awk 'NR==2{print $1}' "$T2014/disy.37" 2>/dev/null)
if [ "$N36" = "1" ] && [ "$N37" = "2" ] && [ "$X2_37" = "1.02" ]; then
  check_ok "dpline_eps (default: $N36 linea; eps_iso=1.0: $N37 lineas con x=$X2_37)"
else
  check_fail "dpline_eps" "eps_iso no discrimina (N36=$N36 N37=$N37 x2=$X2_37)"
fi

# dpline_method: -node_start_refined -> 2 lineas (y=0 y y=0.995);
# -node -> 1 linea (el punto y=0.995 queda sobre el borde deformado 0.99)
N38=$(wc -l < "$T2014/vely.38" 2>/dev/null || echo 0)
N39=$(wc -l < "$T2014/vely.39" 2>/dev/null || echo 0)
Y2_38=$(awk 'NR==2{print $2, $3}' "$T2014/vely.38" 2>/dev/null)
if [ "$N38" = "2" ] && [ "$N39" = "1" ] && \
   awk -v a="$Y2_38" 'BEGIN{split(a,A); d1=A[1]-0.995; d2=A[2]+0.00995;
     exit !(d1<1e-3 && d1>-1e-3 && d2<1e-3 && d2>-1e-3)}'; then
  check_ok "dpline_method (-node_start_refined $N38 lineas vs -node $N39; y2=$Y2_38)"
else
  check_fail "dpline_method" "el metodo no discrimina (N38=$N38 N39=$N39 y2=$Y2_38)"
fi

# dpline_move: con -yes la linea 3 (2o paso) imprime en y=0.499 (la
# particula se movio con vely=-0.01*dt=0.1); sin move y=0.5
Y3_40=$(awk 'NR==3{print $2}' "$T2014/disy.40" 2>/dev/null)
Y3_41=$(awk 'NR==3{print $2}' "$T2014/disy.41" 2>/dev/null)
N40=$(wc -l < "$T2014/disy.40" 2>/dev/null || echo 0)
N41=$(wc -l < "$T2014/disy.41" 2>/dev/null || echo 0)
if [ "$N40" = "4" ] && [ "$N41" = "4" ] && \
   awk -v a="$Y3_40" -v b="$Y3_41" \
   'BEGIN{d1=a-0.499; d2=b-0.5; exit !(d1<1e-3 && d1>-1e-3 && d2<1e-3 && d2>-1e-3)}'; then
  check_ok "dpline_move (linea 3: y=$Y3_40 movida vs y=$Y3_41 fija)"
else
  check_fail "dpline_move" "el move no discrimina (y3_40=$Y3_40 y3_41=$Y3_41)"
fi

# dpline_time / dpoint_time: primera linea = comentario gnuplot "# time 0.1"
if [ "$(head -1 "$T2014/disy.42" 2>/dev/null)" = "# time 0.1" ] && \
   [ "$(head -1 "$T2014/disy.43" 2>/dev/null)" = "# time 0.1" ] && \
   [ "$(wc -l < "$T2014/disy.42" 2>/dev/null || echo 0)" = "3" ] && \
   [ "$(wc -l < "$T2014/disy.43" 2>/dev/null || echo 0)" = "2" ]; then
  check_ok "dpline_time/dpoint_time (primera linea '# time 0.1' + lineas de datos)"
else
  check_fail "dpline_time/dpoint_time" "el comentario de tiempo no se escribe"
fi

# dpoint1: serie temporal en el punto (0.5,0.5): disy = -0.001 (t=0.1) y
# -0.002 (t=0.2); adis.44 (node_dof_calcul: media del desplazamiento)
# = -0.0005 y -0.001
V1_44=$(awk 'NR==1{print $3}' "$T2014/disy.44" 2>/dev/null)
V2_44=$(awk 'NR==2{print $3}' "$T2014/disy.44" 2>/dev/null)
A1_44=$(awk 'NR==1{print $3}' "$T2014/adis.44" 2>/dev/null)
A2_44=$(awk 'NR==2{print $3}' "$T2014/adis.44" 2>/dev/null)
if awk -v a="$V1_44" -v b="$V2_44" -v c="$A1_44" -v d="$A2_44" \
   'BEGIN{exit !(a+0.001<1e-3&&a+0.001>-1e-3 && b+0.002<1e-3&&b+0.002>-1e-3 &&
                 c+0.0005<1e-3&&c+0.0005>-1e-3 && d+0.001<1e-3&&d+0.001>-1e-3)}'; then
  check_ok "dpoint1 (disy.44 $V1_44/$V2_44 serie; adis.44 $A1_44/$A2_44 calcul)"
else
  check_fail "dpoint1" "punto o calcul incorrectos (disy $V1_44/$V2_44, adis $A1_44/$A2_44)"
fi

# ---------------------------------------------------------------------
# Sprint 11 lote 5: verificacion de ARCHIVOS de control_print_node
# (manual Professional 6.330-6.335) y del suavizado de control_print_dof
# (6.271-6.272). control_print_node: "x y <value>" por nodo, un archivo
# por parte seleccionada; el suavizado modifica la salida dof.<index>.
# ---------------------------------------------------------------------

# cpn1: labels -velx -vely (ejemplo del manual) -> vely.10 con el campo
# EXACTO 0/0/-0.01/-0.01 (todos los dofs Dirichlet); partes numericas
# 0 1 -> node_dof_1.11 == vely.10; node_dof_calcul -materi_velocity
# (post_calcul -average) -> avel.12 = (velx+vely)/2 = -0.005 arriba
N10=$(wc -l < "$T2014/vely.10" 2>/dev/null || echo 0)
V10_4=$(awk 'NR==4{print $3}' "$T2014/vely.10" 2>/dev/null)
AV12_4=$(awk 'NR==4{print $3}' "$T2014/avel.12" 2>/dev/null)
if [ "$N10" = "4" ] && \
   awk -v v="$V10_4" 'BEGIN{d=v+0.01; exit !(d<1e-3 && d>-1e-3)}' && \
   awk -v v="$AV12_4" 'BEGIN{d=v+0.005; exit !(d<1e-3 && d>-1e-3)}' && \
   diff -q "$T2014/vely.10" "$T2014/node_dof_1.11" >/dev/null 2>&1; then
  check_ok "cpn1 (labels vely.10; numerico node_dof_1.11 == vely.10; calcul avel.12)"
else
  check_fail "cpn1" "labels/numerico/calcul incorrectos (N10=$N10 v10_4=$V10_4 avel_4=$AV12_4)"
fi

# cpn_angular: angulo en GRADOS desde +x hacia +y con middle (0.5,0.5):
# nodo 4 (1,1) = atan2(0.5,0.5) = 45 grados EXACTO; sin sort el orden de
# nodo es -135/-45/135/45; con sort -angle ascendente -135/-45/45/135
A20_4=$(awk 'NR==4{print $1}' "$T2014/vely.20" 2>/dev/null)
A21_1=$(awk 'NR==1{print $1}' "$T2014/vely.21" 2>/dev/null)
A21_4=$(awk 'NR==4{print $1}' "$T2014/vely.21" 2>/dev/null)
if awk -v a="$A20_4" -v b="$A21_1" -v c="$A21_4" \
   'BEGIN{exit !(a-45<1e-6 && a-45>-1e-6 && b+135<1e-6 && b+135>-1e-6 && c-135<1e-6 && c-135>-1e-6)}'; then
  check_ok "cpn_angular (nodo (1,1)=45 grados; sort -angle: -135..135)"
else
  check_fail "cpn_angular" "angulos incorrectos (a4=$A20_4 s1=$A21_1 s4=$A21_4)"
fi

# cpn_geom: geometry_line sobre el borde inferior -> SOLO los 2 nodos del
# borde (2 lineas) vs 4 lineas sin el record. Icrontrols 54/55: los
# 30-44 estan ocupados por la familia dof_line (dpline1 escribe
# vely.30; dpline_move vely.40/41 - colision de nombres de archivo).
N30=$(wc -l < "$T2014/vely.54" 2>/dev/null || echo 0)
N31=$(wc -l < "$T2014/vely.55" 2>/dev/null || echo 0)
if [ "$N30" = "2" ] && [ "$N31" = "4" ]; then
  check_ok "cpn_geom (filtro geometria: $N30 lineas vs $N31 sin filtro)"
else
  check_fail "cpn_geom" "el filtro de geometria no discrimina (N30=$N30 N31=$N31)"
fi

# cpn_sort: nodos 1=(0,1) 2=(0,0) 3=(1,1) 4=(1,0) -> sin sort la primera
# linea es "0 1 -0.01" (orden de nodo); con sort -y la primera es
# "0 0 0" y la ultima "1 1 -0.01" (y ascendente). Icrontrols 56/57.
S40_1=$(head -1 "$T2014/vely.56" 2>/dev/null)
S41_1=$(head -1 "$T2014/vely.57" 2>/dev/null)
S41_4=$(tail -1 "$T2014/vely.57" 2>/dev/null)
if [ "$S40_1" = "0 1 -0.01" ] && [ "$S41_1" = "0 0 0" ] && [ "$S41_4" = "1 1 -0.01" ]; then
  check_ok "cpn_sort (sort -y ascendente: $S41_1 .. $S41_4 vs $S40_1 sin sort)"
else
  check_fail "cpn_sort" "el sort no discrimina (56_1=$S40_1 57_1=$S41_1 57_4=$S41_4)"
fi

# cpn_zero: -no suprime los valores 0 (2 lineas vs 4 con el default -yes);
# con -velx (TODO 0) el archivo queda VACIO (comparacion de cero EXACTA)
N50=$(wc -l < "$T2014/vely.50" 2>/dev/null || echo 0)
N51=$(wc -l < "$T2014/vely.51" 2>/dev/null || echo 0)
N52=$(wc -l < "$T2014/velx.52" 2>/dev/null || echo 0)
if [ "$N50" = "2" ] && [ "$N51" = "4" ] && [ "$N52" = "0" ]; then
  check_ok "cpn_zero (zero -no: $N50 lineas vs default $N51; velx todo 0: $N52)"
else
  check_fail "cpn_zero" "el filtro de ceros no discrimina (N50=$N50 N51=$N51 N52=$N52)"
fi

# dsmooth1: 1 pasada con -all sobre la cadena 0,1,2,3,4 -> el interior es
# el promedio de vecinos EXACTO (1,1,2,3,3) vs raw (0,1,2,3,4) en dof.61
SEQ60=$(awk '{printf "%s ", $2}' "$T2014/dof.60" 2>/dev/null)
SEQ61=$(awk '{printf "%s ", $2}' "$T2014/dof.61" 2>/dev/null)
if [ "$SEQ60" = "1 1 2 3 3 " ] && [ "$SEQ61" = "0 1 2 3 4 " ]; then
  check_ok "dsmooth1 (1 pasada: interior=promedio de vecinos (1,1,2,3,3) vs raw)"
else
  check_fail "dsmooth1" "suavizado 1 pasada incorrecto (dof.60=$SEQ60 dof.61=$SEQ61)"
fi

# dsmooth_n: A/B 3 pasadas (1.5,1.5,2,2.5,2.5 EXACTO) vs default 10
# pasadas (1.9375,1.96875,2,2.03125,2.0625): todos los valores a menos
# de 0.1 de la MEDIA 2.0 (convergencia) y el nodo central x=2 EXACTO
SEQ70=$(awk '{printf "%s ", $2}' "$T2014/dof.70" 2>/dev/null)
D71_3=$(awk 'NR==3{print $2}' "$T2014/dof.71" 2>/dev/null)
if [ "$SEQ70" = "1.5 1.5 2 2.5 2.5 " ] && \
   awk '{d=$2-2; if (d<0) d=-d; if (d>0.1) exit 1}' "$T2014/dof.71" 2>/dev/null \
   && [ "$D71_3" = "2" ]; then
  check_ok "dsmooth_n (3 pasadas exactas; default 10 convergido a la media 2.0)"
else
  check_fail "dsmooth_n" "pasadas incorrectas (dof.70=$SEQ70 d71_3=$D71_3)"
fi

# ---------------------------------------------------------------------
# Sprint 11 lote 6: control_print_beam_force_moment (6.262-6.264)
# ---------------------------------------------------------------------
# bmom1: 1 linea, dist = 0.5 EXACTO (corte vertical en x=0 sobre el nodo
# fijo, fixed_in_space), fy1 = -F = -0.01, mz1 = -F*L = -0.01, fy2 =
# +F = +0.01, mz2 = 0 (snap del ruido del solver), axial 0 (viga pura
# sin truss). Signo = el del vector de fuerzas internas del elemento
# (ELEMENT_BEAM_MOMENT), verificado empiricamente. Columnas: dist fx1
# fy1 fz1 mx1 my1 mz1 fx2 fy2 fz2 mx2 my2 mz2.
if [ "$(wc -l < "$T2014/beam_force_moment.5" 2>/dev/null)" = "1" ] && \
   awk 'NF==13 && $1>0.4999 && $1<0.5001 && $2==0 && $3>-0.0101 && $3<-0.0099 \
        && $4==0 && $5==0 && $6==0 && $7>-0.0101 && $7<-0.0099 && $8==0 \
        && $9>0.0099 && $9<0.0101 && $10==0 && $11==0 && $12==0 && $13==0 \
        {ok=1} END{exit !ok}' "$T2014/beam_force_moment.5" 2>/dev/null; then
  check_ok "bmom1 (analitico: dist 0.5 fy1=-F mz1=-F*L fy2=+F mz2=0)"
else
  check_fail "bmom1" "valores inesperados: [$(cat "$T2014/beam_force_moment.5" 2>/dev/null)]"
fi

# bmom_switch: -separate_sequential (archivo beam_force_moment.0) con
# switch -yes -> los 12 componentes con el signo INVERTIDO respecto a
# bmom1 (fy1 = +0.01, mz1 = +0.01, fy2 = -0.01)
if [ "$(wc -l < "$T2014/beam_force_moment.0" 2>/dev/null)" = "1" ] && \
   awk 'NF==13 && $1>0.4999 && $1<0.5001 && $3>0.0099 && $3<0.0101 \
        && $7>0.0099 && $7<0.0101 && $9>-0.0101 && $9<-0.0099 && $13==0 \
        {ok=1} END{exit !ok}' "$T2014/beam_force_moment.0" 2>/dev/null; then
  check_ok "bmom_switch (-yes: signos invertidos vs bmom1; archivo sequential .0)"
else
  check_fail "bmom_switch" "valores inesperados: [$(cat "$T2014/beam_force_moment.0" 2>/dev/null)]"
fi

# bmom_truss: la fuerza axial sale del truss (fx1 = +N, fx2 = -N) y el
# resto de columnas de la viga (fy1 = -F, mz1 = -F*L, fy2 = +F, mz2 =
# 0). N se verifica contra control_print_element -element_truss_force
# (element_truss_force_0.7: "x y N" con method -middle): ambos leen
# ELEMENT_TRUSS_FORCE y deben coincidir EXACTAMENTE.
NTR=$(awk 'NR==1{print $3}' "$T2014/element_truss_force_0.7" 2>/dev/null)
if [ "$(wc -l < "$T2014/beam_force_moment.7" 2>/dev/null)" = "1" ] && \
   awk -v ntr="$NTR" 'NF==13 && $1>0.4999 && $1<0.5001 && $2>0.019 && $2<0.021 \
        && $3>-0.0101 && $3<-0.0099 && $7>-0.0101 && $7<-0.0099 \
        && $8>-0.021 && $8<-0.019 && $9>0.0099 && $9<0.0101 && $13==0 \
        && ntr>0.019 && ntr<0.021 {ok=1} END{exit !ok}' \
        "$T2014/beam_force_moment.7" 2>/dev/null; then
  check_ok "bmom_truss (axial fx1=+N fx2=-N con N=$NTR == element_truss_force; fy/mz de la viga)"
else
  check_fail "bmom_truss" "axial o beam incorrectos (NTR=$NTR): [$(cat "$T2014/beam_force_moment.7" 2>/dev/null)]"
fi

# bmom_2d: 2 vigas colineales + corte DIAGONAL 2D (solo x,y) que cruza
# SOLO el elemento 2 (x=0.75): 1 linea, dist = sqrt(1.25) =
# 1.11803398875, mz1 = -F*0.5 = -0.005 (momento en x=0.5), fy2 = +F
if [ "$(wc -l < "$T2014/beam_force_moment.8" 2>/dev/null)" = "1" ] && \
   awk 'NF==13 && $1>1.1179 && $1<1.1182 && $2==0 && $3>-0.0101 && $3<-0.0099 \
        && $7>-0.0051 && $7<-0.0049 && $9>0.0099 && $9<0.0101 && $13==0 \
        {ok=1} END{exit !ok}' "$T2014/beam_force_moment.8" 2>/dev/null; then
  check_ok "bmom_2d (2D: corte diagonal cruza solo el elem 2: dist sqrt(1.25), mz1=-F*0.5)"
else
  check_fail "bmom_2d" "valores inesperados: [$(cat "$T2014/beam_force_moment.8" 2>/dev/null)]"
fi

# bmom_noint: el corte (x=1.5) no cruza ninguna viga -> NO se escribe
# archivo (decision documentada)
if [ ! -f "$T2014/beam_force_moment.9" ]; then
  check_ok "bmom_noint (corte sin vigas: archivo NO creado)"
else
  check_fail "bmom_noint" "beam_force_moment.9 no deberia existir"
fi

# ---------------------------------------------------------------------
# Sub-sprint materi_stress_force, lote 1 (infraestructura): registro +
# dispatch + print. A/B verificado con el calcul.cc de HEAD (binario
# temporal): post_calcul -materi_stress -force moria en
# db_error(POST_CALCUL,0) del dispatch de operats ("Error detected for
# data item : post_calcul"); ahora la rama FORCE existe en calculate()/
# calculate_operat (calcul.cc), los 10 records de configuracion
# post_calcul_materi_stress_force_* y control_print_materi_stress_force
# estan registrados, y el print escribe materi_stress_force.<index>.
# Lote 2: la integracion 2D rellena los valores; msf_parse conserva
# valores 0 porque su reference_point (0.5,0.5) COINCIDE con el
# centroide del elemento (caso degenerado documentado: warning + sin
# fuerzas). El aviso "not yet implemented" queda SOLO para 3D (lote 3).
# ---------------------------------------------------------------------

# msf_parse: 2D con 9 records de configuracion + print -> archivo
# materi_stress_force.100 con 3 lineas de cabecera (#) + 4 nodos; el
# warning de direction_* en 2D y el del centroide (reference_point en
# el centroide -> 0 fuerzas, caso documentado) aparecen en el stdout
if [ -f "$T2014/materi_stress_force.100" ] && \
   [ "$(grep -c '^#' "$T2014/materi_stress_force.100")" -ge 3 ] && \
   [ "$(grep -vc '^#' "$T2014/materi_stress_force.100")" = "4" ] && \
   [ "$(grep -c 'coincides with the centroid' /tmp/msf_parse_safe.out)" = "1" ] && \
   grep -q "direction" /tmp/msf_parse_safe.out; then
  check_ok "msf_parse (2D: cabecera comentada + 4 lineas de nodo + warning centroide unico + warning 2D)"
else
  check_fail "msf_parse" "materi_stress_force.100 o mensajes inesperados ($(ls "$T2014"/materi_stress_force.* 2>/dev/null))"
fi

# msf_parse_3d: 16 items 3D (MCALCUL=20 OK) -> materi_stress_force.200
# con 16 columnas de datos y 8 lineas de nodo
MSF200_NCOL=$(awk '!/^#/ && NF>0 {print NF; exit}' "$T2014/materi_stress_force.200" 2>/dev/null)
if [ -f "$T2014/materi_stress_force.200" ] && \
   [ "$(grep -vc '^#' "$T2014/materi_stress_force.200")" = "8" ] && \
   [ "$MSF200_NCOL" = "17" ]; then
  check_ok "msf_parse_3d (3D: 8 nodos x 16 items + columna nodo = $MSF200_NCOL columnas)"
else
  check_fail "msf_parse_3d" "materi_stress_force.200 inesperado ($MSF200_NCOL columnas)"
fi

# msf_print: -all (300) y -primary (301) -> 4 lineas cada uno (quad4:
# sin promediado, todos los nodos primarios); valores 0 porque el
# reference_point (0.5,0.5) coincide con el centroide del unico
# elemento (caso degenerado documentado) y cabecera comentada
N300=$(grep -vc '^#' "$T2014/materi_stress_force.300" 2>/dev/null || echo 0)
N301=$(grep -vc '^#' "$T2014/materi_stress_force.301" 2>/dev/null || echo 0)
if [ "$N300" = "4" ] && [ "$N301" = "4" ] && \
   awk '!/^#/ && NF==10 { for (i=2;i<=NF;i++) if ($i!=0) bad=1 } END{exit bad}' \
   "$T2014/materi_stress_force.300" 2>/dev/null; then
  check_ok "msf_print (-all 4 lineas == -primary 4 lineas; 9 items 2D todos 0: reference_point en el centroide)"
else
  check_fail "msf_print" "conteo o valores inesperados (N300=$N300 N301=$N301)"
fi

# msf_errors_2dwarn: 2D sin reference_point -> warning y rc=0 (run OK)
if grep -q "default reference point" /tmp/msf_errors_2dwarn_safe.out; then
  check_ok "msf_errors_2dwarn (2D sin reference_point: warning + default documentado, rc=0)"
else
  check_fail "msf_errors_2dwarn" "faltaba el warning del default"
fi

# msf_errors_*: las validaciones fallan con rc=1 y un mensaje claro
# (cada una se ejecuta fuera del bucle: son errores esperados)
msf_error_ok() {
  local t="$1" msg="$2"
  ( cd validation-suite/test-2014 &&
    ulimit -v "$TN_MEMLIMIT_KB" &&
    timeout 60 "$REPO_DIR/build/$TN_BIN" "$t.dat" > "/tmp/${t}_safe.out" 2>&1 )
  local rc=$?
  if [ "$rc" = "1" ] && grep -q "$msg" "/tmp/${t}_safe.out"; then
    check_ok "$t (rc=1 + '$msg')"
  else
    check_fail "$t" "rc=$rc, sin mensaje '$msg'"
  fi
}
msf_error_ok msf_errors_nogroup  "requires post_calcul_materi_stress_force_element_group"
msf_error_ok msf_errors_bothdir  "are mutually exclusive"
msf_error_ok msf_errors_3dnodir  "requires either"
msf_error_ok msf_errors_3dnoref  "requires.*reference_point"
msf_error_ok msf_errors_refcount "reference_point needs one point"
msf_error_ok msf_errors_posttype "NODAL calculation"

# ---------------------------------------------------------------------
# Expected-error tests (generic marker, 2026-10-01). A .dat whose comment
# holds a line "EXPECTED_ERROR: <text>" PASSES only when the run stops
# with rc=1 AND its output contains <text> (fixed string); a clean run
# or another error FAILS. Corpus classifiers must treat these .dat the
# same way (PASS = rc=1 + marker text), not as RUNFAIL.
#   materi_direct_auto (test-2014): single hex8 with
#     _direct_normal_automatic -yes -> zero element_normal input error
#     (Pro aborts on an assertion).
#   element_normal2d_row_error1 (test-2026): single row of 3 quad4,
#     middle element normal zero -> input error naming element 2.
# The two companions of the row test (explicit normal, automatic group
# only on the end elements) must run with rc=0 (targets met).
# ---------------------------------------------------------------------
expected_error_ok() {
  local dir="$1" t="$2" msg
  msg=$(sed -n 's/^[[:space:]]*EXPECTED_ERROR:[[:space:]]*\(.*[^ )]\)[[:space:]]*)*[[:space:]]*$/\1/p' \
    "$dir/$t.dat" | head -1)
  if [ -z "$msg" ]; then
    check_fail "$t" "no EXPECTED_ERROR marker in $t.dat"
    return
  fi
  ( cd "$dir" &&
    ulimit -v "$TN_MEMLIMIT_KB" &&
    timeout 60 "$REPO_DIR/build/$TN_BIN" "$t.dat" > "/tmp/${t}_safe.out" 2>&1 )
  local rc=$?
  if [ "$rc" = "1" ] && grep -qF "$msg" "/tmp/${t}_safe.out"; then
    check_ok "$t (expected error: rc=1 + '$msg')"
  else
    check_fail "$t" "rc=$rc, expected error '$msg' not found"
  fi
}
run_ok() {
  local dir="$1" t="$2"
  ( cd "$dir" &&
    ulimit -v "$TN_MEMLIMIT_KB" &&
    timeout 60 "$REPO_DIR/build/$TN_BIN" "$t.dat" > "/tmp/${t}_safe.out" 2>&1 )
  local rc=$?
  if [ "$rc" = "0" ]; then check_ok "$t (rc=0, targets met)";
  else check_fail "$t" "rc=$rc"; fi
}
T2026="$REPO_DIR/validation-suite/test-2026"
expected_error_ok "$T2014" materi_direct_auto
expected_error_ok "$T2026" element_normal2d_row_error1
run_ok "$T2026" element_normal2d_row_explicit1
run_ok "$T2026" element_normal2d_row_layer1
expected_error_ok "$T2026" element_normal2d_row_mixed_error1
run_ok "$T2026" element_normal2d_row_mixed1

# ---------------------------------------------------------------------
# C16 default group_materi_memory (2026-10-04): an unspecified memory
# behaves as -updated_linear, as in the Professional 02-08-2026.
#   materi_memory_default_shear1: single quad4 simple shear gamma=1,
#     no memory record -> sigxy 0.5, sigxx = sigyy = 0 (= Pro; before
#     C16 GNU gave the -updated values 0.2311/0.4207/-0.2286).
#   materi_memory_default_stretch1: -follow_material uniaxial stretch
#     to twice the length, no memory record -> sigxx 0.695653 (= Pro):
#     the default must keep the current-geometry gradient (a
#     reference-geometry gradient would give 1).
#   materi_memory_default_updated_apply1: same shear, no memory record,
#     control_materi_updated_apply 0 -yes -> -updated values sigxx
#     0.231098, sigxy 0.420739, sigyy -0.228598 (manual Professional
#     6.152). DOCUMENTED DIVERGENCE: the Professional 02-08-2026 binary
#     ignores -yes for a non-specified memory (0 / 0.5 / 0).
# ---------------------------------------------------------------------
run_ok "$T2026" materi_memory_default_shear1
run_ok "$T2026" materi_memory_default_stretch1
run_ok "$T2026" materi_memory_default_updated_apply1

# ---------------------------------------------------------------------
# F3 kinematics (C8, C18, C19), values = Professional 02-08-2026.
#   mesh_default_fixed_in_space1/2: no mesh record -> -fixed_in_space
#     (manual 6.806): velocity-only bar stretched 50% sigxx 0.5 (old GNU
#     -follow_material default 0.407136); rigid rotation stresses 0 (old
#     0.109567).
#   materi_memory_default_displacement1: materi_displacement and no
#     memory record -> -total_linear (manual 6.797): column node-3 disx
#     0.5 (old 0.472222).
#   materi_memory_updated_linear_current1/2: explicit -updated_linear on
#     -follow_material uses the current geometry (manual 2.2.1): bar
#     sigxx 0.407136 (old reference gradient 0.5); column velix 0.52093
#     within 1.e-3 (GNU 0.52048, in-step geometry divergence; old 0.472).
# ---------------------------------------------------------------------
run_ok "$T2026" mesh_default_fixed_in_space1
run_ok "$T2026" mesh_default_fixed_in_space2
run_ok "$T2026" materi_memory_default_displacement1
run_ok "$T2026" materi_memory_updated_linear_current1
run_ok "$T2026" materi_memory_updated_linear_current2
# F3/C20: group_materi_memory -updated_jaumann (manual 2.2.1, 6.686),
#   simple shear gamma 0.5: 1 step 0.125 / 0.234375 / -0.125, 100 steps
#   0.0618335 / 0.2396333 / -0.0618335 (= Professional 02-08-2026; the
#   old GNU rejected the value).
run_ok "$T2026" materi_memory_updated_jaumann_shear1
run_ok "$T2026" materi_memory_updated_jaumann_shear2

# ---------------------------------------------------------------------
# F4 T1 TPL + K10 (2026-10-09): groundflow_total_pressure_limit as in
# manual 6.587 and the Professional 02-08-2026: default limit 0, the
# CALCULATED total pressure is capped at every node (prescribed ones
# included), the pres dof is never cut, and the groundflow velocity is 0
# where the integration-point total pressure reaches the limit.
#   groundflow_total_pressure_limit_default1/high1/negative1: one quad4,
#     consolidation source, free-node pres 1.6667 always; total pressure
#     / gvely: no record 0 / 0, limit 100 1.6667 / 0.16667, limit -0.5
#     -0.5 at every node / 0 (old GNU: gvely 0.16667, total pressure 0,
#     pres -0.5).
#   groundflow_velocity_dry1/2: corpus ground1 (K10), gvely 0 without the
#     record (old GNU -0.05), -0.05 with limit 1.e10.
# ---------------------------------------------------------------------
run_ok "$T2026" groundflow_total_pressure_limit_default1
run_ok "$T2026" groundflow_total_pressure_limit_high1
run_ok "$T2026" groundflow_total_pressure_limit_negative1
run_ok "$T2026" groundflow_velocity_dry1
run_ok "$T2026" groundflow_velocity_dry2
# F4 T2 / C41 (2026-10-09): a tied VELOCITY slave carries its integrated
#   displacement (mpc_node_apply). mohr_coul_direct_mpc_total1 = Pro
#   other/mohr_coul_direct7 (-total_linear default + velocity tie):
#   sigxy 0.0998334 at nodes 1 and 2 (= Pro; old GNU 0.102732 /
#   0.096935, the slave disy lagged one solve).
run_ok "$T2026" mohr_coul_direct_mpc_total1
# F4 T3 / C33 (2026-10-09): set_stress reads DTIME at VERSION_NEW and the
#   principal tension cap honours _visco; each direct law relaxes with its
#   own _visco record (values = Professional 02-08-2026):
#   materi_direct_visco_first_step1: tension_direct, 1 step dt 1: sigxx
#     1.849143e-3 (old 1e-3); _first_step2: compression_direct, dt 0.5
#     then 1: -2.550147e-3 (old -3.024897e-3); _separate1: full MC + tension
#     cap, tm 2 / 4: sigxx -2.654391e-4, sigxy 1.688146e-3 (old -1.2e-3 /
#     1.2e-3); materi_plasti_visco_power_dt1: 1D visco_power, dt 0.5 then 1:
#     sigxx 2.666667 (old 2.818182).
run_ok "$T2026" materi_direct_visco_first_step1
run_ok "$T2026" materi_direct_visco_first_step2
run_ok "$T2026" materi_direct_visco_separate1
run_ok "$T2026" materi_plasti_visco_power_dt1
# F4 T4 / K20 (2026-10-09): no forced 1-point rule for axisymmetric
#   materi_velocity (polynom.cc pol(), calcul_force.cc msf_element_rule);
#   quad4 2x2 / quad9 3x3 Gauss as the Professional 02-08-2026:
#   axisym_quad4_hourglass1: 4x4 quad4 ring patch, point load: post_point
#     disx 2.657533803232e-4, disy -4.380511327008e-4 (= Pro; old GNU
#     hourglass ~1e12, solver did not converge);
#   axisym_quad4_ring1: Akin thick ring: disx 1.081776788020e-2 (= Pro),
#     sigzz 10760.09 vs Pro 10759.22 +- 1 (quad4 recovery, BACKLOG C42;
#     old GNU 1.08238e-2 / 10000);
#   axisym_quad9_uniform1: quad9 no longer refused: sigxx = sigzz = 1000,
#     disx 500 (analytic; old GNU "not available for axisymmetric").
run_ok "$T2026" axisym_quad4_hourglass1
run_ok "$T2026" axisym_quad4_ring1
run_ok "$T2026" axisym_quad9_uniform1
# F4 T4 / C31 (2026-10-09): tet10 Keast 5-point rule (polynom.cc), linear
#   least-squares recovery 1,L1,L2,L3 (miscel.cc), HRZ lumping (general.cc):
#   tet10_body_load1: 48 tet10 under gravity: sigzz -1 (= Pro tet10_2;
#     old GNU -0.894);
#   tet10_mass_hrz1: one tet10, unit point forces, dt 0.01: v = dt/m with
#     the HRZ masses, corner 2.16, mid-edge 0.405 (GNU-specific target,
#     deliberate divergence TET10-HRZ: Pro row sum -1.2 / 0.3; old GNU
#     59999.99 / 0.36);
#   tet10_dynamic_gravity1: 48 tet10 column, sudden gravity with inertia:
#     stable, top disz -0.5286 at t 3 (analytic -0.5 +- 0.05; row sum with
#     the 5-point rule diverges to 5e8, the Pro solver fails).
run_ok "$T2026" tet10_body_load1
run_ok "$T2026" tet10_mass_hrz1
run_ok "$T2026" tet10_dynamic_gravity1
# F4 T7 / C (2026-10-09): one-shot mesh operations (control_mesh_macro,
#   _copy, _mirror, _switch) on the SAME index as control_timestep run on
#   the first step of the control only (top.cc step_start_control_entry),
#   as the Professional 02-08-2026:
#   mesh_macro_timestep1: Pro mesh_activate1 with 20 steps: 202 nodes
#     (= Pro; old GNU 4040, one macro copy per step);
#   mesh_switch_timestep1: two steps, axes swapped once: node 3
#     1 0 -> 0 1 (= Pro; old GNU swapped twice, node 3 unchanged; before
#     the F4 T8 macro numbering the GNU node 3 was 0.5 0 -> 0 0.5).
run_ok "$T2026" mesh_macro_timestep1
n_nodes=$(grep -ac '^node ' "$T2026/mesh_macro_timestep1.dbs" 2>/dev/null)
if [ "$n_nodes" = "202" ]; then check_ok "mesh_macro_timestep1 (202 nodes)";
else check_fail "mesh_macro_timestep1" "nodes=$n_nodes, expected 202"; fi
run_ok "$T2026" mesh_switch_timestep1
if grep -aqE '^node +3 +0 +1 *$' "$T2026/mesh_switch_timestep1.dbs" 2>/dev/null;
then check_ok "mesh_switch_timestep1 (node 3 at 0 1)";
else check_fail "mesh_switch_timestep1" "node 3 not at 0 1"; fi
# F4 T8 / C45 (2026-10-10): control_mesh_copy / _mirror on an index after
#   a control_timestep (mesh.cc mesh_duplicate_node_data /
#   mesh_duplicate_element_state). Old GNU: rc 1 (node_dof_start_refined
#   of the new nodes missing), then 1 % off (integration-point state
#   re-initialised). Targets = Professional 02-08-2026 post_point values.
run_ok "$T2026" mesh_copy_timestep1
run_ok "$T2026" mesh_mirror_timestep1
# F4 T8 / C46d (2026-10-10): control_mesh_macro -bar/-rectangle/-brick
#   numbered as the Professional 02-08-2026 (x fastest, then y, then z;
#   macro.cc). Old GNU: transposed (x slowest).
run_ok "$T2026" mesh_macro_numbering1
if grep -aqE '^node +2 +1 +0 *$' "$T2026/mesh_macro_numbering1.dbs" 2>/dev/null &&
   grep -aqE '^element +1 +-quad4 1 2 4 5 *$' "$T2026/mesh_macro_numbering1.dbs" 2>/dev/null;
then check_ok "mesh_macro_numbering1 (node 2 at 1 0, element 1 1 2 4 5)";
else check_fail "mesh_macro_numbering1" "rectangle numbering differs from the Professional"; fi
run_ok "$T2026" mesh_macro_numbering2
if grep -aqE '^node +2 +1 +0 +0 *$' "$T2026/mesh_macro_numbering2.dbs" 2>/dev/null &&
   grep -aqE '^element +1 +-hex8 1 2 4 5 7 8 10 11 *$' "$T2026/mesh_macro_numbering2.dbs" 2>/dev/null;
then check_ok "mesh_macro_numbering2 (node 2 at 1 0 0, element 1 1 2 4 5 7 8 10 11)";
else check_fail "mesh_macro_numbering2" "brick numbering differs from the Professional"; fi
# F4 T8 / C46a: one-shot mesh records on one index in the Professional
#   order copy, (move), switch, mirror (top.cc step_start).
run_ok "$T2026" mesh_order_timestep1
if grep -aqE '^node +10 +0 +7 *$' "$T2026/mesh_order_timestep1.dbs" 2>/dev/null &&
   grep -aqE '^node +25 +-1 +5 *$' "$T2026/mesh_order_timestep1.dbs" 2>/dev/null;
then check_ok "mesh_order_timestep1 (node 10 at 0 7, node 25 at -1 5)";
else check_fail "mesh_order_timestep1" "copy/switch/mirror order differs from the Professional"; fi
# F4 T8 / C46b: control_mesh_move once per control index, node and
#   node_start_refined moved (manual 6.217; the Professional moves only
#   node_start_refined, DIVERGENCES C46b). Old GNU: once per step.
run_ok "$T2026" mesh_move_timestep1
if grep -aqE '^node +6 +12 +1 *$' "$T2026/mesh_move_timestep1.dbs" 2>/dev/null &&
   grep -aqE '^node_start_refined +6 +12 +1 *$' "$T2026/mesh_move_timestep1.dbs" 2>/dev/null;
then check_ok "mesh_move_timestep1 (node 6 at 12 1, moved once)";
else check_fail "mesh_move_timestep1" "node 6 not at 12 1"; fi
# F4 T8: control_mesh_switch with two entries swaps them (manual 6.244,
#   Professional: -x -y swaps x and y; old GNU: axis order = no change).
run_ok "$T2026" mesh_switch_pair1
if grep -aqE '^node +3 +0 +2 *$' "$T2026/mesh_switch_pair1.dbs" 2>/dev/null;
then check_ok "mesh_switch_pair1 (node 3 at 0 2)";
else check_fail "mesh_switch_pair1" "node 3 not at 0 2"; fi
# F4 T8 / C46d: mesh_activate_gravity_time read at its own index, bottom-
#   to-top activation, default -method2, reduced tangent and stress
#   increment (mesh.cc mesh_activate_gravity_factor, materi.cc). Targets =
#   Professional values; old GNU: full gravity from the first step.
run_ok "$T2026" mesh_activate_gravity_time1
# F4 T9 / C43 (2026-10-10): -matrix_superlu route (so_suplu.c). The old
#   code read L and U (dQuerySpace) after destroying them and free()d the
#   new[] matrix arrays (ASan alloc-dealloc-mismatch). Old GNU: segmentation
#   fault (rc 139) after the first solve. Analytic plane-strain targets.
run_ok "$T2026" superlu_solve1
# F4 T9 / C43 (2026-10-10): -matrix_pardiso is solved with SuperLU (so.cc;
#   old GNU: always the Bi-CG) and falls back to the Bi-CG when the SuperLU
#   factorization fails. pardiso_bicg_fallback1 forces the failure with the
#   memory guard budget TN_PARDISO_SUPERLU_MB=4 (4 MB above the current
#   address space; SuperLU needs ~130 MB). Same analytic targets.
run_ok "$T2026" pardiso_superlu1
if grep -aq "Solution Found" /tmp/pardiso_superlu1_safe.out 2>/dev/null &&
   ! grep -aq "falling back to the iterative Bi-CG" /tmp/pardiso_superlu1_safe.out 2>/dev/null;
then check_ok "pardiso_superlu1 (-matrix_pardiso solved by SuperLU)";
else check_fail "pardiso_superlu1" "SuperLU did not solve -matrix_pardiso"; fi
( cd "$T2026" &&
  ulimit -v "$TN_MEMLIMIT_KB" &&
  TN_PARDISO_SUPERLU_MB=4 timeout 60 "$REPO_DIR/build/$TN_BIN" \
    pardiso_bicg_fallback1.dat > /tmp/pardiso_bicg_fallback1_safe.out 2>&1 )
RC_PBF=$?
if [ "$RC_PBF" = "0" ] &&
   grep -aq "SuperLU failed for -matrix_pardiso - falling back" /tmp/pardiso_bicg_fallback1_safe.out 2>/dev/null;
then check_ok "pardiso_bicg_fallback1 (SuperLU out of memory -> Bi-CG, targets met)";
else check_fail "pardiso_bicg_fallback1" "rc=$RC_PBF or no Bi-CG fallback"; fi
# F4 T9 / C43: the iterative solver multiplies with an assembled CSR matrix
#   (so_bicg.cc bicg_csr_build/bicg_csr_sys) instead of the element-by-
#   element products; two threads (options_processors 2). print_solver -yes
#   reports the operator; old GNU: no CSR line. Analytic targets.
run_ok "$T2026" bicg_csr_threads1
if grep -aq "assembled CSR operator" /tmp/bicg_csr_threads1_safe.out 2>/dev/null;
then check_ok "bicg_csr_threads1 (assembled CSR operator, 2 threads)";
else check_fail "bicg_csr_threads1" "no assembled CSR operator"; fi
# F4 T10 / C48 (2026-10-10): element_interface_intpnt_materi_tension_status
#   is an INTEGER record; interface.cc interface_element() PUT it from the
#   wrong buffer (idum[1] with length n_intpnt), so the statuses were stack
#   garbage (old GNU: "2 -total_linear 3"). Values = Professional
#   02-08-2026 on this input (-closed x3 / -opened x3); two runs must give
#   byte-identical .dbs files.
run_ok "$T2026" interface_tension_status1
ITS_DBS="$T2026/interface_tension_status1.dbs"
if grep -aqE '^element_interface_intpnt_materi_tension_status +1 +-closed -closed -closed *$' "$ITS_DBS" 2>/dev/null &&
   grep -aqE '^element_interface_intpnt_materi_tension_status +2 +-opened -opened -opened *$' "$ITS_DBS" 2>/dev/null;
then check_ok "interface_tension_status1 (-closed / -opened per integration point)";
else check_fail "interface_tension_status1" "wrong element_interface_intpnt_materi_tension_status"; fi
cp "$ITS_DBS" /tmp/interface_tension_status1_run1.dbs 2>/dev/null
run_ok "$T2026" interface_tension_status1
if cmp -s "$ITS_DBS" /tmp/interface_tension_status1_run1.dbs;
then check_ok "interface_tension_status1 (two runs byte-identical)";
else check_fail "interface_tension_status1" "two runs give different .dbs"; fi
# F4 T5 follow-up (2026-10-10): the hardsoil yield function reads gamma_p
#   = kappa at the step-start value plus the increment of its own return
#   mapping (stress.cc set_stress pl_unknowns), not the nodal kappa of the
#   previous material evaluation of the step. Same result for 1 and 2
#   equilibrium iterations (kap 0.00104991, sigxx -3.11898); old GNU
#   0.0030 / 0.0006 (period-2 oscillation in the iteration count).
run_ok "$T2026" hardsoil_kappa_iterations1
run_ok "$T2026" hardsoil_kappa_iterations2
# F4 T5 / C35 (2026-10-10): post-solve material state pass (top.cc): the
#   material is evaluated once more on the final iterate of every step, so
#   the reported/stored state follows the solved velocity (= Professional
#   02-08-2026). One equilibrium iteration per step in all three:
#   stretch1 node 2 velx 0.5 in every step (old GNU 0.5 1 1 0.5 0 0),
#   plastic1 eppxx 0.01 / sigxx 0.015 at node 1 (old GNU 0 / 0),
#   rigid1 eptxx = sigxx = 0 (old GNU -0.01).
rm -f "$T2026/velx2.his"
run_ok "$T2026" post_solve_state_stretch1
if [ -f "$T2026/velx2.his" ] &&
   [ "$(wc -l < "$T2026/velx2.his")" = "6" ] &&
   awk '{ d=$2-0.5; if (d<0) d=-d; if (d>1.e-10) bad=1 } END { exit bad }' "$T2026/velx2.his";
then check_ok "post_solve_state_stretch1 (node 2 velx 0.5 in all 6 steps)";
else check_fail "post_solve_state_stretch1" "node 2 velx history is not 0.5 in every step"; fi
run_ok "$T2026" post_solve_state_plastic1
run_ok "$T2026" post_solve_state_rigid1

# ---------------------------------------------------------------------
# K03 (2026-10-06): default convection_apply -no (manual Professional
# 6.395, = Professional 02-08-2026). sfnet condif2 model in the Pro
# dialect: without the record the heat is not convected by the material
# velocity, middle temp = 0.5 (pure diffusion); with -yes temp = 1.
# Before K03 the GNU default -yes gave 1 for both.
# ---------------------------------------------------------------------
run_ok "$T2026" convection_apply_default1
run_ok "$T2026" convection_apply_yes1

# ---------------------------------------------------------------------
# C24 (2026-10-05): control_print_frd / control_print_gmsh loop over the
# PRIMARY unknowns (as print_vt.cc and print_dof after C23). Single quad4
# with derivatives: the FRD must hold exactly VELO DISP STRESS TOSTRAIN
# NDTEMP (once each), node 2 D1 = 1e-2 and D3 = 0 (2D), SXX = 10; the
# gmsh file 17 node_ fields (2+2+6+6+1 components), no heap garbage.
# Before C24: 68 FRD blocks and 324 node_ fields (shifted duplicates).
# The .frd is appended to when it exists, so both outputs are removed
# before the run.
# ---------------------------------------------------------------------
rm -f "$T2026/print_frd_gmsh_dof_loop1.frd" "$T2026/print_frd_gmsh_dof_loop1.msh"
run_ok "$T2026" print_frd_gmsh_dof_loop1
FRDL="$T2026/print_frd_gmsh_dof_loop1.frd"
MSHL="$T2026/print_frd_gmsh_dof_loop1.msh"
if [ "$(awk '$1=="-4"{printf "%s ", $2}' "$FRDL" 2>/dev/null)" = "VELO DISP STRESS TOSTRAIN NDTEMP " ] && \
   awk '$1=="-4"{b=$2} b=="DISP" && $1=="-1" && $2==2 { if ($3!="1.00000E-02" || $5!="0.00000E+00") bad=1; n++ }
        b=="STRESS" && $1=="-1" && $2==2 { if ($3!="1.00000E+01") bad=1; n++ }
        END { exit (bad || n!=2) }' "$FRDL" && \
   [ "$(grep -c '^"node_' "$MSHL" 2>/dev/null)" = "17" ] && \
   [ "$(grep '^"node_' "$MSHL" | sort -u | wc -l)" = "17" ] && \
   ! rg -q 'e-3[0-9][0-9]' "$FRDL" "$MSHL"; then
  check_ok "print_frd_gmsh_dof_loop1 (C24: 5 FRD blocks, 17 gmsh node_ fields, D3=0, SXX=10)"
else
  check_fail "print_frd_gmsh_dof_loop1" "FRD/gmsh blocks not one per primary dof (C24)"
fi

# ---------------------------------------------------------------------
# K04 condif7 (2026-10-06): the _node / _element_node restriction of the
# edge families (condif heat/convection/radiation, groundflow flux,
# force_edge*, support_edge_normal) selects the SIDES whose nodes are
# all in the list (= Professional 02-08-2026). Before: node by node, so
# a corner node also loaded the adjacent side (condif7 flux doubled).
#   edge_node_restriction_heat1: nodes 1 3 4 -> temp 12 / 18 (Pro);
#     before 20 / 20.
#   edge_node_restriction_force1: force_edge_node 1 3 -> disx 1 / 1
#     (Pro); before 2 / 2.
# ---------------------------------------------------------------------
run_ok "$T2026" edge_node_restriction_heat1
run_ok "$T2026" edge_node_restriction_force1

# ---------------------------------------------------------------------
# K04 force18 (2026-10-06): force_point is an external load along
# +force (manual Professional 6.501, = Professional 02-08-2026). Unit
# quad4 bar, force_point (1, 0.5) force (2, 0) -> disx(2) = disx(4) = 2
# (before: -2, the force was subtracted from node_rhside).
# ---------------------------------------------------------------------
run_ok "$T2026" force_point_sign1

# ---------------------------------------------------------------------
# K04 ground9 (2026-10-06): group_interface_groundflow_permeability and
# group_interface_condif_conductivity are per unit length/area (manual
# Professional 6.628/6.625): each facing pair carries w_i*L (as the
# mechanical interface springs); the interface capacity likewise. Middle pres/temp 0.5 (= Pro 02-08-2026);
# before: 0.6 (ground9 3 wide) and 1/3 (thermal ground9).
# ---------------------------------------------------------------------
run_ok "$T2026" interface_flux_measure_ground1
run_ok "$T2026" interface_flux_measure_condif1
# interface capacity: per unit length, split over the pair, implicit;
# groundflow_interface at t = 3 -> pres(5) 1.75 (Pro 1.749998578),
# before 0.8146.
run_ok "$T2026" interface_capacity_transient1

# ---------------------------------------------------------------------
# K04 stress_pressure_factor (2026-10-06): a step starting below the
# pressure history estimates its end pressure with the FACTORED stiffness
# and stays unloading/reloading only strictly below the history (rule
# measured on the Professional 02-08-2026, dt 0.1/0.15/0.25/0.3/0.35/0.4,
# factors 0.5 and 3). reload1 = Pro test (sigxx -3.05; before -3.10),
# reload2 = factor 3, partial unloading, dt 0.3 (-2.8; before -3.4).
# ---------------------------------------------------------------------
run_ok "$T2026" stress_pressure_factor_reload1
run_ok "$T2026" stress_pressure_factor_reload2

# ---------------------------------------------------------------------
# K04 incnav3 (2026-10-06): group_materi_expansion_volume reduces the
# density in the momentum inertia (general.cc) as in the gravity force
# (materi.cc): Pro incnav3 tvely -4.990415335 (= Pro; before -2.4976).
# ---------------------------------------------------------------------
run_ok "$T2026" expansion_volume_inertia1

# ---------------------------------------------------------------------
# K13 (2026-10-06): the direct cut-offs (tension/MC/compression, full and
# plane mode) store their cut as plastic strain C^-1 (sig_trial - sig)
# and the direct-law yield value (0 after a cut; = Pro 02-08-2026, except
# the shear plastic strain: tensor value, Pro twice it).
# ---------------------------------------------------------------------
run_ok "$T2026" direct_plastic_strain_tension1
run_ok "$T2026" direct_plastic_strain_mc1
run_ok "$T2026" direct_yield_value1

# ---------------------------------------------------------------------
# K13 follow-up (2026-10-06): an element without a yield law reports
# materi_plasti_f 0 (Pro), not -1e10 (plasti_f_no_law1, from Pro
# tension_direct1).
# ---------------------------------------------------------------------
run_ok "$T2026" plasti_f_no_law1

# ---------------------------------------------------------------------
# K11 (2026-10-06): the strain tensors are recovered with the stress
# operator (strain_recovery_clamped1: clamped-edge eptyy 0 exact = Pro;
# before -8.6e-5).
# ---------------------------------------------------------------------
run_ok "$T2026" strain_recovery_clamped1

# ---------------------------------------------------------------------
# C28 (2026-10-06): node_force is an external load along +force (dof.cc
# adds it, as force_point): node_force_sign1 disx 2 (= Pro; before -2);
# control_mesh_cut_node_force stores the substitute with that sign
# (delete.cc): mesh_cut_node_force_sign1 node_force 6 = +1 (= Pro mesh_cut_1;
# before -1), remaining bar sigxx 1.
# ---------------------------------------------------------------------
run_ok "$T2026" node_force_sign1
run_ok "$T2026" mesh_cut_node_force_sign1

# ---------------------------------------------------------------------
# K06 (2026-10-06): group_materi_elasti_poisson_power replaces the
# poisson ratio of the elastic C (stress.cc clears the C of
# group_materi_elasti_young first). Oedometer sigxx -0.4389 / sigyy
# -1.4101 (Pro 02-08-2026 -0.4417 / -1.4123); before -1.4389 / -3.2924
# (C(young,0.3) + C(young,nu(p)), twice the stiffness).
# ---------------------------------------------------------------------
run_ok "$T2026" elasti_poisson_power_replace1

# ---------------------------------------------------------------------
# K17 support_edge_normal_3 (2026-10-06): support damping c and density
# d/dt enter the iteration matrix (area.cc). Pro support_edge_normal_3
# without the freeze record: disy 5.00998004e-4 (= Pro); before 124.75.
# ---------------------------------------------------------------------
run_ok "$T2026" support_edge_normal_damping_tangent1

# ---------------------------------------------------------------------
# K17 truss14 / truss2 (2026-10-06): truss -total is geometrically
# nonlinear (current axis, N = EA (L-L0)/L0) and the nonlinear truss
# memories carry the tension geometric stiffness N/L (I - n n) in the
# matrix (truss.cc). truss_total_rotation1 (= Pro truss14): disx
# 0.707107, disy 1.707107; truss_updated_transverse1 (= Pro truss2):
# reactions 0.1/0.1, veliy 0.82358. Before: solver breakdown, rc 1.
# ---------------------------------------------------------------------
run_ok "$T2026" truss_total_rotation1
run_ok "$T2026" truss_updated_transverse1

# ---------------------------------------------------------------------
# C39 truss8 (2026-10-06): a node whose elements are all inactive by
# group_time is not driven by its bounda_dof (bounda.cc
# bounda_node_dormant). Node 1 disx 0.1 (= Pro truss8); before 0.2.
# ---------------------------------------------------------------------
run_ok "$T2026" group_time_bounda_dormant1

# ---------------------------------------------------------------------
# K05 / C5 (2026-10-07): contact_target_geometry/_switch (manual
# Professional 6.106/6.107) is read for every index, read-only inside
# the parallel contact loop (contact_geometry_get, contact.cc). Before:
# index 0 was copied to contact_geometry with a db PUT inside the
# parallel loop (abort "Data is allocated in a parallel loop"), index
# > 0 was dropped. 1D Pro other/contac2: node 2 x = 1.0 (index 1 and
# index 0), 1.02 with the switch -yes (all = Professional 02-08-2026).
# ---------------------------------------------------------------------
run_ok "$T2026" contact_target_geometry_index1
run_ok "$T2026" contact_target_geometry_index0
run_ok "$T2026" contact_target_geometry_switch1

# ---------------------------------------------------------------------
# K05 (2026-10-07): element-target contact (contact_target_element_group,
# manual Professional 6.105) needs node_boundary (6.875); generated by
# area.cc contact_node_boundary_set() when absent. Before: no boundary,
# no contact. 1D Pro other/contac1: velx 0.02 / 0.01; 2D quad4 block:
# vely -0.01 (contacter) / -0.02 (target side) (= Professional 02-08-2026).
# ---------------------------------------------------------------------
run_ok "$T2026" contact_element_target_boundary1
run_ok "$T2026" contact_element_target_boundary2d1

# ---------------------------------------------------------------------
# F4 T12 / BACKLOG C47 (2026-10-11): control_mesh_extrude quad9 -> hex27
# (mesh.cc mesh_extrude). The inner z planes (mid-planes and layer tops)
# hold ONE node per source node, shared by the neighbouring hex27, as in
# the Professional. Before: per-element plane nodes, the elements were
# disconnected at the inner planes (force11 ring 7.9x too soft). Two
# quad9 x two hex27 layers: 75 node records (old GNU 84); nodal
# displacements = Professional 02-08-2026 (targets, 1e-9).
# ---------------------------------------------------------------------
run_ok "$T2026" extrude_quad9_shared_nodes1
n_ext=$(grep -c '^ *node ' "$T2026/extrude_quad9_shared_nodes1.dbs" 2>/dev/null)
if [ "$n_ext" = "75" ]; then
  check_ok "extrude_quad9_shared_nodes1 (75 nodes: 15 per z plane)"
else
  check_fail "extrude_quad9_shared_nodes1" "node records $n_ext, expected 75"
fi
# F4 T6f C57: post_point in the quad9 -> hex27 extruded mesh (inner-plane
# node_start_refined z; Professional post_point_dof targets)
run_ok "$T2026" post_point_extrude_hex27_1
# C71 (F4 T6f, 2026-10-13): axisymmetric -lobatto element on the axis
# r=0 integrated with Gauss (= Professional); before NaN.
run_ok "$T2026" axisym_lobatto_axis1

# ---------------------------------------------------------------------
# F4 T6a (2026-10-11), cheap batch.
# C55: control_print_dof_rhside with derivatives (print_un.cc): the nodal
#   right-hand side is indexed per principal unknown. Old GNU: column of
#   another dof (vely_rhside = 0) and a heap read past NODE_RHSIDE.
#   Values = Professional 02-08-2026.
run_ok "$T2026" dof_rhside_derivatives1
if awk '$1==1 && $2==1 { d=$3-1.04166666667; if (d<0) d=-d; f=(d<1.e-6) } END { exit !f }' \
     "$T2026/vely_rhside.31" 2>/dev/null &&
   awk '{ n++; d=$3+2.08333333333; if (d<0) d=-d; if (d>1.e-6) bad=1 } END { exit (bad || n!=4) }' \
     "$T2026/sigyy_rhside.31" 2>/dev/null; then
  check_ok "dof_rhside_derivatives1 (vely_rhside 1.0416667, sigyy_rhside -2.0833333)"
else
  check_fail "dof_rhside_derivatives1" "rhside columns differ from the Professional"
fi
# C56: words longer than MCHAR-1 (input.cc input_word): split inside a
#   comment (rc 0), an input error elsewhere. Old GNU: buffer overflow
#   (the data-word case crashed, rc 139).
run_ok "$T2026" input_long_comment_word1
expected_error_ok "$T2026" input_long_word_error1
# F4 T6f ranges (range.cc, Professional user manual "Ranges"): index
#   -step (old GNU: wrong keyword sign, junk indices), -from a -to a,
#   two -from/-to blocks, defines inside a range (old raw cin reads);
#   descending -from/-to and a 9.0 in a value range are errors as in
#   the Professional; value -step is a guard (already right).
run_ok "$T2026" range_index_step1
run_ok "$T2026" range_index_single1
run_ok "$T2026" range_index_two_blocks1
run_ok "$T2026" range_index_define1
run_ok "$T2026" range_value_step1
expected_error_ok "$T2026" range_index_descending_error1
expected_error_ok "$T2026" range_value_descending_error1
expected_error_ok "$T2026" range_value_double_error1
# -size_tot / -size_dev (calcul.cc): sqrt(A:A) and sqrt(B:B) as the
#   manual and the Professional calcul1 (targets 1 and 0.816497); old GNU
#   sqrt(0.5 A:A).
run_ok "$T2026" post_calcul_size_tot1
# K14a: hyhis1/2/3/5/7 written back for the Masin laws (hypoplas.cc);
#   targets = Professional hypo7 at t 0.1. Old GNU: 0.
run_ok "$T2026" hypo_masin_history_slots1
# C6: the iface_gen family inputs carry no generated interface element
#   any more; the run must generate it (element 3, element_group 3 10).
for t in iface_gen iface_gen_method iface_gen_geom iface_gen_quad9; do
  if ! grep -qE '^element +3 ' "$T2014/$t.dat" &&
     grep -aqE '^element +3 +-quad4 ' "$T2014/$t.dbs" 2>/dev/null &&
     grep -aqE '^element_group +3 +10 *$' "$T2014/$t.dbs" 2>/dev/null; then
    check_ok "$t (interface element 3 generated, not in the input)"
  else
    check_fail "$t" "interface element 3 not generated by the run"
  fi
done
# K24 (F4 T6b, deliberate divergence class A, 2026-10-11): 3D Biot u-p
#   consolidation with solver_matrix_symmetric -yes. The GNU keeps the
#   non-symmetric matrix (Bi-CG), the Professional symmetrizes it and
#   drifts. Targets = Professional 02-08-2026 with -no at node 1:
#   k24_consol_3d1/2 (-yes/-no, t 2.2): pres -0.7085296, disx 0.0542011
#   (Pro -yes -0.678332); k24_consol_3d_undrained1/2 (k 1e-8, one step):
#   pres -0.847677 (Pro -yes -0.80318). Same targets for -yes and -no.
run_ok "$T2026" k24_consol_3d1
run_ok "$T2026" k24_consol_3d2
run_ok "$T2026" k24_consol_3d_undrained1
run_ok "$T2026" k24_consol_3d_undrained2
# F4 T6c (2026-10-11): control_restart -yes resets the integration-point
#   state (manual 6.368; the Professional deletes element_dof). Unit quad4,
#   10 steps to t 1, restart, 5 steps to t 0.5 = a fresh 0.5 run. Targets
#   at node 4 = Professional 02-08-2026: _ip1 Drucker-Prager sigyy
#   -24.18351 (GNU before the fix +19.65), _ip2 elastic sigyy -52.08333
#   (GNU before the fix +52.08333).
run_ok "$T2026" control_restart_reset_ip1
run_ok "$T2026" control_restart_reset_ip2
# F4 T6d C11 (2026-10-11): truss memory default and control_materi_dynamic
#   (truss.cc, general.cc materi_dynamic_effective, so.cc gate). Targets =
#   Professional 02-08-2026, node 2 disx unless noted (old GNU in brackets):
#   default -total_linear with materi_displacement (Pro dynamic3) 60.44025
#   [32.6583]; explicit -total_linear, inverting truss 60.44025 [-5962.6];
#   only control_materi_dynamic 10 0. 0.0958890 [6.4477], 0.5 0.932818
#   [1.45599]; no dynamic record 0.995505 [1.91225]; default -updated_linear
#   without materi_displacement, node 2 rhside velx 50 [-4.677].
run_ok "$T2026" truss_memory_default_displacement1
run_ok "$T2026" truss_total_linear_inversion1
run_ok "$T2026" control_materi_dynamic_only1
run_ok "$T2026" control_materi_dynamic_half1
run_ok "$T2026" truss_dynamic_no_record1
run_ok "$T2026" truss_memory_default_velocity1
# F4 T6d C59 (2026-10-11, owner decision): the Professional input checks
#   of group_truss_memory (truss.cc): -updated/-updated_linear with
#   materi_displacement and -total/-total_linear without it stop with an
#   input error naming both records.
expected_error_ok "$T2026" truss_memory_updated_displacement_error1
expected_error_ok "$T2026" truss_memory_total_velocity_error1
# F4 T6d C52 (2026-10-11): interface Mohr-Coulomb history = the CLAMPED
#   tangential force (interface.cc), so a slip reversal unloads
#   elastically. Quasi-static reversal, element 3 shear: analytic 4.0,
#   Professional 02-08-2026 4.000000083 (GNU before the fix 5.0).
run_ok "$T2026" interface_mc_slip_reversal1
# F4 T6e C53 (2026-10-11): post-solve material state pass for slide and
#   contact models (parallel_contact + slide() re-applied on the final
#   iterate) plus the slide law fixes (absent slide_plasti_tension = no
#   tension limit, cap c + |Fn|*tan(phi)).
#   slide_friction_cohesion1: slide4 with cohesion 0.002, end time 10:
#     node_slide_force 1 tangential -2.065482e-3 (Professional 02-08-2026
#     -2.065482222595e-3; HEAD -0.002, slide fix without the pass
#     -2.066157e-3).
#   contact_block_pass1: Pro-syntax contact_block (inertia, gravity,
#     contact penalty): node_dof vely/disy/sigyy = Professional (<= 7e-13),
#     node_rhside 1 vely ~0 (HEAD -0.180239).
run_ok "$T2026" slide_friction_cohesion1
run_ok "$T2026" contact_block_pass1
# F4 T6f WU1 (2026-10-11): nodal stress/strain recovery of the simplex
#   elements and the 1-point bar2 = ARITHMETIC MEAN of the adjacent
#   elements' values (Professional 02-08-2026, measured exact; GNU was
#   area/volume weighted). Two elements of different size and stiffness
#   sharing nodes; targets = Professional node_dof:
#   tria3 1-point (mean1) and 3-point nodal rule (mean2), tria6 nodal
#   rule (tolerance 2e-4 of the family scale: the tria6 element itself
#   differs ~1e-5), bar2 1-point (1D, gravity), tet4 1-point.
#   validation_14_pass1 (test-2014): the Professional validation_14 first
#   pass (DP, hole, tria3 after control_mesh_split); hole nodes 4 and 120
#   sigyy/eppyy/sigxx/eppxx = Professional within 1e-3 (GNU before the
#   fix: node 4 sigyy -11.23 vs -8.5195).
run_ok "$T2026" tria3_nodal_recovery_mean1
run_ok "$T2026" tria3_nodal_recovery_mean2
run_ok "$T2026" tria6_nodal_recovery_mean1
run_ok "$T2026" bar2_nodal_recovery_mean1
run_ok "$T2026" tet4_nodal_recovery_mean1
run_ok "$T2014" validation_14_pass1
# F4 T6f WU2 (2026-10-11): the hardening parameter kappa is an
#   integration-point state (ELEMENT_DOF), not the h-lumped nodal value
#   interpolated back; its nodal value uses the strain recovery weight.
#   kappa_ip_state1: two quad4 von Mises + kap hardening, reactions =
#   Professional within 1e-3 (before 5e-3 off), nodal kap within 5 %
#   (quad4 recovery operator, C42/C64). validation_3_kappa (test-2014):
#   the Professional validation_3 with Pro targets post_node_result vely
#   -129.302 (+-0.13; before -128.835), disx/disy within 0.1 %, kap node 1
#   within 3 % (before 11.6 % off).
run_ok "$T2026" kappa_ip_state1
run_ok "$T2014" validation_3_kappa
# C68 F4 T6f (2026-10-12): kappa_shear, the cap1 pc, the damage d, the
#   kinematic hardening rho and the Maxwell stresses are integration-
#   point states too (ELEMENT_DOF, nodal value with the strain recovery
#   weight), as kappa above. Four tria3 (1-point rule) on a non-uniform
#   field, targets = Professional 02-08-2026 (RED on the previous binary:
#   damage reactions 1.9e-1 off, rho 1.1e-1, kapsh 1.3e-2, Maxwell nodal
#   msigyy 8.4e-2). cap1 has no Professional test (C68b, model diverges).
run_ok "$T2026" damage_ip_state1
run_ok "$T2026" kappa_shear_ip_state1
run_ok "$T2026" rho_ip_state1
run_ok "$T2026" maxwell_ip_state1
# C68d F4 T6f (2026-10-12): the old-state read loop of materi() reused
#   the Maxwell chain index for old_rho / old_epi / old_epp (two or more
#   chains: wrong slots, reads past the 9-entry arrays). Two identical
#   chains must equal one chain with the summed stiffness (exact in the
#   Professional 02-08-2026 too); targets = the one-chain GNU run, tight
#   (1e-6 relative). Before the fix: node 1 reaction -13.32 vs -19.17.
run_ok "$T2026" maxwell_two_chains_state1
# F4 T6f (2026-10-12), class A divergences (owner decision 2026-10-12),
#   guard tests: GNU already follows the continuum / manual answer.
#   mmc_tension_vertex1/2 (MC-VERTEX-ZZ): test-2014 mmc_tension at dt 0.05
#   and 0.047; classic MC uniaxial tension, plane strain, nu 0: sigxx
#   46.19 (2c cos phi/(1+sin phi)), sigzz 0 +- 0.05, eppzz 0 +- 1e-4 (MC
#   vertex). The Professional gives a dt-dependent artefact (sigzz 0.974 /
#   eppzz -9.7e-4 at dt 0.05, 0.367 at 0.047).
#   membrane_eptzz_total1/2 (MEMBRANE-EPTZZ): plane-stress von Mises
#   uniaxial, dt 0.01 / 0.1: eptzz = accumulated elastic + plastic zz
#   strain, analytic -0.0088249 (+- 3e-4 / 2e-3: GNU eppyy/eppzz split is
#   first order in dt, BACKLOG C69). The Professional stores only the last
#   step's zz correction (-9.6e-5 / -7.9e-4).
run_ok "$T2026" mmc_tension_vertex1
run_ok "$T2026" mmc_tension_vertex2
run_ok "$T2026" membrane_eptzz_total1
run_ok "$T2026" membrane_eptzz_total2
#   hypo_wolff_igs_oedometer1/2 (K14b): test-2014 hypo3 (von Wolffersdorff
#   + intergranular strain, strain reversal) at dt 1e-3 / 1e-5; sigxx
#   converged reference -0.10466 (Python integration of the law), +- 2e-3
#   / 5e-4. GNU -0.1034 / -0.104235; the Professional -0.1629 at dt 1e-3
#   (time-discretisation error at the reversal; hyhis slots: K14 open).
#   damage_mazars_eq_two_positive1 (C7): hex8 diag(2e-4,-3e-4,1e-4), dam
#   0.28610207515 (manual eps_eq = sqrt(sum <eps_i>+^2)); the Professional
#   0.36409 (weights 4,2,1 on the ascending positive principal strains).
#   damage_mazars_vonmises1 (C7b): damage_stress1 bar + von Mises 3e6, dam
#   0.925372225695, eppxx 1.87863e-4; the Professional 0.96946 / 1.368e-4.
run_ok "$T2026" hypo_wolff_igs_oedometer1
run_ok "$T2026" hypo_wolff_igs_oedometer2
run_ok "$T2026" damage_mazars_eq_two_positive1
run_ok "$T2026" damage_mazars_vonmises1

# ---------------------------------------------------------------------
# K21 WU1 (F4 T6f, 2026-10-12): integration-point -> node recovery of the
# 2D quad4/quad9 Gauss rules = Professional 02-08-2026 (per-element value,
# arithmetic mean over the adjacent elements). Targets = Professional
# node_dof; whole stress field = Professional to <= 1e-7 of its scale.
#   recovery_examp1_quad4_point_mean1 (sfnet examp1, quad4 1 point): the
#     point value per element; node 19 sigxx -0.2350738 (before -0.256815,
#     volume weighted).
#   recovery_incnav4_quad9_minimal1 (sfnet incnav4, quad9 2x2 Gauss):
#     physical linear fit; nodes 1/3 sigxy 0.05 (before 0.02 / 0.08).
#   recovery_tube1_axis1 (Pro tube1, axisymmetric quad4): on-axis elements
#     (a node at r=0) sum h sig / sum h; node 3 sigxx -159.5805 (before
#     -66.40), node 2 -324.826, node 33 -14.745, node 40 -165.8009.
#   recovery_quad4_distorted_plane1 (distorted plane quad4 2x2): physical
#     linear least squares 1,x,y; node 21 sigxx 2.0014257 (before 1.1829),
#     node 24 5.6815309 (before 5.0343).
# ---------------------------------------------------------------------
run_ok "$T2026" recovery_examp1_quad4_point_mean1
run_ok "$T2026" recovery_incnav4_quad9_minimal1
run_ok "$T2026" recovery_tube1_axis1
run_ok "$T2026" recovery_quad4_distorted_plane1
#   C70 (F4 T6f, 2026-10-13): the other Gauss rules = Professional too
#     (per-element linear least squares 1,x,y(,z) in PHYSICAL coordinates;
#     before: natural-coordinate fits / volume-weighted h):
#     recovery_hex8_distorted1 (hex8 2x2x2) node 11 sigxx 1.3071691
#     (before 1.1216); recovery_hex8_point_mean1 (hex8 1 point: point value,
#     arithmetic mean) node 2 sigxx 2.3571531 (before 2.9277);
#     recovery_hex27_distorted1 (3x3x3) node 41 sigxx 3.6670886 (before
#     3.4134); recovery_tet10_curved1 (Keast 5) node 14 sigxx 2.6187305
#     (before 2.5316); recovery_prism6_linear1 (linear stress, exact) node 1
#     sigxx 0.55 (before 1.0575); recovery_quad9_curved1 (3x3) node 11
#     sigxx 2.4332610 (before 2.2832); recovery_quad9_axis1 (axisymmetric
#     quad9 3x3 / 2x2 at r=0: sum h sig / sum h) node 10 sigxx 0.9518217
#     (before -0.1123).
run_ok "$T2026" recovery_hex8_distorted1
run_ok "$T2026" recovery_hex8_point_mean1
run_ok "$T2026" recovery_hex27_distorted1
run_ok "$T2026" recovery_tet10_curved1
run_ok "$T2026" recovery_prism6_linear1
run_ok "$T2026" recovery_quad9_curved1
run_ok "$T2026" recovery_quad9_axis1
#   C79 (F4 T6f, 2026-10-13): the Lobatto rules = Professional (per
#     element the point value at the node when the points sit at the
#     nodes, the physical linear fit for the quad9/hex27 -minimal corner
#     rules; arithmetic mean over the elements; before: volume weighted h,
#     and no contribution at the quad9/hex27 -minimal non-corner nodes):
#     recovery_lobatto_quad4_mean1 node 8 sigxx 10.4179437 (before
#     11.3589); recovery_lobatto_quad4_axisym_mean1 node 8 sigxx 13.4596104
#     (before 14.4006); recovery_lobatto_hex8_mean1 node 26 sigxx
#     11.6814286 (before 12.5087); recovery_lobatto_quad9_minimal1 node 20
#     sigxx 13.9811919 (before 0); recovery_lobatto_hex27_minimal1 node 100
#     sigxx 16.4575752 (before 0). msf_axisym_lobatto_axis1: the C71
#     on-axis Gauss switch mirrored in post_calcul -materi_stress -force
#     (momx_sig 1/6 = the -gauss run; before 0.5).
run_ok "$T2026" recovery_lobatto_quad4_mean1
run_ok "$T2026" recovery_lobatto_quad4_axisym_mean1
run_ok "$T2026" recovery_lobatto_hex8_mean1
run_ok "$T2026" recovery_lobatto_quad9_minimal1
run_ok "$T2026" recovery_lobatto_hex27_minimal1
run_ok "$T2026" msf_axisym_lobatto_axis1

# ---------------------------------------------------------------------
# K21 WU2 group_materi_stokes (F4 T6f, 2026-10-12): default -yes (manual
# 6.752, = Professional 02-08-2026), read PER ELEMENT GROUP as the manual
# says (owner decision; the Professional reads index 0 for every group:
# DIVERGENCES STOKES-IDX, class A). validation_1 flow.
#   stokes_default1: no stokes record; target = Professional post_point
#     velx 0.1465713 +- 0.003 (C72 residual: GNU 0.144159; before WU2
#     -0.0570386, -no default).
#   stokes_group1: post_point elements in group 1 with -no, group 0 no
#     record -> mixed per-group flow 0.124105 (Pro 0.1465713 fails it).
#   stokes_group2: every element in group 1 with -no, group 0 -yes ->
#     = the one-group -no run -0.0570386 (Pro 0.1465713 fails it).
# ---------------------------------------------------------------------
run_ok "$T2026" stokes_default1
run_ok "$T2026" stokes_group1
run_ok "$T2026" stokes_group2
run_ok "$T2026" convected_stress_recovery1
run_ok "$T2026" validation1_steady1
# F4 T6f C74 (2026-10-13): a breakdown that a direct retry (so.cc) solves
#   is reported with ONE warning line per solve; the full "Error ...
#   Probably you have a wrong input file" block is kept for a breakdown no
#   direct solver can rescue (or solver_bicg_stop -no). Old GNU: the error
#   block for every rescued solve (9 here, 599 on validation1_steady1).
#   F4 T6f C77 (2026-10-15): the 9 breakdowns were false (absolute 1e-16
#   test on tiny vectors); with the scale-relative test the iterative
#   solver converges on every solve: no breakdown warning, no error block.
#   Target = the SuperLU solution.
run_ok "$T2026" bicg_breakdown_lu_retry1
if ! grep -aq "breakdown, system solved with the direct" \
     /tmp/bicg_breakdown_lu_retry1_safe.out 2>/dev/null &&
   ! grep -aq "wrong input file" /tmp/bicg_breakdown_lu_retry1_safe.out 2>/dev/null;
then check_ok "bicg_breakdown_lu_retry1 (no false breakdown, no error block)";
else check_fail "bicg_breakdown_lu_retry1" "breakdown warning or error block printed"; fi
# F4 T6f C77 (2026-10-15): relative stop of the iterative solver. The same
#   flow with the inflow scaled by 1e-5: the old absolute floors (|r|^2 <
#   1e-16) accepted iterates 4.0e-5 relative off the direct solution
#   (velx 9.77555956992e-07); target = SuperLU 9.77517060814e-07 +- 1e-13.
run_ok "$T2026" bicg_relative_stop1

# ---------------------------------------------------------------------
# Sub-sprint materi_stress_force, lote 2 (integracion 2D quad4/quad9):
# verificacion ANALITICA por estatica. Decisiones de diseno (manual-
# developer/post_calcul_materi_stress_force.md): caras extremas = las
# 2 aristas mas perpendiculares a la direccion de espesor t (definida
# por el reference_point); mom con dt = (x-C)*t (distancia en direccion
# de espesor respecto al centroide, CONSISTENTE entre las 2 caras de la
# seccion - el parametro de arista s daria signos opuestos); fuente de
# tension = la nodal (incognita resuelta); cuadratura de arista
# Gauss(2)/Lobatto(3) (el integrando del momento es cuadratico/cubico);
# promedio quad9 (average -yes default) en el plano medio; bandas FE
# documentadas: quad9 mom 1-3% y she +/-30% (condicion de superficie
# libre debil de la formulacion mixta), quad4 1-en-espesor con shear
# locking (mom 0.23x medido) - por eso los tests de flexion usan
# quad9. GOTCHA cazado con ASan: leer coords nodales con MDIM (3) en
# vez de ndim (2) lee 8 bytes FUERA del record NODE (heap OOB) y
# corrompe el solve posterior (campo de tensiones localizado en las
# ultimas columnas) - corregido (i<ndim).
# ---------------------------------------------------------------------

# msf_beam2d (410): mensula 8xquad9, L=8, P=1e-2 en la punta (x=8),
# empotrada en x=0. Secciones primarias (file nodes 0,2,4,...,16 =
# x=0..8): moms = P*(8-x) dentro de 3e-3 (medido 0.0820, 0.0708, 0.0600,
# 0.0500, 0.0400, 0.0300, 0.0200, 0.0099, 0.0003); nors < 2.5e-3 (sin
# axial); shes en la banda documentada [0.5P, 1.5P] (medido 0.0109..0.0133).
MSF_BEAM_OK=1
for i in 0 2 4 6 8 10 12 14 16; do
  x=$((i/2))
  if awk -v i="$i" -v x="$x" '!/^#/ && $1==i { d=$10-0.01*(8-x); if (d<0) d=-d; ok=(d<0.003 && $7>0.005 && $7<0.015 && $4<0.0025) } END{exit !ok}' \
    "$T2014/materi_stress_force.410"; then :; else MSF_BEAM_OK=0; fi
done
# promediado quad9 (LOT 4): node 1 (x=0.5, plano medio del elem 1) =
# la media de las DOS caras del PROPIO elemento (manual 6.906). Con la
# fuente IP (ELEMENT_DOF) las caras de los elementos adyacentes en una
# seccion compartida DIEREN (el campo sigma es discontinuo entre
# elementos) y los nodos primarios de la seccion muestran la MEDIA de
# las dos caras -> el valor promediado ya no es la media EXACTA de los
# valores de salida: banda FE medida |dev| < 0.002 (1.1% en el
# empotramiento, 9.6% en el elemento de la punta cargada).
if awk '!/^#/ && ($1==0||$1==1||$1==2) { v[$1]=$10 } END{ d=v[1]-0.5*(v[0]+v[2]); if (d<0) d=-d; exit !(d<0.002) }' \
   "$T2014/materi_stress_force.410"; then :; else MSF_BEAM_OK=0; fi
if [ "$(grep -vc '^#' "$T2014/materi_stress_force.410")" = "51" ] && \
   [ "$MSF_BEAM_OK" = "1" ]; then
  check_ok "msf_beam2d (mensula 8xquad9: moms=P*(8-x) 1-3%, promediado EXACTO, shes banda FE documentada, nors~0)"
else
  check_fail "msf_beam2d" "moms/she/nors fuera de banda o conteo inesperado"
fi

# msf_beam2d_pure (420): FLEXION PURA 4 puntos (viga apoyada, cargas en
# x=1 y x=3): tramo central x in [1,3] con M = P*1 = 1e-2 cte y V = 0.
# node 4 (x=2): moms ~ 1e-2 (medido 0.01015, 1.5%; SIGNED desde el
# lote 1 del sprint 12: la convencion del Professional, medido -0.01
# con el ref point de este modelo - el check compara MAGNITUDES) y
# shes = 0 EXACTO
# (LOT 4: la fuente IP da la polucion de las caras por elemento, medido
# 2.3e-4 = 2.3% de P - la cancelacion casi exacta del campo nodal
# promediado se pierde; banda < 5e-4 documentada); node 8 (x=4, apoyo):
# moms ~ 0 (medido 0.00044). DECISION: el par de extremos se descarto
# (par nodal no excita flexion en el FE discreto: trabajo nulo sobre el
# modo de flexion - verificado: responde con cizalla pura; quad4 Y
# quad9).
MSF_PURE_OK=1
awk '!/^#/ && ($1==4||$1==8) { v[$1]=$10; if (v[$1]<0) v[$1]=-v[$1]; s[$1]=$7 } END{ d4=v[4]-0.01; if (d4<0) d4=-d4; ok=(d4<5.e-4 && s[4]<5.e-4 && v[8]<1.e-3); exit !ok }' \
  "$T2014/materi_stress_force.420" || MSF_PURE_OK=0
if [ "$(grep -vc '^#' "$T2014/materi_stress_force.420")" = "27" ] && \
   [ "$MSF_PURE_OK" = "1" ]; then
  check_ok "msf_beam2d_pure (flexion pura: M cte 1e-2 en el tramo central, shes=0 EXACTO)"
else
  check_fail "msf_beam2d_pure" "momento/cizalla del tramo central inesperados"
fi

# msf_quad9 (320/321): mensula 4xquad9, L=4, P=1e-2: -all 27 lineas vs
# -primary 15 (los 12 nodos del plano medio son AVERAGED y se omiten);
# secciones primarias x=0,2,4: moms ~ 4e-2, 2e-2, 0 (medido 0.0406,
# 0.0200, 0.0003); el promediado = la media de las caras del PROPIO
# elemento (LOT 4; banda FE medida |dev| < 0.002, ver msf_beam2d).
MSF_Q9_OK=1
if awk '!/^#/ && ($1==0||$1==1||$1==2) { v[$1]=$10 } END{ d=v[1]-0.5*(v[0]+v[2]); if (d<0) d=-d; exit !(d<0.002) }' \
  "$T2014/materi_stress_force.320"; then :; else MSF_Q9_OK=0; fi
if [ "$(grep -vc '^#' "$T2014/materi_stress_force.320")" = "27" ] && \
   [ "$(grep -vc '^#' "$T2014/materi_stress_force.321")" = "15" ] && \
   [ "$MSF_Q9_OK" = "1" ] && \
   awk '!/^#/ && $1==1 {found=1} END{exit found}' "$T2014/materi_stress_force.321"; then
  check_ok "msf_quad9 (-all 27 vs -primary 15: nodos promediados omitidos; moms 1-3%; promedio EXACTO)"
else
  check_fail "msf_quad9" "conteos, promediado o filtro -primary inesperados"
fi

# msf_quad9_noavg (330/331): A/B con average -no: los nodos del plano
# medio NO reciben nada (0) y no se marcan como averaged -> -all ==
# -primary (27 lineas cada uno); los primarios conservan sus valores.
MSF_Q9NA_OK=1
awk '!/^#/ && ($1==1||$1==10||$1==19) { if ($4!=0 || $7!=0 || $10!=0) bad=1 }
     END{exit bad}' "$T2014/materi_stress_force.330" || MSF_Q9NA_OK=0
if [ "$(grep -vc '^#' "$T2014/materi_stress_force.330")" = "27" ] && \
   [ "$(grep -vc '^#' "$T2014/materi_stress_force.331")" = "27" ] && \
   [ "$MSF_Q9NA_OK" = "1" ]; then
  check_ok "msf_quad9_noavg (average -no: nodos medios a 0, -all == -primary 27 lineas)"
else
  check_fail "msf_quad9_noavg" "average -no inesperado (nodos medios no nulos o conteos distintos)"
fi

# msf_nor (340): traccion axial pura N=2e-2 (2 cargas de 1e-2 en x=8,
# nu=0.3): nors = N en todas las secciones (medido 0.0209 en x=0 por la
# perturbacion de Poisson del empotramiento y 0.0200 en el resto),
# shes ~ 0 y moms ~ 0. GOTCHA documentado: con nu=0.0 el solve
# multi-elemento LOCALIZA la deformacion (el dof de tension sale 0 en
# la mitad izquierda - quirk del solver mixto iterativo, no del
# integrador).
MSF_NOR_OK=1
awk '!/^#/ && ($1==0||$1==4||$1==8) { v[$1]=$4; s[$1]=$7; m[$1]=$10 } END{ d0=v[0]-0.02; if (d0<0) d0=-d0; d4=v[4]-0.02; if (d4<0) d4=-d4; d8=v[8]-0.02; if (d8<0) d8=-d8; ok=(d0<2.e-3 && d4<1.e-3 && d8<1.e-3 && s[0]<1.e-6 && m[0]<1.e-3 && s[4]<1.e-6 && m[4]<1.e-3); exit !ok }' \
  "$T2014/materi_stress_force.340" || MSF_NOR_OK=0
if [ "$(grep -vc '^#' "$T2014/materi_stress_force.340")" = "18" ] && \
   [ "$MSF_NOR_OK" = "1" ]; then
  check_ok "msf_nor (axial N=2e-2: nors=N 0-5%, shes~0, moms~0)"
else
  check_fail "msf_nor" "axial inesperado (nors fuera de N o shes/moms no nulos)"
fi

# msf_shear (430): CORTE SIMPLE (borde inferior fijo, superior con
# velx=1e-2 prescrito): campo uniforme EXACTO sigma_xy = G*gamma =
# 384.615*1e-3 = 0.384615 -> she = 0.384615 EXACTO (nor = 0). LOT 5
# (fuerzas internas del elemento): moms = 0.1923 = el momento del PAR de
# reacciones del bloque de corte puro respecto al punto de la seccion
# (las reacciones +-0.3846 de las caras fija/prescrita con brazo 0.5
# respecto al centro de la cara extrema: 0.3846*0.5 = 0.1923 - el
# momento de las fuerzas internas del metodo, no el int sigma_nn*dt del
# manual que es 0). Banda < 0.2 documentada.
MSF_SHEAR_OK=1
awk '!/^#/ && NF==10 { d=$7-0.3846153846; if (d<0) d=-d; if (d>1.e-4 || $4>1.e-6 || $10>0.2) bad=1 } END{exit bad}' \
  "$T2014/materi_stress_force.430" || MSF_SHEAR_OK=0
if [ "$(grep -vc '^#' "$T2014/materi_stress_force.430")" = "4" ] && \
   [ "$MSF_SHEAR_OK" = "1" ]; then
  check_ok "msf_shear (corte simple: she = G*gamma = 0.384615 EXACTO, nor=0, moms = par de reacciones 0.1923)"
else
  check_fail "msf_shear" "corte simple inesperado (she != G*gamma o nor/mom fuera de banda)"
fi

# ---------------------------------------------------------------------
# Sub-sprint materi_stress_force, lote 3 (integracion 3D hex8/hex27).
# GOTCHA MAYOR del GNU (documentado en manual-developer): el solve
# mixto 3D materi_stress + BiCG es DEGENERADO - la matriz v-v del
# hex8/hex27 con la integracion default (Lobatto 2x2x2/3x3x3) da
# A*b ~ 0 para muchas configuraciones de carga (el cantilever con
# carga en punta se detiene en la iteracion 0 con x=0, o diverge con
# hex27: error final 2.16e+13 tras 720 iteraciones) y el ring con
# presion convergia mal (meseta ~1e-12, tensiones contaminadas). Por
# eso los tests 3D usan DEFORMACION PRESCRITA (Dirichlet puro, sin
# solve): campos lineales/cuadraticos exactos para el hex8/hex27, con
# expectativas analiticas EXACTAS. Los .dat documentan cada campo.
# Decisiones de diseno 3D (manual-developer): caras extremas = las 2
# candidatas con |n*t_hat| minima (t_hat = normalizado(centroide -
# reference_point)); el check "not opposing" usa el producto escalar
# de las normales < 0.1 (los +-theta de un anillo de 45 grados no son
# paralelos: dot = -cos(22.5) = -0.92); el cuadratura de cara =
# Gauss(2)xGauss(2) (hex8) / Lobatto(3)xLobatto(3) (hex27) con el
# factor 4 (los dos intervalos [-1,1]); los momentos se miden desde el
# PUNTO MEDIO DE LA CARA (el centroide por esquinas de los elementos
# curvos queda en la cuerda y anade un momento espurio, verificado);
# el hex8 1-en-espesor no puede representar el campo de flexion
# (u_x ~ y^2) -> shear locking (documentado; el test hex8 usa corte
# puro, campo lineal EXACTO).
# ---------------------------------------------------------------------

# msf_sheet3d (500/501): FLEXION PURA PRESCRITA (hex27, 4 elementos,
# 81 nodos): u_y = kappa*(x-0.5)*y, u_x = -kappa*y^2/2, nu=0 ->
# sigma_yy = E*kappa*(x-0.5) EXACTO -> mom1 = E*kappa/12 = 0.0833333
# const en todas las secciones (+-y; SIGNED desde el lote 1 del sprint
# 12: medido -0.0833 con el ref point de este modelo - el check compara
# MAGNITUDES), nor = 0, she = 0. -all 81 vs
# -primary 45 (los 4 planos medios x 9 nodos = 36 promediados).
MSF_S3_OK=1
awk '!/^#/ && NF==17 { m=$13; if (m<0) m=-m; d=m-0.0833333333; if (d<0) d=-d; if (d>1.e-5 || $5>1.e-8 || $9>1.e-8) bad=1 } END{exit bad}' \
  "$T2014/materi_stress_force.500" || MSF_S3_OK=0
if [ "$(grep -vc '^#' "$T2014/materi_stress_force.500")" = "81" ] && \
   [ "$(grep -vc '^#' "$T2014/materi_stress_force.501")" = "45" ] && \
   [ "$MSF_S3_OK" = "1" ]; then
  check_ok "msf_sheet3d (hex27 flexion prescrita: mom1 = E*kappa/12 = 0.0833333 EXACTO en los 81 nodos, nor=she=0, -all 81 vs -primary 45)"
else
  check_fail "msf_sheet3d" "mom1/nor/she o conteos inesperados"
fi

# msf_sheet3d_hex8 (510/511): CORTE PURO PRESCRITO (hex8, campo lineal
# u_x = gamma*y EXACTO para el hex8): sigma_xy = G*gamma -> she =
# G*gamma*t = 500*1e-3 = 0.5 EXACTO en los 20 nodos, nor = 0. LOT 5:
# mom1s = 0.25 = el par de reacciones 3D del corte puro (las reacciones
# +-0.5 de las caras fija/prescrita con brazo 0.5 respecto al centro de
# la cara extrema), banda < 0.3; el valor del manual (int sigma_nn*dt)
# es 0. hex8 sin promedio: -all == -primary (20 lineas cada uno).
MSF_S8_OK=1
awk '!/^#/ && NF==17 { d=$9-0.5; if (d<0) d=-d; if (d>1.e-5 || $5>1.e-8 || $13>0.3) bad=1 } END{exit bad}' \
  "$T2014/materi_stress_force.510" || MSF_S8_OK=0
if [ "$(grep -vc '^#' "$T2014/materi_stress_force.510")" = "20" ] && \
   [ "$(grep -vc '^#' "$T2014/materi_stress_force.511")" = "20" ] && \
   [ "$MSF_S8_OK" = "1" ]; then
  check_ok "msf_sheet3d_hex8 (hex8 corte puro prescrito: she = G*gamma*t = 0.5 EXACTO, nor=mom1=0, -all == -primary)"
else
  check_fail "msf_sheet3d_hex8" "she/nor/mom1 o conteos inesperados"
fi

# msf_hex27_avg (530/531): promedio del plano medio hex27 (manual
# 6.908, average -yes default): placa 2xhex27 con la flexion prescrita:
# mom1 = 0.0833333 en todos; los 18 nodos de los planos medios (9 por
# elemento) reciben la media EXACTA de las 2 caras extremas y se
# marcan averaged -> -all 45 lineas vs -primary 27 (los promediados
# omitidos); los nodos promediados = los 18 que NO estan en las caras
# +-y de los elementos.
MSF_HA_OK=1
awk '!/^#/ && NF==17 { m=$13; if (m<0) m=-m; if (m<0.0833 || m>0.0834) bad=1 } END{exit bad}' \
  "$T2014/materi_stress_force.530" || MSF_HA_OK=0
# los nodos promediados aparecen en -all y faltan en -primary
MSF_HA_AVG=1
for pos in 3 4 5 9 10 11 18 19 20 24 25 26 33 34 35 39 40 41; do
  awk -v p="$pos" '!/^#/ && $1==p {found=1} END{exit !found}' "$T2014/materi_stress_force.530" || MSF_HA_AVG=0
  awk -v p="$pos" '!/^#/ && $1==p {found=1} END{exit found}' "$T2014/materi_stress_force.531" || MSF_HA_AVG=0
done
if [ "$(grep -vc '^#' "$T2014/materi_stress_force.530")" = "45" ] && \
   [ "$(grep -vc '^#' "$T2014/materi_stress_force.531")" = "27" ] && \
   [ "$MSF_HA_OK" = "1" ] && [ "$MSF_HA_AVG" = "1" ]; then
  check_ok "msf_hex27_avg (promedio plano medio hex27: -all 45 vs -primary 27, 18 nodos promediados con la media exacta)"
else
  check_fail "msf_hex27_avg" "conteos, valores o filtro -primary inesperados"
fi

# msf_tunnel3d (520/521): TUNEL - anillo de 8xhex27 con EXPANSION
# RADIAL PROPORCIONAL prescrita u_r = u0*r/R (u0=1e-3, nu=0):
# sigma_theta = E*u0/R = 1.0 const sobre el espesor EXACTO -> nor =
# E*u0*t/R = 0.1 (= p*R con la presion equivalente p = E*u0*t/R^2 =
# 0.1, vaso de presion) en los 144 nodos. LOT 5 (fuerzas internas del
# elemento): nor = 0.0998068, shes = 0.0137271, mom1s = 0.0002028 -
# los valores DISCRETIZADOS FE del campo (la hoop varia en el anillo
# curvo), IDENTICOS a los del Tochnog Professional en el mismo modelo
# (arness §4.2: nor 0.09980686, shes 0.013727, mom1 -2.0e-4) - el
# metodo de equilibrio reproduce el comportamiento del Professional
# (la integracion de campo del L4 daba el valor analitico 0.1 EXACTO
# porque interpolaba el campo prescrito). direction_exclude = eje
# (0 0 1); caras extremas = las 2 +-circunferenciales. -all 144 vs
# -primary 72 (8 planos medios x 9).
MSF_TN_OK=1
awk '!/^#/ && NF==17 { d=$5-0.0998068; if (d<0) d=-d; if (d>2.e-4 || $13>3.e-4 || $9>0.02) bad=1 } END{exit bad}' \
  "$T2014/materi_stress_force.520" || MSF_TN_OK=0
if [ "$(grep -vc '^#' "$T2014/materi_stress_force.520")" = "144" ] && \
   [ "$(grep -vc '^#' "$T2014/materi_stress_force.521")" = "72" ] && \
   [ "$MSF_TN_OK" = "1" ]; then
  check_ok "msf_tunnel3d (hex27 anillo expansion prescrita: nor = E*u0*t/R = p*R = 0.1 EXACTO en los 144 nodos, mom1=she=0, -all 144 vs -primary 72)"
else
  check_fail "msf_tunnel3d" "nor/mom1/she o conteos inesperados"
fi

# msf_errors_3dbothdir: 3D con exclude+include juntos -> error
 msf_error_ok msf_errors_3dbothdir "are mutually exclusive"
 # msf_errors_3drefcount: 3D con 2 grupos y 1 solo reference_point
 msf_error_ok msf_errors_3drefcount "reference_point needs one point"

# ---------------------------------------------------------------------
# Sub-sprint materi_stress_force, lote 4 (fuente de sigma = puntos de
# integracion del elemento + tests nuevos) y lote 5 (fuerzas internas
# del elemento - la estatica del equilibrio: la seccion lee los
# resultantes de f_elem = int B^T sigma dV de los nodos de la cara = la
# estatica del cuerpo libre de las cargas; ver
# manual-developer/post_calcul_materi_stress_force.md LOT 5).
# ---------------------------------------------------------------------

# msf_cant3d_hex27 (600): CANTILEVER 3D hex27 con carga REAL - el test
# que el L3 no pudo hacer (el solve mixto 3D degeneraba; el fix C/D lo
# desbloqueo). 4 hex27 a lo largo de z, seccion 1x1, empotrado en z=0,
# carga total P=1e-2 en -x lumped en las 4 esquinas de la cara z=4.
# ESTATICA: mom1(z) = P*(L-z) (0.04/0.03/0.02/0.01/0 en z=0..4).
# LOT 5 (fuerzas internas del elemento): mom1s = 0.0403, 0.0297, 0.0199,
# 0.00997, 0.00016 (dentro del 1% - la estatica del cuerpo libre de la
# solucion equilibrada; SIGNED desde el lote 1 del sprint 12: medidos
# negativos con el ref point de este modelo - el check compara
# MAGNITUDES) y shes = 0.0097..0.0103 = P dentro del 4% (la
# cizalla contaminada 0.08*P del esquema mixto del L4 desaparece: la
# resultante de las fuerzas internas de los nodos de la cara es la del
# equilibrio); nors ~ 0.
MSF_C3D_OK=1
awk '!/^#/ && NF==17 { m=$13; if (m<0) m=-m; if ($1==1) d=m-0.04; else if ($1==19) d=m-0.03;
       else if ($1==37) d=m-0.02; else if ($1==55) d=m-0.01; else if ($1==73) d=m;
       else next; if (d<0) d=-d; if (d>0.002 || $9>0.011 || $9<0.009 || $5>0.001) bad=1 }
     END{exit bad}' "$T2014/materi_stress_force.600" || MSF_C3D_OK=0
if [ "$(grep -vc '^#' "$T2014/materi_stress_force.600")" = "81" ] && \
   [ "$MSF_C3D_OK" = "1" ]; then
  check_ok "msf_cant3d_hex27 (cantilever 3D hex27 con carga REAL: mom1 = P*(L-z) dentro del 9%, shes/nors en banda FE documentada)"
else
  check_fail "msf_cant3d_hex27" "mom1/she/nors o conteo inesperados"
fi

# msf_axisym (610): AXISIMETRICO - el pendiente del L2. Pared vertical
# (anillo) r in [9.75,10.25], quad4 Z, vely prescrito (Dirichlet, 10
# pasos): eps_zz = -0.001 -> sigma_zz = E*eps = -1.0 EXACTO uniforme.
# La convencion del manual 6.911 (l = 2*PI*r, fuerzas por unidad de
# circunferencia): nor = [int sigma_zz*(2*PI*r) ds]/(2*PI*r) =
# sigma_zz*t = 0.5 EXACTO (el L4 anade el peso 2*PI*r al integrando:
# antes la integral sin pesar daba sigma*t/(2*PI*r) = 0.00796, un
# factor 2*PI*r pequeno - bug documentado; SIGNED desde el lote 1 del
# sprint 12: la pared esta a compresion, nors = -0.5 con la convencion
# tension-positiva del Professional - el check compara MAGNITUDES).
# she = 0, mom ~ 0 (F4 T4 / K20: with the 2x2 rule the face moment is
# 1/960 ~ 0.00104, the analytic value of the 2*PI*r-weighted section;
# the former 1-point rule gave 0).
MSF_AX_OK=1
awk '!/^#/ && NF==10 { n=$4; if (n<0) n=-n; d=n-0.5; if (d<0) d=-d; m=$10; if (m<0) m=-m;
       if (d>0.005 || $7>1.e-6 || m>0.005) bad=1 }
     END{exit bad}' "$T2014/materi_stress_force.610" || MSF_AX_OK=0
if [ "$(grep -vc '^#' "$T2014/materi_stress_force.610")" = "4" ] && \
   [ "$MSF_AX_OK" = "1" ]; then
  check_ok "msf_axisym (anillo axisimetrico: nor = sigma_zz*t = 0.5 EXACTO por unidad de circunferencia, she=0, mom~0)"
else
  check_fail "msf_axisym" "nor/she/mom o conteo inesperados"
fi

# ---------------------------------------------------------------------
# Shear locking del quad4 + integracion reducida selectiva (SRI, Hughes):
# familia qsri. FENOMENO (ver manual-developer/group_element_selective_
# reduced_integration.md): la flexion pura requiere u_y ~ -kappa*x^2/2,
# no representable por el bilineal -> cizalla parastia gamma_xy != 0 ->
# respuesta en flexion DEMASIADO RIGIDA. MEDIDO en qsri_beam2d (A/B):
# sin SRI mom = 0.231*P*(8-x) y she = 0.741*P (reproduce el lock
# documentado del lote MSF 2); con SRI mom = 0.312*P*(8-x) y she =
# 0.65*P (mejoria real; el valor clasico del SRI ~1.0* requiere el
# solve de desplazamiento puro - ver la nota del solver en el manual).
# PATCH TESTS (estados constantes, EXACTOS con SRI): qsri_patch_s she =
# G*gamma = 0.384615384615 EXACTO (12 digitos), qsri_patch_t nor
# identico con/sin SRI (el SRI integra constantes exactamente).
# MODOS: qsri_modes_rigid -> rotacion rigida prescrita da tensiones 0
# EXACTAS (los 3 modos rigidos son modos cero); qsri_modes (hourglass
# prescrito) -> tensiones NO nulas (el modo espurio del quad4 reducido
# TOTAL tiene energia positiva con el SRI: sin modo hourglass).
# ---------------------------------------------------------------------

# qsri_beam2d (440): A/B sin SRI - el LOCK del quad4 1-en-espesor.
# LOT 5 (fuerzas internas del elemento): el momento y la cizalla de
# seccion son la ESTATICA DEL CUERPO LIBRE de la solucion EQUILIBRADA
# (K*u = P, el punto fijo del esquema) - mom en x=0 = 0.08 = P*8 EXACTO
# y she = P = 0.01 EXACTO en TODAS las secciones, con o sin SRI: las
# fuerzas internas de una solucion en equilibrio cumplen la estatica de
# las cargas independientemente de la formulacion del elemento. El lock
# (mom 0.231x/she 0.741x del L2-L4) era una propiedad del CAMPO sigma
# (la tension de flexion del lockeado), ya no visible en las fuerzas de
# seccion (el lock queda en la deflexion y en el campo sigma). El
# discriminador SRI>OFF del A/B se traslada a los valores del campo
# (los targets de sigma del .dat) - aqui se verifican los valores del
# equilibrio.
QSR_BEAM_OFF_OK=1
awk '!/^#/ && $1==0 { d=$10-0.08; if (d<0) d=-d; if (d>5.e-4) bad=1 }
     !/^#/ && $1==4 { d=$7-0.01; if (d<0) d=-d; if (d>2.e-4) bad=1 }
     END{exit bad}' "$T2014/materi_stress_force.440" || QSR_BEAM_OFF_OK=0
# qsri_beam2d_sri (441): A/B con SRI - el punto fijo del esquema
# corregido (fix C/D) es la solucion del ELEMENTO (K_SRI*u = P); con el
# metodo de fuerzas internas la seccion da la estatica del cuerpo libre
# de ESA solucion: mom en x=0 = 0.08 = P*8 EXACTO (los 0.9375x del L4
# eran la integracion SOLO de la tension de flexion del campo sigma) y
# la cizalla interior = P EXACTA (la polucion 8.4*P del sigma_xy crudo
# de los puntos de Gauss desaparece: la resultante de las fuerzas
# internas del nodo de la cara es la del equilibrio, no el integral del
# campo). La correccion SRI de la fuerza interna (la cizalla reducida
# de 1 punto del feedback, ver calcul_force.cc LOT 5) es necesaria: sin
# ella la fuerza interna leeria K_full*u_sri (el lock), no K_SRI*u_sri.
QSR_BEAM_SRI_OK=1
awk '!/^#/ && $1==0 { d=$10-0.08; if (d<0) d=-d; if (d>5.e-4) bad=1; if ($10<0.05) bad=1 }
     !/^#/ && $1==4 { d=$7-0.01; if (d<0) d=-d; if (d>2.e-3) bad=1 }
     END{exit bad}' "$T2014/materi_stress_force.441" || QSR_BEAM_SRI_OK=0
if [ "$QSR_BEAM_OFF_OK" = "1" ] && [ "$QSR_BEAM_SRI_OK" = "1" ]; then
  check_ok "qsri_beam2d (LOT 5: estatica del cuerpo libre EXACTA mom=P*8 she=P con y sin SRI - el lock queda en el campo sigma)"
else
  check_fail "qsri_beam2d" "estatica del equilibrio inesperada"
fi

# qsri_patch_s (442) / qsri_patch_s_off (443): corte simple prescrito
# gamma=1e-3 -> she = G*gamma = 0.384615384615 EXACTO en los 4 nodos,
# nor = 0; moms = 0.1923 = el par de reacciones (igual que msf_shear,
# LOT 5); IDENTICO con y sin SRI (el SRI integra constantes
# exactamente) -> los archivos deben ser byte-identicos salvo cabecera
QSR_PATCH_S_OK=1
awk '!/^#/ { d=$7-0.384615384615; if (d<0) d=-d; if (d>1.e-6 || $4>1.e-8 || $10>0.2) bad=1 }
     END{exit bad}' "$T2014/materi_stress_force.442" || QSR_PATCH_S_OK=0
awk '!/^#/ { d=$7-0.384615384615; if (d<0) d=-d; if (d>1.e-6 || $4>1.e-8 || $10>0.2) bad=1 }
     END{exit bad}' "$T2014/materi_stress_force.443" || QSR_PATCH_S_OK=0
if [ "$QSR_PATCH_S_OK" = "1" ]; then
  check_ok "qsri_patch_s (she = G*gamma = 0.384615384615 EXACTO, nor=0, moms = par de reacciones 0.1923; identico con/sin SRI)"
else
  check_fail "qsri_patch_s" "el estado constante no es exacto o difiere entre SRI/no-SRI"
fi

# qsri_patch_t (444) / qsri_patch_t_off (445): traccion uniaxial
# prescrita eps=1e-3 -> nor = 1.025641 (valor del esquema mixto,
# IDENTICO con y sin SRI: el SRI no degrada el estado constante)
QSR_PATCH_T_OK=1
awk '!/^#/ { d=$4-1.025641; if (d<0) d=-d; if (d>1.e-5) bad=1 }
     END{exit bad}' "$T2014/materi_stress_force.444" || QSR_PATCH_T_OK=0
awk '!/^#/ { d=$4-1.025641; if (d<0) d=-d; if (d>1.e-5) bad=1 }
     END{exit bad}' "$T2014/materi_stress_force.445" || QSR_PATCH_T_OK=0
if [ "$QSR_PATCH_T_OK" = "1" ]; then
  check_ok "qsri_patch_t (nor identico con/sin SRI = 1.025641: estado constante no degradado)"
else
  check_fail "qsri_patch_t" "el estado constante difiere entre SRI/no-SRI"
fi

# qsri_modes (448): modo hourglass prescrito u_y = a*x*y -> con el SRI
# las tensiones son NO nulas (el modo tiene energia positiva: NO hay
# modo espurio de energia cero). she en [0.2, 0.35] y nor en [0.15, 0.3]
QSR_MODES_OK=1
awk '!/^#/ { if ($7<0.2 || $7>0.35 || $4<0.15 || $4>0.3) bad=1 }
     END{exit bad}' "$T2014/materi_stress_force.448" || QSR_MODES_OK=0
# qsri_modes_rigid (449): rotacion rigida prescrita -> tensiones 0
# EXACTAS (los 3 modos rigidos son modos cero exactos de la matriz SRI)
QSR_MODES_RIGID_OK=1
awk '!/^#/ { if ($4!=0 || $7!=0 || $10!=0) bad=1 }
     END{exit bad}' "$T2014/materi_stress_force.449" || QSR_MODES_RIGID_OK=0
if [ "$QSR_MODES_OK" = "1" ] && [ "$QSR_MODES_RIGID_OK" = "1" ]; then
  check_ok "qsri_modes (hourglass con energia positiva; rotacion rigida con tension 0 EXACTA: 3 modos cero, sin hourglass)"
else
  check_fail "qsri_modes" "modos espurios o rigidos incorrectos"
fi

# ---------------------------------------------------------------------
# SRI hex8 (3D) - extension de la integracion reducida selectiva al
# trilineal (familia qsri3d, archivos 700-711). HALLAZGO documentado
# (2026-08-29, manual-developer/group_element_selective_reduced_
# integration.md): el SRI hex8 de cizalla sola (Hughes) retiene modos
# de energia cero - 9 en el elemento aislado (6 rigidos + 3 twist) y
# modos de warping de seccion en malla (u_y = A(x)*(2z-1),
# u_z = A(x)*(2y-1): deformaciones normales cero, cizallas que se
# anulan en el centroide de seccion). La mensula empotrada 3D
# (qsri3d_beam) sigue siendo resoluble (la carga es ortogonal al
# null-space y el CG honesto converge) y el SRI hex8 recupera ~0.90x
# de la deflexion de Euler-Bernoulli frente al 0.221x del hex8 full
# (el lock documentado del L3) - medido via el target de velocidad de
# la punta en el .dat. Los patch tests (estados constantes) son
# EXACTOS e IDENTICOS con y sin SRI (el SRI integra constantes
# exactamente); la rotacion rigida da tensiones 0 EXACTAS (los 6
# modos rigidos son modos cero); el modo twist prescrito es invisible
# a la seccion SRI (la regla reducida de 1 punto no lo ve - el modo de
# energia cero, documentado).
# ---------------------------------------------------------------------

# qsri3d_patch_s (700/701) / _off (702/703): corte puro prescrito
# (placa de msf_sheet3d_hex8, 4 hex8, u_x = gamma*y, nu=0):
# she = G*gamma*t = 0.5 EXACTO, nor = 0, mom1 = 0.25 (el par de
# reacciones); IDENTICO con y sin SRI -> los 3 archivos deben dar los
# mismos valores.
QSR3D_PATCH_S_OK=1
for f in 700 702; do
  awk '!/^#/ { d=$9-0.5; if (d<0) d=-d; if (d>1.e-6 || $5>1.e-8) bad=1 }
       END{exit bad}' "$T2014/materi_stress_force.$f" || QSR3D_PATCH_S_OK=0
done
if [ "$QSR3D_PATCH_S_OK" = "1" ]; then
  check_ok "qsri3d_patch_s (she = G*gamma*t = 0.5 EXACTO, nor = 0, mom1 = par 0.25; identico con/sin SRI: el SRI integra constantes)"
else
  check_fail "qsri3d_patch_s" "el corte puro 3D no es exacto o difiere entre SRI/no-SRI"
fi

# qsri3d_patch_t (704/705) / _off (706/707): traccion confinada
# prescrita u_x = eps*x (eps = 1e-3): nor = sigma_xx = 1.346154
# EXACTO, identico con/sin SRI (el estado constante no se degrada)
QSR3D_PATCH_T_OK=1
for f in 704 706; do
  awk '!/^#/ { d=$5-1.346154; if (d<0) d=-d; if (d>1.e-5 || $9>1.e-8) bad=1 }
       END{exit bad}' "$T2014/materi_stress_force.$f" || QSR3D_PATCH_T_OK=0
done
if [ "$QSR3D_PATCH_T_OK" = "1" ]; then
  check_ok "qsri3d_patch_t (nor = 1.346154 EXACTO, identico con/sin SRI: estado constante no degradado)"
else
  check_fail "qsri3d_patch_t" "la traccion 3D difiere entre SRI/no-SRI"
fi

# qsri3d_modes_rigid (708): rotacion rigida prescrita -> seccion ~0
# (ruido numerico 1e-17; los 6 modos rigidos son modos cero exactos
# del SRI hex8)
QSR3D_RIGID_OK=1
awk '!/^#/ { if ($5>1.e-12 || $9>1.e-12 || $13>1.e-12) bad=1 } END{exit bad}' \
  "$T2014/materi_stress_force.708" || QSR3D_RIGID_OK=0
if [ "$QSR3D_RIGID_OK" = "1" ]; then
  check_ok "qsri3d_modes_rigid (rotacion rigida 3D -> seccion 0 EXACTA: 6 modos rigidos)"
else
  check_fail "qsri3d_modes_rigid" "modos rigidos no exactos"
fi

# qsri3d_modes (709): modo twist prescrito u_x = a*(y-0.5)*(z-0.5) ->
# la seccion SRI lee ~0 (el modo tiene energia CERO bajo la regla
# reducida: deformaciones normales cero y cizallas que se anulan en el
# centroide de seccion; el analisis de eigenvalues documenta los 9
# modos cero del SRI hex8 aislado). La ley constitutiva COMPLETA si ve
# las cizallas en los puntos 2x2x2 (los dofs sigma de salida son no
# nulos) - el modo no es invisible al campo sigma, solo a la rigidez
# SRI.
QSR3D_TWIST_OK=1
awk '!/^#/ { if ($5>1.e-8 || $9>1.e-8 || $13>1.e-8) bad=1 } END{exit bad}' \
  "$T2014/materi_stress_force.709" || QSR3D_TWIST_OK=0
if [ "$QSR3D_TWIST_OK" = "1" ]; then
  check_ok "qsri3d_modes (modo twist prescrito -> seccion ~0: el modo de energia cero del SRI hex8, documentado)"
else
  check_fail "qsri3d_modes" "el modo twist no es invisible a la regla reducida"
fi

# qsri3d_beam (710) / _off (711): MENSULA 3D 1-en-espesor hex8
# (L=8, seccion 1x1, P=1e-2 en la punta). El discriminador SRI>OFF es
# la DEFLEXION de la punta (targets de velocidad en los .dat: SRI
# 0.18379 vs OFF 0.04528 = 0.897x vs 0.221x de la Euler-Bernoulli -
# el lock del hex8 full y la recuperacion del SRI). Las secciones
# 710/711 leen las fuerzas internas del elemento (LOT 5): con el fix
# del ELEMENT_DOF 3D (2026-08-29 - el write del sigma constitutivo al
# slot nder-correcto) el campo sigma de los IPs ya NO es cero y el
# fallback a la nodal quedo obsoleto: la seccion lee el estado en
# equilibrio -> mom1s = P*L = 0.08 EXACTO para AMBOS (la estatica del
# cuerpo libre es independiente de la formulacion del elemento; el
# discriminador SRI es la deflexion, no el momento de seccion).
QSR3D_BEAM_OK=1
M710=$(awk '!/^#/ {print $13; exit}' "$T2014/materi_stress_force.710" 2>/dev/null)
M711=$(awk '!/^#/ {print $13; exit}' "$T2014/materi_stress_force.711" 2>/dev/null)
if [ -n "$M710" ] && [ -n "$M711" ] && \
   awk -v a="$M710" -v b="$M711" \
   'BEGIN{ d1=a-0.08; d2=b-0.08;
          exit !(d1<0.005 && d1>-0.005 && d2<0.005 && d2>-0.005) }'; then
  check_ok "qsri3d_beam (SRI: deflexion 0.897x EB vs OFF 0.221x - targets .dat; mom1s 710 ~ 711 = P*L = 0.08: la estatica del cuerpo libre del estado equilibrado)"
else
  check_fail "qsri3d_beam" "momens de seccion o targets inesperados (710=$M710 711=$M711)"
fi

# ---------------------------------------------------------------------
# Sprint 12 lote 3: familia solver_* global + timestep_*/tochnog_version/
# volume_factor_x/zip (manual Professional 6.1047-6.1056, 6.1089-6.1091,
# 6.1094-6.1095). El zip corre en directorio AISLADO: su glob gzipearia
# todos los *.dbs del directorio compartido.
# ---------------------------------------------------------------------
TSLV_OK=1
# tslv_override: el banner de SuperLU prueba que el record GLOBAL
# solver piso al control_solver -matrix_iterative_bicg
if rg -q "Solution Found" /tmp/tslv_override_safe.out 2>/dev/null; then
  check_ok "tslv_override (solver -matrix_superlu GLOBAL pisa a control_solver bicg: banner SuperLU en stdout)"
else
  check_fail "tslv_override" "sin banner de SuperLU - el override global no tomo efecto"; TSLV_OK=0
fi
# tslv_bsym: el warning de solver_matrix_symmetric -yes sobre sistema
# medido NO simetrico (el beam) - el record se lee y se consume
if rg -q "solver_matrix_symmetric -yes but the measured" /tmp/tslv_bsym_safe.out 2>/dev/null; then
  check_ok "tslv_bsym (warning de -yes sobre sistema NO simetrico medido: el record se consume)"
else
  check_fail "tslv_bsym" "sin warning de matrix_symmetric"; TSLV_OK=0
fi
# tslv_vfx: la estatica de serie con espesor 2 en la mitad derecha
if v9=$(awk '$1=="node_dof"&&$2==9{print $3; exit}' "$T2014/tslv_vfx.dbs" 2>/dev/null) && \
   awk -v v="$v9" 'BEGIN{exit !(v>1.0e-3 && v<1.3e-3)}'; then
  check_ok "tslv_vfx (volume_factor_x: velx punta 1.15e-3 = N*4/(E*1)+N*4/(E*2), A/B 1.59e-3 sin el record)"
else
  check_fail "tslv_vfx" "velx punta fuera de banda (esperada ~1.15e-3)"; TSLV_OK=0
fi
# tslv_ver: el record tochnog_version con el ano de build (target del .dat)
if rg -aq "^tochnog_version" "$T2014/tslv_ver.dbs" 2>/dev/null; then
  check_ok "tslv_ver (tochnog_version presente en el .dbs con la fecha de build)"
else
  check_fail "tslv_ver" "sin record tochnog_version en el .dbs"; TSLV_OK=0
fi
# tslv_zip: corrida AISLADA - el .dbs.gz aparece y no hay warning
ZDIR="/tmp/tslv_zip_isolated"
rm -rf "$ZDIR" && mkdir -p "$ZDIR" && cp "$T2014/tslv_zip.dat" "$ZDIR/" && \
  ( cd "$ZDIR" && ulimit -v "$TN_MEMLIMIT_KB" && timeout 120 "$REPO_DIR/build/$TN_BIN" tslv_zip.dat > run.out 2>&1 )
if [ -f "$ZDIR/tslv_zip.dbs.gz" ] && ! rg -q "gzipping.*failed" "$ZDIR/run.out" 2>/dev/null; then
  check_ok "tslv_zip (zip -yes: .dbs.gz creado al final, sin warning)"
else
  check_fail "tslv_zip" "sin .dbs.gz o con warning de gzip"; TSLV_OK=0
fi
rm -rf "$ZDIR"
if [ "$TSLV_OK" = "1" ]; then
  echo "==> Sprint 12 lote 3 (solver_* global + timestep_* + tochnog_version + volume_factor_x + zip): TODAS OK"
fi

# ---------------------------------------------------------------------
# Sprint 13 lotes 1-3: support_edge_normal (Winkler distribuido, manual
# Professional 6.1067-6.1078) + gforce7 (arness canary). Chequeos de los
# records de salida de los .dbs.
# F0 T4 (2026-10-05): this is the ONE set of tsup checks (the former
# reduced-mode set, which passes on the regenerated inputs). The old
# full-mode duplicate was removed: its targets described the pre-incident
# inputs (e.g. tsup_tcap fx=-0.005) and its `rg -a "<rec>" 1 0 <file>`
# calls never matched. tsup_gap / tsup_tcap now really check (the
# reduced versions exited on the first record / had no END exit).
# ---------------------------------------------------------------------
tsup_check() {
  local t="$1" desc="$2" prog="$3"
  if awk "$prog" "$T2014/$t.dbs" 2>/dev/null; then check_ok "$t ($desc)";
  else check_fail "$t" "$desc: fuera de banda"; fi
}
tsup_check tsup_winkler "Winkler prescrito: fuerza nodal (-0.05, 0.5) EXACTA" \
  '$1=="node_support_edge_normal_force" && ($2==1||$2==2) { d=$4-0.5; if (d<0) d=-d; if (d>1.e-4) bad=1; d=$3+0.05; if (d<0) d=-d; if (d>1.e-4) bad=1 } END { exit bad }'
tsup_check tsup_solve "columna sobre muelle: uy = -3e-3" \
  '$1=="node_dof" && $2==3 { d=$10+0.003; if (d<0) d=-d; exit !(d<3.2e-5) }'
tsup_check tsup_3d "cara hex8: fuerza nodal z = 0.025 por esquina" \
  '$1=="node_support_edge_normal_force" && NF==6 { d=$5-0.025; if (d<0) d=-d; if (d>1.e-4) bad=1 } END { exit bad }'
tsup_check tsup_gap "tension gap -yes: fx=0 y status=1 en el nodo 1" \
  '$1=="node_support_edge_normal_force" && $2==1 { d=$3; if (d<0) d=-d; f=(d<1.e-4) } $1=="node_support_edge_normal_plasti_tension_status" && $2==1 { st=($3==1) } END { exit !(f && st) }'
tsup_check tsup_tcap "tension_double cap: fx=0 en el nodo 1" \
  '$1=="node_support_edge_normal_force" && $2==1 { n++; d=$3; if (d<0) d=-d; if (d>1.e-4) bad=1 } END { exit (bad || !n) }'
tsup_check tsup_fric "friccion sin clip: fx=0.005 en el nodo 1" \
  '$1=="node_support_edge_normal_force" && $2==1 { d=$3-0.005; if (d<0) d=-d; exit !(d<1.e-4) }'
tsup_check tsup_residual "factor=0 default: fx=0.005 en el nodo 1" \
  '$1=="node_support_edge_normal_force" && $2==1 { d=$3-0.005; if (d<0) d=-d; exit !(d<1.e-4) }'
tsup_check tsup_damp "damping+factor: fx=-0.0075 en el nodo 1" \
  '$1=="node_support_edge_normal_force" && $2==1 { d=$3+0.0075; if (d<0) d=-d; exit !(d<1.e-4) }'
tsup_check tsup_auto "damping automatico del grupo: fy=2.73607 en el nodo 1" \
  '$1=="node_support_edge_normal_force" && $2==1 { d=$4-2.73607; if (d<0) d=-d; exit !(d<1.e-3) }'
tsup_check tsup_init "force_initial + time factor: fy=0.1 en el nodo 1" \
  '$1=="node_support_edge_normal_force" && $2==1 { d=$4-0.1; if (d<0) d=-d; exit !(d<1.e-4) }'
tsup_check tsup_dens "density: fy=0.525 en los nodos 1 y 2" \
  '$1=="node_support_edge_normal_force" && ($2==1||$2==2) { d=$4-0.525; if (d<0) d=-d; exit !(d<1.e-4) }'
tsup_check gforce7 "arness canary: nory=-12.34" \
  '$1=="node_dof_calcul" && $2==2 { d=$4+12.34; if (d<0) d=-d; exit !(d<0.05) }'

# Gate de runs (F0 T4): un run del bucle con rc!=0 (crash, timeout, ASan
# abort, target no alcanzado) es un fallo aunque sus file-checks pasen.
if [ "$HIPO_OK" != "$HIPO_TOTAL" ]; then
  echo "!! runs fallidos: $((HIPO_TOTAL-HIPO_OK))/$HIPO_TOTAL:$RUN_FAILED (detalle en /tmp/<test>_safe.out)"
  CHECK_FAIL=1
fi

echo "==> Checks: $N_CHECK_OK OK / $N_CHECK_FAIL FALLO; runs: $HIPO_OK/$HIPO_TOTAL OK"
if [ "$CHECK_FAIL" = "1" ]; then
  echo "==> ALGUNAS VERIFICACIONES DE ARCHIVOS FALLARON"
  exit 1
else
  echo "==> Verificacion de archivos de salida (Sprint 11 lotes 1-6 + sub-sprint materi_stress_force lotes 1-5 + familia qsri + familia qsri3d + lote 3 sprint 12): TODAS OK"
fi

echo "==> Log de compilacion completo en /tmp/tn_build_safe.log"
