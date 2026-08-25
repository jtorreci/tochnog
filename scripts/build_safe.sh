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

MAKE_FLAGS=( "SYS_FILE=sysposix" "OBJ=o" "BCPP=" "VCPP="
  "COMPILER_C=gcc" "COMPILER_CPP=g++"
  "COMPILER_FLAGS=-c -O1 -Wall -D_REENTRANT $SUPERLU_INC $SQLITE_INC"
  "LINK_FLAGS_BEFORE=" )

LINK_FLAGS_AFTER="$NUMLIB_LINK -Wl,--start-group $SUPERLU_A -l:liblapack.so.3 -l:libblas.so.3 $NUMLIB_FORTRAN -Wl,--end-group $SQLITE_LIB -lm -lpthread -o build/tochnog"

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
# + groundflow_pressure_atm/_def (groundflow_pressure_atmospheric, keyword
#   heredado del GNU: cap de la presion estatica en phreatic_coord. Con
#   pa=0.5 la succion +1 se mantiene en 0.5 y la compresion -1 pasa intacta;
#   con default 0 la succion se anula -> 0).
# + groundflow_total_pressure_limit/_dry (manual Professional 6.588: cap del
#   pres resuelto tras el solve; nodos Dirichlet exentos. Con limit 0.5 los
#   nodos libres quedan en 0.5 (sin limit ~1.67); con limit 0 y pres 0 el
#   elemento esta seco y se salta el termino de consolidacion -> pres 0).
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
#   node_inertia calculado m*a+m*g PUT en VERSION_NEW: -19.83 ~= m*g=-20
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
#   sigxx=-1.7143/sigyy=-2.5714; el codigo incremental (C evaluada con la
#   tension del paso previo) aterriza en nu~0.43, p~2.06 -> sigxx -1.4389
#   y sigyy -3.2924 (84% del analitico; ventana que discrimina la base
#   nu=0.3 -> -0.5769). control_timestep_iterations 8).
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
for t in hypo1 hypo2 hypo3 hypo4 smooth1 dof1 mlx1 vtk_dof1 gen1 genbeam1          reset1 cda1 cda_arith cda_copy cda_activate cdist_normal cdist_corr cdist_clamp cd_method cd_geom \
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
         contact_block contact_ctrl_apply contact_heatgen \
         groundflow_consolidate_off groundflow_vangenuchten groundflow_nonsaturated_off \
         groundflow_total_pressure_tension groundflow_interface groundflow_flux_edge \
         groundflow_phreatic_multiple groundflow_seepage \
         groundflow_pressure_atm groundflow_pressure_atm_def \
         groundflow_total_pressure_limit groundflow_total_pressure_limit_dry \
         condif_heat_edge condif_heat_vol condif_heat_vol2 \
         condif_convec condif_rad condif_convec_el \
         aeg_node aeg_seq bt_factor \
         iface_condif iface_expansion iface_tangref \
         node_force_inertia node_slide node_pressure \
         creset_geom creset_iface \
         fedge_alias fedge_restrict fvol_elem cmat_gate fproj_tunnel \
         dsmall dignore          mdirect_comp mdirect_gate mdp_shear mfactor mmc_tension mmchs_soft \
         mcap2 mcap_legacy mcrunch mcrunch_low mvoid mvoid_low mpower mshf mshf_nof mk0 mk0_off; do
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
echo "==> Resumen: $HIPO_OK/$HIPO_TOTAL runs OK (13 tests: 12 preexistentes + familia iface_mc en 10 runs + familia 3D en 5 runs + familia generate_interface en 6 runs + familia materi_direct en 5 runs + materi_displacement_relative en 2 runs + slide/reset_value en 2 runs + cda_arith/copy/activate en 3 runs + cdist_normal/corr/clamp en 3 runs + cd_method/cd_geom en 2 runs + gravity/settlement en 4 runs + contact en 1 run + contact_block/ctrl_apply/heatgen en 3 runs + groundflow_consolidate_off en 1 run + groundflow_vangenuchten/groundflow_nonsaturated_off en 2 runs + groundflow_total_pressure_tension/groundflow_interface en 2 runs + groundflow_flux_edge en 1 run + groundflow_phreatic_multiple en 1 run + groundflow_seepage en 1 run + groundflow_pressure_atm/_def en 2 runs + groundflow_total_pressure_limit/_dry en 2 runs + condif_heat_edge/vol/vol2 en 3 runs + condif_convec/rad/convec_el en 3 runs + aeg_node/aeg_seq/bt_factor en 3 runs + iface_condif/expansion/tangref en 3 runs + node_force_inertia/slide/pressure en 3 runs + creset_geom/iface en 2 runs + fedge_alias/restrict, fvol_elem y cmat_gate en 4 runs + fproj_tunnel en 1 run + dsmall/dignore en 2 runs + mdirect_comp/gate en 2 runs + mdp_shear/mfactor en 2 runs + mmc_tension en 1 run + mmchs_soft en 1 run + mcap2/mcap_legacy en 2 runs + mcrunch/mcrunch_low en 2 runs + mvoid/mvoid_low en 2 runs + mpower en 1 run + mshf/mshf_nof en 2 runs + mk0/mk0_off en 2 runs)."
echo "==> Log de compilacion completo en /tmp/tn_build_safe.log"
