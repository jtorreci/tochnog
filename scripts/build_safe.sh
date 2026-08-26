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
# Sprint 11 lote 1 (control_print_*): los archivos de salida que se
# APPENDEN entre ejecuciones se limpian antes del bucle para que las
# verificaciones de la seccion posterior partan de archivos frescos.
rm -f validation-suite/test-2014/dof.20 validation-suite/test-2014/dof.21 \
      validation-suite/test-2014/coord.20 validation-suite/test-2014/coord.21 \
      validation-suite/test-2014/disx1.his validation-suite/test-2014/disx2.his
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
#   desviadora a sigxx -1.225 y kappa crece a 0.000293; el estado final
#   cumple f ~ 0: 0.687/(1000*(1-0.687/38.49)) - 2*0.687/3000 -
#   0.000293 = -0.000061). mhardsoil_plast (sigxx=-4): kappa 0.0006 =
#   2x gp0_off -> el endurecimiento ESCALA con q inicial. GOTCHAS: el
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
         mcap2 mcap_legacy mcrunch mcrunch_low mvoid mvoid_low mpower mshf mshf_nof mk0 mk0_off \
         myoung6 myoung6_e2 myoung6_e3 myoung6_apply msph msph_flat \
         mcap1 mcap1_elast mcap1_comb \
         mhardsoil_elast mhardsoil_elast2 mhardsoil_unload mhardsoil_unload_flat \
         mhardsoil_plast mhardsoil_plast_elast mhardsoil_gp0 mhardsoil_gp0_off \
         mstrain_cap mstrain_cap_elast mstrain_compression mstrain_compression_elast \
         mstrain_diprisco mstrain_diprisco_elast \
         mstrain_druckprag mstrain_druckprag_elast \
         mdiprisco_hist mc_pressure_min mc_pressure_min_off mrepeat_save \
         dbmeth partialname meshdoff dofrhside elmethod hreltime numit dofid_no; do
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
echo "==> Resumen: $HIPO_OK/$HIPO_TOTAL runs OK (13 tests: 12 preexistentes + familia iface_mc en 10 runs + familia 3D en 5 runs + familia generate_interface en 6 runs + familia materi_direct en 5 runs + materi_displacement_relative en 2 runs + slide/reset_value en 2 runs + cda_arith/copy/activate en 3 runs + cdist_normal/corr/clamp en 3 runs + cd_method/cd_geom en 2 runs + gravity/settlement en 4 runs + contact en 1 run + contact_block/ctrl_apply/heatgen en 3 runs + groundflow_consolidate_off en 1 run + groundflow_vangenuchten/groundflow_nonsaturated_off en 2 runs + groundflow_total_pressure_tension/groundflow_interface en 2 runs + groundflow_flux_edge en 1 run + groundflow_phreatic_multiple en 1 run + groundflow_seepage en 1 run + groundflow_pressure_atm/_def en 2 runs + groundflow_total_pressure_limit/_dry en 2 runs + condif_heat_edge/vol/vol2 en 3 runs + condif_convec/rad/convec_el en 3 runs + aeg_node/aeg_seq/bt_factor en 3 runs + iface_condif/expansion/tangref en 3 runs + node_force_inertia/slide/pressure en 3 runs + creset_geom/iface en 2 runs + fedge_alias/restrict, fvol_elem y cmat_gate en 4 runs + fproj_tunnel en 1 run + dsmall/dignore en 2 runs + mdirect_comp/gate en 2 runs + mdp_shear/mfactor en 2 runs + mmc_tension en 1 run + mmchs_soft en 1 run + mcap2/mcap_legacy en 2 runs + mcrunch/mcrunch_low en 2 runs + mvoid/mvoid_low en 2 runs + mpower en 1 run + mshf/mshf_nof en 2 runs + mk0/mk0_off en 2 runs + myoung6/myoung6_e2/myoung6_e3/myoung6_apply en 4 runs + msph/msph_flat en 2 runs + mcap1/mcap1_elast/mcap1_comb en 3 runs + mhardsoil_elast/elast2/unload/unload_flat/plast/plast_elast/gp0/gp0_off en 8 runs + mstrain_cap/_elast, mstrain_compression/_elast, mstrain_diprisco/_elast, mstrain_druckprag/_elast en 8 runs + mdiprisco_hist en 1 run + mc_pressure_min/_off en 2 runs (Sprint 10 lote 9) + mrepeat_save en 1 run (Sprint 10 lote 10) + dbmeth/partialname/meshdoff/dofrhside/elmethod/hreltime/numit/dofid_no en 8 runs (Sprint 11 lote 1))."

# ---------------------------------------------------------------------
# Sprint 11 lote 1: verificacion de ARCHIVOS y STDOUT de los 8 keywords
# control_print_* (los targets de los .dat verifican el modelo; aqui se
# verifican los archivos generados y el stdout capturado en /tmp).
# ---------------------------------------------------------------------
T2014="$REPO_DIR/validation-suite/test-2014"
CHECK_FAIL=0
check_ok()   { echo "    $1: OK"; }
check_fail() { echo "    $1: FALLO ($2)"; CHECK_FAIL=1; }

# dof1: control_print_dof_id DEFAULT -yes -> dof.20 con 4 columnas (x y dof node)
if awk 'NF!=4{exit 1}' "$T2014/dof.20" 2>/dev/null; then
  check_ok "dof.20 (control_print_dof_id default -yes: lineas x y dof node)"
else
  check_fail "dof.20" "esperaba 4 columnas por linea (x y dof node)"
fi

# dofid_no: control_print_dof_id -no -> dof.21 con 3 columnas (x y dof)
if awk 'NF!=3{exit 1}' "$T2014/dof.21" 2>/dev/null; then
  check_ok "dof.21 (control_print_dof_id -no: lineas x y dof)"
else
  check_fail "dof.21" "esperaba 3 columnas por linea (x y dof)"
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

if [ "$CHECK_FAIL" = "1" ]; then
  echo "==> ALGUNAS VERIFICACIONES DE ARCHIVOS FALLARON"
  exit 1
else
  echo "==> Verificacion de archivos de salida (Sprint 11 lote 1): TODAS OK"
fi

echo "==> Log de compilacion completo en /tmp/tn_build_safe.log"
