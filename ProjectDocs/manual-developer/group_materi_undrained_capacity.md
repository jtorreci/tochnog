# group_materi_undrained_capacity — implementation

Undrained groundwater analysis on the element level (manual Professional
6.760 + theory 2.2.7). Family implemented 2026-09-05 to unblock the corpus
tests ground17, undrained1 and undrained2.

## Keywords

- `group_materi_undrained_capacity` — DOUBLE_PRECISION, 1 value,
  data_class MATERI, data_required GROUP_TYPE (database.cc).
- `element_intpnt_materi_undrained_pressure` — DOUBLE_PRECISION, one slot
  per integration point (data_length npointmax, fixed_length 0 =
  variable record length, **version_all = 1** like ELEMENT_DOF so the
  step-end version copy promotes the converged per-step value to
  VERSION_NORMAL). data_class ELEMENT.
- `element_intpnt_materi_undrained_pressure_average` — registered (1 value,
  not yet filled; cosmetic record of the Professional .dbs).
- `control_materi_undrained_apply` — already registered (parse-only since an
  earlier batch); now consumed in materi.cc.
- `young_apparent` / `poisson_apparent` — INTEGER name entries of the
  post_calcul record (operator values); naming + value computation in
  calcul.cc.
- alias `-tpres` -> GROUNDFLOW_PRESSURE in db_number() (database.cc): the
  corpus spelling of `-topres` (undrained1 uses `-tpres`, ground13 uses
  `-topres`; the manual writes `-topres`).

## Physics (verified against the Professional binary 25-10-2023)

Storage equation without permeability, solved element-wise:
`C * p_dot = div(v)`. Per integration point `p_u += dvol/C` per step. The
total stress of the momentum equilibrium is `sigma_eff + (p_fixed +
p_u)*I` and the momentum tangent gains the volumetric term
`dt/C * (div w)(div v)` (the linear system must carry the undrained
stiffness: a fixed-point that lags p_u needs many iterations; measured on
undrained2: `eps = F/(E + 1/C)`).

Measured on the corpus models:
- undrained2 (1D bar, E=3, C=1, F=1): sigma_dof = E*eps = -0.75,
  p_u = -0.25, eps = -0.25 EXACT (the Professional .dbs).
- undrained1 (2D quad4 + groundflow bounded at -tpres 0 + phreatic level
  above the mesh): p_u = -49999.99875 per IP (Professional -49997.5,
  0.005% from the discrete load distribution; target tolerance 1e2),
  effective sigyy -49999.99875 vs Professional -50002.5 at the mid point
  (same 0.005%, in tolerance), to_pres profile 10*z EXACT (0 at the top,
  -10 at z=-1) after the -topres fix below.
- ground17 (2D, prescribed volumetric strain -1e-5, E=1e7, nu=0, C=1e-7):
  young_apparent = 1e7 and poisson_apparent = 0 EXACT at every node
  (= Professional digit for digit); the undrained pressure p_u = -100 does
  not enter the apparent moduli (they use the EFFECTIVE incremental
  stresses).

## -topres conversion (bounda.cc)

The head-to-dof conversion of `bounda_dof ... -topres` depends on the source
of the static pressure (`groundflow_phreatic_coord` gained the out-param
`level_source`; groundfl.cc):
- static-pressure-height region (level_source 2): `pres = load - static`
  (unchanged; ground13 verified digit by digit).
- phreatic level (level_source 1): MEASURED on the Professional binary the
  head keeps the load (`h = load`) and `p_total = h - rho*g*z`, so
  `pres = p_dynamic = load - rho*g*level` (uniform under a constant level;
  probes: undrained1 +10, groundflow load -100 -> -90/-105 inside/above).
  The static subtraction of the previous batch was only verified for the
  static-height branch (ground13); for a phreatic level ABOVE the mesh it
  produced p_fixed = 0 instead of the measured -rho*g*z profile.
- no level: `pres = load + rho*g*z` (unchanged).

## Stress recovery / effective output

The stress DOFs keep the EFFECTIVE constitutive value: in materi.cc the
undrained pressure joins `total_new_sig` (the internal-force total stress,
same place as the groundflow pressure) while `new_sig` (constitutive, -> the
sig dofs and element stresses) is untouched. The record is written
read-modify-write per integration point: p_old basis from VERSION_NORMAL
(step start), the buffer preserving the other IP slots from VERSION_NEW
(without this the last IP call clobbers the previous slots of the element
pass - measured [0,0,0,val] on undrained1).

## Records

- `node_dof_previous_step` (database.cc): INTERNAL (external 0) per-node
  snapshot of the node dofs at the start of the step, captured by top() at
  every solving step (before the equilibrium loop, after the control data of
  the step). Consumed by the apparent-moduli operators.
- `element_intpnt_materi_undrained_pressure`: pre-allocated with zeros in
  top() (input meshes) and in step_start (macro-created meshes appear after
  the top() initialization - the parallel element loop cannot allocate).

## Apparent moduli operators (calcul.cc)

`post_calcul -materi_stress -young_apparent/-poisson_apparent` (6.903):
from the INCREMENTAL stress and strain of the last time step
(current node dofs minus node_dof_previous_step):
`K = dp/dvol`, `G = dq/(3*deps_q)` with q and eps_q the deviatoric
measures, `E = 9KG/(3K+G)`, `nu = (3K-2G)/(2(3K+G))`, 0 when the
determination is not possible (zero incremental strain). Needs
`materi_strain_total` initialized (the strain basis). ground17 verified
E = 1e7, nu = 0 EXACT.

## Pending

- `element_intpnt_materi_undrained_pressure_average` registered but not
  filled.
- p_u is stored per integration point; the theory ("element-by-element")
  and the Professional .dbs (identical values at all IPs) are consistent
  with an element-level pressure for uniform states; a non-uniform strain
  field distinguishes per-IP from element-average (untested).
- control_materi_undrained_apply -no mid-run switches (the gate reads the
  control of the current step; switching off mid-run stops the pressure
  increments but keeps the accumulated values).
