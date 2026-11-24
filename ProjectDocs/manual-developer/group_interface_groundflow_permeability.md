# group_interface_groundflow_permeability

## Files and functions

- `interface.cc` — in `interface_element()` (lines 390-409). When
  `groundflow_pressure` is active and
  `GROUP_INTERFACE_GROUNDFLOW_PERMEABILITY` is present, a through-interface
  flux couples the pressure dofs of the facing node pairs:
  ```c
  if ( has_pe ) {
    long int ns1 = nnol/2;
    for ( inol=0; inol<ns1; inol++ ) {
      long int jn1 = inol*npuknwn + pres_indx/nder;
      long int jn2 = (inol+ns1)*npuknwn + pres_indx/nder;
      double pres1 = new_dof[inol*nuknwn+pres_indx];
      double pres2 = new_dof[(inol+ns1)*nuknwn+pres_indx];
      double q = pe_iface * ( pres1 - pres2 );
      element_rhside[jn1] -= q;
      element_rhside[jn2] += q;
      element_matrix[jn1*nnol*npuknwn+jn1] += pe_iface;
      element_matrix[jn1*nnol*npuknwn+jn2] -= pe_iface;
      element_matrix[jn2*nnol*npuknwn+jn1] -= pe_iface;
      element_matrix[jn2*nnol*npuknwn+jn2] += pe_iface;
    }
  }
  ```
- `database.cc` — keyword registration (alphabetical, between
  `GROUP_INTERFACE_GROUNDFLOW_CAPACITY` and
  `GROUP_INTERFACE_GROUNDFLOW_TOTAL_PRESSURE_TENSION`): type
  `DOUBLE_PRECISION`, `data_length = 1`, `data_class = GROUNDFLOW`,
  `data_required = GROUP_INTERFACE`.
- Enum `GROUP_INTERFACE_GROUNDFLOW_PERMEABILITY` in `tochnog.h`.

## Implementation details

- The flux is `q = pe * (pres_side1 - pres_side2)` per unit length (2D) or area
  (3D); positive flows out of side 1 into side 2.
- The stiffness matrix entries `+-pe` on the pair give the Newton tangent of
  the flux (symmetric, 4-off-diagonal pattern).
- The interface separates two flow domains; without this record the two sides
  are decoupled (unless connected through the volume elements).
- Works for the interface element types (bar2/quad4/prism6/hex8); the facing
  pairs are the first half vs second half of the node list.

## External dependencies

- Core `db()` accessor; globals `groundflow_pressure`, `pres_indx`, `nnol`.

## Hardcoded parameters / pending refactorings

- No gap dependence: the flux is applied even for an opened interface. A
  refinement would scale `pe` with the opening (crack width).

## K04 fix (2026-10-06): tributary measure

`interface.cc` `interface_element()`, groundflow block: the pair
conductance is `pe_i = pe * w_ip[inol] * iface_measure` (the Lobatto/even
weights and the side length/face area computed for the mechanical springs;
`w_ip` is now released at the end of the function). The same weighting in
the `group_interface_condif_conductivity` block (`k_i`). The capacity
block (`group_interface_groundflow_capacity`): `C_i = 0.5 * C * w_ip[inol %
ns1] * iface_measure` per node (pair storage split over its two nodes),
`element_rhside -= C_i * dp/dt`, `element_matrix` diagonal and
`element_lhside += C_i/dt` (implicit, as `general.cc` for the continuum
capacity). Evidence (Professional 02-08-2026): `ground9` 1/3 -> 0.5, the
same model 2 and 3 wide 0.5; thermal analog 1/3 -> 0.5;
`groundflow_interface` transient 1.0/1.75/1.998 at t = 1/3/10 and
interface capacity 0/2 (1.3333/0.8) all equal, and independent of the
interface length (2 x tall model). Without the implicit capacity the
weighted permeability made `groundflow_interface` diverge (the explicit
storage term was stable only with the doubled conductance).
Pending: `group_interface_groundflow_total_pressure_tension` still adds an
explicit `static_pres - pres` residual without a matrix term (unchanged).
