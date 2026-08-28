# post_calcul_materi_stress_force

## Description

`post_calcul -materi_stress -force` (manual Professional 6.913)
calculates the normal force, shear force and moment(s) of
isoparametric elements with a SINGLE element over the structure
thickness (sheet piles, tunnel shells, ...). Since LOT 5 (2026-08-28)
the section forces are computed from the ELEMENT INTERNAL FORCES
f_elem = ∫Bᵀσ dV (the equilibrium-based method, see the notes): the
section resultant over an end face = the sum of the internal forces of
the face nodes = the free-body statics of the applied loads (EXACT for
the converged solution, like the Professional's node_dof_calcul). The
results are written per node to `node_dof_calcul`.

The 2D result is a set of 9 items per node:

| item | meaning |
|------|---------|
| `norx_sig` `nory_sig` | normal force per unit length, GLOBAL PLOT vector components (drawn in the structure thickness direction) |
| `nors_sig` | normal force per unit length, the PHYSICAL size (design value) |
| `shex_sig` `shey_sig` | shear force per unit length, plot components |
| `shes_sig` | shear force per unit length, physical size (always positive - only the size is available, manual 6.913) |
| `momx_sig` `momy_sig` | moment per unit length, plot components |
| `moms_sig` | moment per unit length, physical size |

The 3D result has 16 items (nor, she, mom1, mom2 with 3 components +
size each):

| item | meaning |
|------|---------|
| `norx_sig` `nory_sig` `norz_sig` | normal force per unit length, plot components |
| `nors_sig` | normal force per unit length, physical size |
| `shex_sig` `shey_sig` `shez_sig` | shear force per unit length, plot components |
| `shes_sig` | shear force per unit length, physical size |
| `mom1x_sig` `mom1y_sig` `mom1z_sig` | moment in THICKNESS direction per unit length, plot components |
| `mom1s_sig` | that moment, physical size |
| `mom2x_sig` `mom2y_sig` `mom2z_sig` | moment in LENGTH direction per unit length, plot components |
| `mom2s_sig` | that moment, physical size |

The 3D numerical integration (hex8/hex27) is implemented (lot 3):
`mom1` = the thickness bending moment (radial bending moment in a
tunnel shell, thickness bending moment in a sheet pile) and `mom2` =
the length-direction bending moment (manual 6.913), both from the
σ_nn moment contributions with a distance relative to the middle of
the element.

Definitions (manual 6.913): `nor` = normal stresses σ_nn integrated
over the thickness (a positive value = tension); `she` = shear stresses
σ_nt integrated over the thickness (Tochnog always outputs a positive
value, only the size); `mom` (2D) = the moment contributions of σ_nn
with a distance in the thickness direction dt relative to the middle of
the element, integrated over the thickness. All results are per unit
length l of the element (2D plane: l = 1; axisymmetric: l = 2π·radius).

## Input syntax

```
post_calcul -materi_stress -force
post_calcul_materi_stress_force_element_group group_0 group_1 ...
post_calcul_materi_stress_force_reference_point x_0 y_0 x_1 y_1 ...
[post_calcul_materi_stress_force_average -yes | -no]
[post_calcul_materi_stress_force_outer -yes | -no]
[post_calcul_materi_stress_force_plot_switch sw_nor sw_she sw_mom]
[post_calcul_materi_stress_force_thickness_switch sw_0 sw_1 ...]
[post_calcul_materi_stress_force_direction_exclude dir_x dir_y dir_z]   (3D)
[post_calcul_materi_stress_force_direction_include dir_x dir_y dir_z]   (3D)
```

Key records (see their individual pages for details):

- `post_calcul_materi_stress_force_element_group` (mandatory): the
  target element groups.
- `post_calcul_materi_stress_force_reference_point` (mandatory in 3D;
  in 2D a warning + the default (0,0) when absent): one point per
  element group. It defines the thickness direction of the structure:
  the plot vectors are oriented outward/inward consistently from this
  point (manual 6.914: for a sheet pile use a point at a large
  perpendicular distance).
- `post_calcul_materi_stress_force_average` (default -yes, quad9 only
  in 2D): the middle-plane nodes receive the average of the two
  opposing end faces.
- `post_calcul_materi_stress_force_outer` (default -no): only the
  nodes at the maximum distance from the reference point receive
  values (nicer vector plots).
- `post_calcul_materi_stress_force_plot_switch` (3 switches in 2D):
  -yes inverts the drawing direction of each vector item.
- `post_calcul_materi_stress_force_thickness_switch`: 3D concept
  (shortest/longest element direction), validated but not consumed in
  2D.
- `post_calcul_materi_stress_force_direction_exclude` /
  `_include`: 3D concepts; in 2D they produce a warning and are
  ignored.

## Example

Cantilever 2D beam (8 quad9 elements, L = 8, h = 1) with a vertical
tip load P = 1e-2 at node 51, clamped at x = 0 (nodes 1, 18, 35):

```
post_calcul_materi_stress_force_element_group 0
post_calcul_materi_stress_force_reference_point 4.0 1000.0
post_calcul -materi_stress -force
...
control_print_materi_stress_force 410 -all
```

The file `materi_stress_force.410` contains the section forces per
node: at the section at distance (8−x) from the tip, `moms_sig` =
P·(8−x) (within 1-3% with quad9), `shes_sig` ≈ P and `nors_sig` ≈ 0
(no axial load). The middle-plane nodes of each element receive the
average of the two adjacent sections.

## Notes

- 2D supports `-quad4` and `-quad9`; 3D supports `-hex8` / `-hex27`
  (error for other element types in the target groups, manual 6.913).
- 3D face selection (manual 6.909/6.911): specify
  `post_calcul_materi_stress_force_direction_exclude` (typically the
  tunnel length axis) or `_include` (typically the sheet pile height
  direction) so that Tochnog knows on which element sides the forces
  act; both together are an error. The direction must leave exactly 4
  element sides consistent with it (manual 6.913), otherwise the
  element is skipped with a warning. The two END faces where the
  forces are primarily calculated (manual 6.908) are the 2 sides most
  perpendicular to the thickness direction defined by the reference
  point.
- 3D thickness direction in a face: the SHORTEST element direction by
  default; `post_calcul_materi_stress_force_thickness_switch -yes`
  switches to the LONGEST (manual 6.917).
- `average` is available for quad9 (2D) and hex27 (3D) elements.
- **LOT 5 (2026-08-28): the section forces are the EQUILIBRIUM statics
  of the loads.** The section resultant over an end face = the sum of
  the ELEMENT INTERNAL FORCES (f_elem = ∫Bᵀσ dV, the consistent nodal
  forces of the integration-point stresses) of the face nodes. For the
  converged solution the internal forces are in equilibrium with the
  applied loads by construction, so the section forces are the EXACT
  free-body statics: `nors` = N, `shes` = V and `moms` = M of the
  loads (verified in the arness: gforce7q4 N/V EXACT, M 99.84%; the
  fixed-fixed beam identity |M_end|+|M_center| = pL²/8 EXACT). The
  previous stress-field integration (lots 2-4) carried the pollution
  of the mixed u-σ scheme (the section values depended on the raw σ
  field, which is not in equilibrium in the GNU's coarse runs). The
  sign convention: `nor` = n̂·R (positive = tension), `she` = |t̂·R|
  (always positive), `mom` = Σ (n̂·f)·arm about the middle of the face.
- The moments of the section include the full weak-form moment of the
  internal forces about the section point: for a PURE-SHEAR state the
  `moms` = the moment of the reaction couple of the shear block
  (e.g. 0.1923 = 0.3846·0.5 in msf_shear/qsri_patch_s), while the
  manual's σ_nn·dt integral is 0 - the internal-force definition is
  the consistent one for the equilibrium method (documented in the
  tests).
- **Equilibrium requirement**: the section forces are the free-body
  statics ONLY when the internal forces are in equilibrium with the
  loads (the converged quasi-static solution). For the GNU's coarse
  multi-step quad9/hex8 runs the recovered σ state is NOT in
  equilibrium (a pre-existing property of the staggered scheme,
  documented in the L4 baseline as the N/V pollution) - in those cases
  the equilibrium section values deviate from the targets (the arness
  gforce7 shows 5×; the Professional's exact values require the
  equilibrium σ state, which its solver produces). The section forces
  of a NON-equilibrium σ are not meaningful - use the fine meshes /
  converged runs for design values.
- 3D force-loaded validation (LOT 4-5): the cantilever hex27 with a
  real tip load (msf_cant3d_hex27): mom1 = P·(L−z) within 1% (LOT 5)
  and the section shear = P within 4% (the polluted 0.08·P of the
  integration-point field is gone - the equilibrium shear is the free
  body of the loads). The tunnel ring (msf_tunnel3d): the section
  forces now match the Tochnog Professional's node_dof_calcul values
  digit-for-digit (nor 0.0998068, shes 0.013727, mom1 2.0e-4 - the
  FE-discretized values, NOT the analytic 0.1 of the prescribed
  field).
- The results are per unit length l; the x/y (and z) components are
  ONLY convenient global plot vectors - the physical design value is
  the size (`nors_sig`, `shes_sig`, `moms_sig`).
- The stresses are read from the ELEMENT integration points
  (`element_dof`, the constitutive stresses the element used in the
  step - "the element forces needed for this option are setup in a
  timestep", manual 6.913), so `materi_stress` must be present in the
  initia section AND `options_element_dof -yes` (the default; an error
  is raised otherwise). For the solved 3D models with the
  `derivatives` keyword, where the staggered element loop does not
  propagate the strain into the element integration points (measured
  in gforce10/gforce13), the recovered nodal stresses are used instead
  (documented fallback with a warning).
- At least 1 timestep must be done (the stresses come from the
  solution of a step; manual 6.913: "At least 1 timestep should be
  done").
- Accuracy note (2D, verified analytically): with `-quad9` the section
  forces of a converged run are the EXACT free-body statics (msf_beam2d:
  moms = P·(8−x) and shes = P to 6 digits at every section; the shear
  pollution band of the stress-field integration is gone). The
  quad4 1-in-thickness shear locking shows in the DEFLECTION and the σ
  field, NOT in the section forces (the internal forces of the
  equilibrated locked solution still satisfy the statics of the loads
  - the qsri family verifies mom = P·L EXACT with and without the SRI).
  A dedicated test (msf_shear) verifies the exact shear value on a
  uniform simple-shear field.
- Axisymmetric 2D (LOT 4-5, verified): l = 2π·radius (the radial
  coordinate of the element centroid, manual 6.911); the element
  internal forces carry the physical circumference 2π·r in the IP
  volumes, so the per-unit-length values are the section forces per
  unit CIRCUMFERENCE. Verified with msf_axisym: nor = σ_zz·t EXACT,
  she = 0, mom = 0 (the exact free-body statics of the ring).
- `outer -yes` and `plot_switch` are implemented in 2D (3 switches)
  and 3D (4 switches, manual 6.916).
