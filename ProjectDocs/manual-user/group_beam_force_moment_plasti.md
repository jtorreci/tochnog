# group_beam_force_moment_plasti

## Description

`group_beam_force_moment_plasti index fx1 fy1 fz1 mx1 my1 mz1 fx2 fy2 fz2 mx2 my2 mz2`
(manual 6.589, index = element group): ideally plastic limits of the
forces and moments at the two ends of the beam elements of the group, in
the local beam axes (local x along the beam, see `group_beam_direction_z`),
first node first. All values must be zero or positive; a zero limit is a
hinge from the start.

## Effect on the calculation

- Plastic hinges form at the element ends. While an end force or moment
  is below its limit the beam is elastic; when it reaches the limit it
  stays at the limit (ideally plastic, no hardening) and the element
  deforms plastically at that end. Unloading is elastic.
- The limits act per component, without interaction (no N-M or M-M
  interaction):
  - `mz1`, `mz2`: bending moment about local z (plastic rotation in the
    local x-y plane);
  - `my1`, `my2`: bending moment about local y;
  - `fy`, `fz`: shear force; the end shear is (m1 + m2)/L, so the smaller
    of the two node limits bounds both ends;
  - `mx`: torsion; the smaller of the two node limits applies;
  - `fx`: no effect (the beam has no axial stiffness; combine with a truss
    for the axial force).
- Equilibrium is kept: element forces and nodal reactions stay consistent
  (end shear = (m1 + m2)/L, root moment of a cantilever = shear times
  length). The state is the end forces `element_beam_force_moment`.
- A beam with this record always uses the local-axes beam, also in 2D
  (same results as the classic 2D beam while elastic).
- Hinges only form at element ends (nodes). Mesh densely where hinges
  are expected: with element length h and a distributed load q, the
  moment at the nearest node differs from an interior maximum by at most
  q h^2 / 8.
- Solution procedure: in the time step where a hinge forms, the element
  keeps its elastic stiffness in the system matrix (the forces are
  already plastic), so that step converges linearly; use enough
  `control_timestep_iterations` or small steps near the hinge formation.
  From the next step on the matrix has the hinge. Beyond the collapse
  load a force-controlled calculation has no solution (singular matrix);
  prescribe a displacement to follow the collapse.
- Error: "Error: GROUP_BEAM_FORCE_MOMENT_PLASTI limits should not be
  negative."

Differences with the Professional 02-08-2026: DIVERGENCES
BEAM-PLASTIC-HINGE (the Professional only clips the end force components
and keeps the elastic stiffness, the element is out of equilibrium).

## Example

Propped cantilever (length 4, Mp 1, point load at the middle):

```
group_beam_force_moment_plasti 0  1e10 1e10 1e10 1e10 1e10 1  1e10 1e10 1e10 1e10 1e10 1
```

The hinge at the fixed end forms at P = 16 Mp/(3L) = 1.3333, the collapse
load is 6 Mp/L = 1.5.

Self-checking inputs in `validation-suite/test-2026`:
`group_beam_force_moment_plasti1` (first hinge, redistribution),
`group_beam_force_moment_plasti2` (elastic unloading, residual state),
`group_beam_force_moment_plasti3` (collapse under prescribed displacement),
`group_beam_force_moment_plasti4` (hinge at a node of two elements),
`group_beam_force_moment_plasti5` (shear and torsion limits),
`group_beam_force_moment_plasti_2d1` (2D),
`group_beam_force_moment_plasti_error1` (negative limit).
