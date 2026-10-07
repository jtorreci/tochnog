# groundflow_velocity_sign

## Description

Sign convention of the groundflow (Darcy) velocity when the
`groundflow_velocity` unknown is initialized.

The Professional manual (section 8.2, "Groundflow velocities") defines

```
v_i^g = k_ij * dh/dx_j
```

that is, the groundflow velocity follows from the permeability times the
**positive** gradient of the hydraulic pressure head `h` (the
`groundflow_pressure` dof, label `-pres`).

This matters when the groundflow equation is coupled with heat transport
(`-groundflow -condif`): the convection velocity of the temperature
equation is exactly the groundflow velocity (manual 8.3, "Heat transport
in ground water flow"). If the sign of the velocity is inverted, the
convection direction is inverted and an upwind node gets the downstream
temperature instead of the upstream one.

GNU (from the sfnet 2014 sources) had this term negated by a 2013 change
(`-pe*grad(h)`). The sign is now the Professional one (`+pe*grad(h)`), so
a coupled `-groundflow -condif` model produces the same upwind direction
as the Professional binary.

## Usage

The velocity is a solved dof; it is switched on by listing
`groundflow_velocity` in the initialisation block. No record has to be
set for the sign: it is the built-in convention.

```
groundflow_pressure
groundflow_velocity
condif_temperature
end_initia

group_type 0 -groundflow -condif
group_groundflow_permeability 0 10.
group_condif_density 0 1.0
group_condif_capacity 0 1.0
group_condif_conductivity 0 0.1
```

## Parameters

None. The velocity value is determined by `group_groundflow_permeability`
and the head gradient.

## Example

`condif5` (sfnet corpus): a 1D column with nodes at x=0,1,2. The head is
prescribed as `-2, -1, 0` at the three nodes, so `dh/dx = +1`. With
`k = 10` the groundflow velocity is `+10` and the temperature
prescribed `1` at x=0 is convected to the middle node, which reaches the
target `temp = 1.0`. Before the fix the velocity was `-10`, the upwind
node was x=2 (`temp = 0`), and the middle node settled at `1.1e-16`.

```
bounda_dof 0  1 -temp
bounda_time    0  0. 1. 100. 1.
bounda_dof 1  3 -temp -pres
bounda_time    1  0. 0. 100. 0.
bounda_dof 2  1 -pres
bounda_time    2  0. -2. 100. -2.
bounda_dof 3  2 -pres
bounda_time    3  0. -1. 100. -1.

group_groundflow_permeability 0  10.
group_condif_conductivity     0  0.1
```

Result: node dofs `temp/gvelx/pres` = `1,10,-2` / `1,10,-1` / `0,10,0`,
identical to the Professional binary.

## Differences with the Professional version

- **groundflow velocity output (`ground1`)** (deliberate divergence decided by
  the project owner on 2026-10-07; AUDIT K10; SEGUIMIENTO Descarte). GNU: v = k dh/dy =
  -0.05 (manual), at the nodes and at the integration points. Professional
  02-08-2026: 0 at the nodes and at the integration points (post-processing
  only; the pressure field agrees). Status: deliberate divergence (project
  decision 2026-10-07): the GNU follows the manual; the Professional
  behaviour is also Pro-suspect. Root cause measured 2026-10-07: the
  Professional zeroes the groundflow velocity where the calculated total
  pressure reaches `groundflow_total_pressure_limit`, whose default is 0
  (manual 6.587); with `groundflow_total_pressure_limit 1.e10` the
  Professional gives -0.05 = GNU (see `groundflow_total_pressure_limit.md`,
  DIVERGENCES TPL, pending owner decision). Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **`gvelx` sign in the interface-coupled groundflow path
  (`groundflow_interface`)** (known difference, not decided; BACKLOG C2).
  GNU: gvelx -1.99962e-6. Professional 02-08-2026: +1.99960e-6 (pressure
  field identical). Status: Open, LOW; root cause not established. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **groundflow `pres` dof datum** (known difference, not decided; AUDIT
  K09). GNU: Stores the pressure (`ground4` pres 0). Professional
  02-08-2026: Stores the head h (500), as its manual says. Status:
  Convention difference (audit BOTH-plausible); total pore pressure and
  stresses identical. No decision. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).

- **Decision 2026-10-07 (project owner):** align the GNU with manual 6.587 and the Professional on `groundflow_total_pressure_limit` (default 0, cap the calculated total pressure, zero the groundflow velocity where the limit is reached). This reverses the earlier K10 decision. Implementation scheduled for F4; see `ProjectDocs/DIVERGENCES-PROFESSIONAL.md` (TPL, K10).
