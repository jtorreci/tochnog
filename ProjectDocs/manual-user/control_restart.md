# control_restart

## Description

`control_restart index -yes` restarts the calculation at the start of
control index `index` from the initial state (manual Professional
6.368):

- the mesh is the undeformed mesh (`node` reset to the start
  coordinates of the current, possibly refined, mesh);
- the nodal dofs (`node_dof`) are reset to their initial values;
- `time_current` is reset to the initial time;
- the integration-point state (stresses, strains, plastic strains,
  history variables kept per integration point in `element_dof`) is
  rebuilt from the reset nodal dofs (since F4 T6c, 2026-10-11).

This lets you calculate path dependent behaviour (plasticity, damage)
completely from the start, typically on a mesh refined by an earlier
`control_refine_locally` pass.

Before F4 T6c the GNU kept the old integration-point state after a
restart: a run restarted after loading continued from the loaded
stresses (an elastic probe gave sigyy +52.08 instead of -52.08).

## Usage

```
control_restart  index  -yes
```

| parameter | meaning |
|---|---|
| `index` | control index at whose start the restart is done |
| `-yes` / `-no` | restart or not |

## Example

Unit quad4, top edge pushed down with vely -1e-2, Drucker-Prager
(`validation-suite/test-2026/control_restart_reset_ip1.dat`):

```
control_timestep 10 0.1 1.0
control_restart 20 -yes
control_timestep 30 0.1 0.5
```

The result at t 0.5 equals a fresh 0.5 run (Professional 02-08-2026):
sigyy -24.1835 at node 4. The elastic copy (`control_restart_reset_ip2`)
gives sigyy -52.0833, eptyy -5e-3.

## Differences with the Professional version

- The Professional sets `number_timesteps` to 0 on a restart
  (`number_timesteps_total` is kept). The GNU has no `number_timesteps`
  record, so there is nothing to reset (output-only difference).
- Element records outside the integration-point dofs: `element_truss_force`
  is kept by both binaries (probe: 1.5 in both), the undrained pressure
  equals the Professional (probe). Interface, spring, contact-spring,
  tendon and damage effective-stress records are not verified yet
  (DIVERGENCES C58, BACKLOG C58).
