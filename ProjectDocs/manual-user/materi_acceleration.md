# materi_acceleration

Acceleration dofs added to the node_dof records (manual Professional 4.11).

## When to use

Declare `materi_acceleration` in the initialization part (together with
`materi_velocity` and usually `materi_displacement`) when the calculation
needs accelerations as nodal results, or when an acceleration is prescribed
as a boundary condition (base excitation / earthquake input):

```
echo -no
number_of_space_dimensions 2
materi_displacement
materi_velocity
materi_acceleration
materi_stress
end_initia
```

## Meaning

The GNU dynamics scheme solves **velocities** as the primary unknowns. The
acceleration dofs are *derived records*: at the end of every time step

```
acc = ( v_new - v_old ) / dt
```

is written to the `accx`, `accy`, `accz` (3D) slots of the node_dof record.
Because the GNU updates displacements as `dis = dis_old + v*dt`, declaring
`materi_acceleration` does not add equations to the system: it only adds the
record and its labels.

## Prescribing an acceleration (base excitation)

A `bounda_dof` on `-accx`/`-accy`/`-accz` does **not** constrain the derived
record directly. It imposes the *velocity* bound

```
v_new = v_old + a*dt
```

on the same direction, with `a` taken from the `bounda_time` diagram of the
record. The acceleration then integrates into the velocity and the
displacement exactly like a base motion:

```
node 1  0
node 2  1
element 1  -truss 1 2
...
bounda_dof   10  1 -accx
bounda_time  10  1.234
control_timestep 10  1.e-1 10.
```

Here node 1 moves with constant acceleration 1.234; after 10 s its
displacement is `61.7` (`0.5*a*t^2`, Euler integration of the velocity
bounds, 62.317 with dt = 0.1 — matches the Professional binary run to the
digit).

## Output

`-accx` and friends are regular dof labels: they can be used in
`control_print_history`, `post_point_dof`, `target_item`, etc. For free
nodes the printed value is the backward difference of the converged
velocities; for a node with a prescribed acceleration it equals the
prescribed value.

## Notes

- The order of the dof groups in the initia part defines the column order
  of the node_dof record, as usual (the labels resolve independently of
  that order).
- SMC accelerogram files (`bounda_time_smc`) are not implemented yet; see
  `bounda_time_smc.md` / the PENDING notes in the seguimiento.
