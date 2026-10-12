# materi_plasti_rho

## Description

`materi_plasti_rho` (initialisation part, Professional manual 4.27) adds
the plastic kinematic hardening tensor (back stress) rho to the
`node_dof` records, with components `rhoxx rhoxy rhoxz rhoyy rhoyz rhozz`.
It is the state variable of `group_materi_plasti_kinematic_hardening`:
the yield function is evaluated on `sigma - rho`, and rho grows with the
plastic strain increment, `inc_rho = H * inc_epp`.

rho is an integration-point history state: it is kept per integration
point between steps (since C68, 2026-10-12), as the hardening parameter
`kap`. The node values are a recovery (averaging) of the
integration-point values with the same weights as the strains; the
constitutive update always uses the integration-point value. Before
2026-10-12 GNU fed the averaged node value back into the update (11 %
off the Professional in the reactions of a four-element test).

## Usage

```
materi_plasti_rho          (initialisation part)
group_materi_plasti_kinematic_hardening <group> <H>
```

Also `materi_stress` and `materi_strain_plasti` should be initialised.

## Parameters

- `H` (`group_materi_plasti_kinematic_hardening`): kinematic hardening
  modulus, `inc_rho = H * inc_epp`.

## Example

```
number_of_space_dimensions 2
materi_velocity
materi_displacement
materi_stress
materi_strain_total
materi_strain_plasti
materi_plasti_rho
end_initia
...
group_materi_plasti_vonmises 0 24.3
group_materi_plasti_kinematic_hardening 0 2000.
```

Complete input: `validation-suite/test-2026/rho_ip_state1.dat`.

## Differences with the Professional version

None known for the integration-point state (rho reactions within 7e-6
relative of the Professional 02-08-2026 on `rho_ip_state1`).
