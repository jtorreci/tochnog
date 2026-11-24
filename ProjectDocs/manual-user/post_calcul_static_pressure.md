# post_calcul_static_pressure

## Description

`post_calcul_static_pressure static_pressure` (manual Professional 6.920).
The static pressure of the complete domain for postprocessing. It replaces
the static pressure of the phreatic level (or of
`post_calcul_static_pressure_height`) in:

- `post_calcul -groundflow_pressure -static_pressure` (the value itself),
- `-dynamic_pressure` (total minus this static pressure),
- the p_dynamic of `post_calcul -materi_stress -safety_piping`.

`-total_pressure`, `-safety_lifting` and the calculation itself are not
changed. `node_static_pressure` / `node_dynamic_pressure` still override
per node.

## Usage

```
post_calcul -groundflow_pressure -static_pressure -groundflow_pressure -dynamic_pressure
post_calcul_static_pressure -7.
```

## Example

`validation-suite/test-2026/post_calcul_static_pressure1.dat` (Professional
`ground14`): static -7 everywhere, dynamic = total + 7, piping at the
bottom 4.333333.

## Differences with the Professional version

None measured (same node values).
