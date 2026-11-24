# truss_rope_apply

## Description

`truss_rope_apply switch` (manual Professional 6.1092, global) and
`control_truss_rope_apply index switch` (6.393, per timestep index): with
`-no`, any truss rope data (`group_truss_rope -yes`) is ignored, so the
trusses also carry compression. Convenient to test an input "just linear".

## Usage

```
truss_rope_apply switch
control_truss_rope_apply index switch
```

| Parameter | Meaning |
|-----------|---------|
| `switch` | `-yes` (default) applies the rope data, `-no` ignores it |
| `index` | (control form) the `control_timestep` index the switch applies to |

The control record of the current timestep overrules the global one.

## Example

```
group_truss_rope 0 -yes
truss_rope_apply -no
```

Self-checking inputs in `validation-suite/test-2026/`: `truss_rope_apply1`
(global -no: compressed truss -5.5), `control_truss_rope_apply1` (control
10 -no: -5.5), `truss_rope_apply2` (global -no, control 10 -yes: rope, 0).

## Differences with the Professional version

- The Professional 02-08-2026 accepts both records but keeps the rope
  (force 0 in compression) with `-no`; the GNU follows the manual and
  ignores the rope data (DIVERGENCES ROPE-APPLY, class D).
