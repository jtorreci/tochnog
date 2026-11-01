# post_apply

## Description

`post_apply switch` (manual Professional 6.900, global) and
`control_post_apply index switch` (6.259, per control index): with `-no`
the post-processing commands (records with "post" in the name:
`post_point`, `post_node`, `post_line`, `post_quadrilateral`,
`post_calcul`, `post_integrate`, `post_element_force`, the `post_*_summed`
records, ...) are not evaluated. Only the `post_node_rhside_ratio` is
evaluated always. Default `-yes`. Useful to save time in large runs, or in
the control indices whose post results are not needed.

## Usage

```
post_apply switch
control_post_apply index switch
```

| Parameter | Meaning |
|-----------|---------|
| `switch` | `-yes` (default) evaluates the post records, `-no` skips them |
| `index` | (control form) the control index the switch applies to |

## Behaviour (measured on the Professional 02-08-2026)

- `post_apply -no`: the post result records are never created (also not
  in the final `.dbs`, no history file).
- `control_post_apply index -no`: in the timesteps of that index the post
  results keep the value of the last evaluated step (a record that did not
  exist yet stays absent).
- The control record of the current index overrules the global one
  (`post_apply -no` with `control_post_apply 1 -yes`: evaluated in index 1
  only).
- Unless `post_apply -no` is set, the post records skipped in the last
  timesteps are evaluated once more at the end of the run, so the final
  `.dbs` (and the `target_*` checks) hold the end state.

## Example

```
post_point 0 0.5 0.5
control_timestep 0 0.1 0.5
control_timestep 1 0.1 0.5
control_post_apply 1 -no
```

Self-checking inputs in `validation-suite/test-2026/`: `post_apply1`
(global -no: no post record), `post_apply2` (global -no, control 1 -yes),
`control_post_apply1` (index 1 frozen at 0.0015, final 0.0055). The frozen
history is checked by `scripts/build_safe.sh` on the `control_print_history`
file. All equal the Professional 02-08-2026.

See also `control_post_element_force` and `control_post`.

## Differences with the Professional version

None measured. The GNU has no `post_bounda_force_summed` /
`post_solver_diagonal_*` records (the Professional keeps them with
`post_apply -no`).
