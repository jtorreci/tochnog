# control_post_element_force

## Description

`control_post_element_force index switch` (manual Professional 6.260):
with `-no`, the `post_element_force` records are not evaluated in the
timesteps of that index (saves CPU time); the other post records still
are. Default `-yes`. A `post_apply` / `control_post_apply -no` also skips
them.

## Usage

```
control_post_element_force index switch
```

## Example

```
post_element_force 0 0. 1. 0. 0. 0.5 0.
post_element_force_geometry 0 -geometry_line 1
control_post_element_force 1 -no
```

Self-checking input `validation-suite/test-2026/control_post_element_force1`
(force frozen at 3.2967 in index 1, final 12.0879, = Professional
02-08-2026). As for `control_post_apply`, the skipped result is evaluated
once more at the end of the run.

## Differences with the Professional version

None measured.
