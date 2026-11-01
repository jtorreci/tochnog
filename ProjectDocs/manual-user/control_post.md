# control_post

## Description

`control_post index switch` (manual Professional 6.258): with `-yes` the
post-processing commands are evaluated for control index `index`, even if
it performs no timesteps. Default: not set.

## Usage

```
control_post index switch
```

## Behaviour

The GNU, as the Professional 02-08-2026, evaluates the post records in
every control index, also in an index without timesteps; `control_post
-yes` therefore does not change the output, and the Professional gives the
same result with `-no` or without the record. Self-checking input
`validation-suite/test-2026/control_post1` (`control_post_apply 1 -no`,
then `control_post 2 -yes`: index 2 evaluates `post_element_force`
12.0879).

## Differences with the Professional version

None measured.
