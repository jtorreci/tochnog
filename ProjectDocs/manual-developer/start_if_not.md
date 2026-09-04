# start_if_not / end_if_not

Developer notes for the negated conditional data-part blocks.

## Where

`input.cc`, `input_read_string()`:

- the `start_if` handler was generalized to also accept `start_if_not`; the
  module-level flag `using_if_not` (`str[9]=='n'`) remembers the block kind
  so the skip loop and the end handler accept the matching terminator
  (`end_if` / `end_if_not`). `apply_if` is inverted for the negated form.
- the data-part value loop (record reading) stops on all four tokens
  (`start_if`, `start_if_not`, `end_if`, `end_if_not`) so they are never
  swallowed as record values.
- the include-file and data_ignore skip loops were left unchanged (the
  four tokens are not expected inside those constructs).

## Behavior

`start_if_not X`: looks up the `start_define` word `X` (must exist, values
`true`/`false`); the block is applied when `X` is `false`. Skipped content
is consumed token by token up to the matching end marker. Nested blocks and
comments inside a block are rejected as before.

## Verification

dynamic6/dynamic7 of the corpus use `start_if use_damping_boundary` and
`start_if_not use_damping_boundary` on the same define; both parse now.
dynamic6 reaches the time loop and diverges in the velocity solve (physics
blocker of the wave family); dynamic7 is blocked one record later by
`post_calcul -materi_stress -young_apparent` (pending operator family).
