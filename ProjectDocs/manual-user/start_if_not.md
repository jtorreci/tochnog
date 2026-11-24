# start_if_not / end_if_not

Conditional blocks in the data part (Professional syntax).

## Syntax

```
start_define
  use_damping_boundary true
end_define

bounda_dof 10  -all -vely
start_if use_damping_boundary
  support_edge_normal                    10  1.e6 0.5e6
  support_edge_normal_geometry           10  -right
end_if
start_if_not use_damping_boundary
  bounda_dof 20  -right -velx
end_if_not
```

The records between `start_if <define> ... end_if` are read only when the
`start_define` word is `true`; `start_if_not <define> ... end_if_not` only
when it is `false`. Blocks cannot be nested, and comments are not allowed
inside a block (same rules as `start_if`).

## Why

The GNU parser already supported `start_if ... end_if` (sfnet 2014). The
Professional test suite also uses the negated form
(`start_if_not use_damping_boundary ... end_if_not`, e.g. dynamic6/dynamic7
of the corpus), which failed with "I do not know what to do with:
start_if_not".

## Implementation status

Both forms are supported. See the developer note for the parser details.
