# group_materi_history_variable_user / _parameters

## Description

With `-yes` the user supplied routine `user_materi_history_variable`
(`user.cc`) updates the material history variables of the group
(`materi_history_variable` in the initialization part) after the stress
update of every integration point (manual Professional 6.671, 6.672).
`group_materi_history_variable_user_parameters` gives the routine any
number of parameters.

The distributed routine is a stub that stops with an error, as in the
Professional; program it and rebuild Tochnog.

## Syntax

```
group_materi_history_variable_user index -yes|-no
group_materi_history_variable_user_parameters index p0 p1 ...
```

## Example

```
materi_history_variable 2
...
group_materi_history_variable_user 0 -yes
group_materi_history_variable_user_parameters 0 1. 2. 3.
```

Tests: `validation-suite/test-2026/group_materi_history_variable_user_parameters1`,
`group_materi_history_variable_user_error1`.

## Differences with the Professional version

- None measured: both stop with "Error: routine
  user_materi_history_variable is not programmed."; the parameters
  alone run unchanged. The GNU argument list is its own.
