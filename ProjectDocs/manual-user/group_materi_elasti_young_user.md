# group_materi_elasti_young_user

## Description

With `-yes` the user supplied routine `user_young` (`user.cc`) gives the
young modulus of the group (manual Professional 6.662), for example to
program a stiffness degradation in cyclic loading. The routine receives
the dof values, the old history variables and the old stresses, and the
`group_materi_elasti_young` value (0 when not given). The poisson ratio
of the group is used with it.

The distributed `user_young` is a stub that stops with an error, as in
the Professional; program it and rebuild Tochnog.

## Syntax

```
group_materi_elasti_young_user index -yes|-no
```

## Example

```
group_materi_elasti_young 0 1000.
group_materi_elasti_young_user 0 -yes
```

Tests: `validation-suite/test-2026/group_materi_elasti_young_user1`, `_error1`.

## Differences with the Professional version

- None measured: both stop with "Error: routine user_young is not
  programmed." when the stub is not programmed. The argument list of
  the GNU routine is its own (the Professional source is not available).
