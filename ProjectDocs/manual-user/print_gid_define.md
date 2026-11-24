# print_gid_define

Manual Professional 6.983.

## Description

With `-yes` (default) the Professional uses the names of `start_define ... end_define` instead of
the geometry names in its GiD files. The GNU GiD files hold no geometry names, so the record is
accepted (`-yes` or `-no`) and has no effect.

## Syntax

```
print_gid_define  switch
```

## Differences with the Professional version

None visible: no effect was measured in the Professional 02-08-2026 either.

## Tests

`validation-suite/test-2026/print_gid_define1` (accepted, mesh unchanged).
