# check_warning

## Description

`check_warning -no` suppresses the warning messages printed by Tochnog (the
ones starting with "Warning"). It is useful to run without certain warning
messages being shown, for example while probing variations of an input file.

The suppression only affects messages that are printed through `pri()`;
messages written directly with `cout` are not filtered. Features that emit
warnings through `pri()` automatically honor this switch.

## Usage

Place it in the data part, as a keyword line:

```
check_warning -no
```

## Parameters

| Parameter | Meaning                                                       |
|-----------|---------------------------------------------------------------|
| `-yes`    | Keep warning messages (default).                              |
| `-no`     | Suppress messages that start with "Warning".                  |

## Example

Minimal input that suppresses the "Warning" messages:

```
control_geometry
   cartesian
control_time
   end 0.1
   dt 0.01
control_print
   history 0
   step 10
materi_elasti_young
   0 2.e7
materi_elasti_poisson
   0 0.3
materi_density
   0 2500.
check_warning
   -no
```

Run it; any message starting with "Warning" is no longer written to the output.
