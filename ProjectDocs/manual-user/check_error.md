# check_error

## Description

`check_error -no` suppresses the error messages printed by Tochnog (the ones
starting with "Error"). It is useful to run without certain error messages
being shown, for example while probing variations of an input file.

The suppression only affects messages that are printed through `pri()`;
messages written directly with `cout` are not filtered.

## Usage

Place it in the data part, as a keyword line:

```
check_error -no
```

## Parameters

| Parameter | Meaning                                                     |
|-----------|-------------------------------------------------------------|
| `-yes`    | Keep error messages (default).                              |
| `-no`     | Suppress messages that start with "Error".                  |

## Example

Minimal input that suppresses the "Error" messages:

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
check_error
   -no
```

Run it; any message starting with "Error" is no longer written to the output.
