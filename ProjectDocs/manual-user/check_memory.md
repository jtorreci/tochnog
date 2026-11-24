# check_memory

## Description

`check_memory -yes` reports the peak memory usage of the calculation just
before the run finishes. The value is read from the OS with `getrusage` and
printed to the standard output in GB. It is useful to verify that a mesh and
its solver data fit in RAM, or to see how close a run came to the machine's
limit.

## Usage

Place it in the data part, as a keyword line:

```
check_memory -yes
```

## Parameters

| Parameter | Meaning                                                    |
|-----------|------------------------------------------------------------|
| `-yes`    | Enable the peak-memory report at the end of the run.       |
| `-no`     | Disable the report (default).                              |

## Example

Minimal input that reports the peak memory used by the calculation:

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
check_memory
   -yes
```

Run it and the log will end with a line like
`Peak memory usage: 0.12 GB.`
