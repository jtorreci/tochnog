# check_memory_usage

## Description

`check_memory_usage -yes` records the maximum memory used by the calculation
and stores it in the data record `check_memory_usage_result` (in GB), so the
peak consumption is available to other features or for post-processing. It is
intended for monitoring memory usage on large calculations. Like
`check_memory`, the value is read from the OS with `getrusage`.

## Usage

Place it in the data part, as a keyword line:

```
check_memory_usage -yes
```

## Parameters

| Parameter | Meaning                                                    |
|-----------|------------------------------------------------------------|
| `-yes`    | Enable peak-memory tracking and store the result in `check_memory_usage_result`. |
| `-no`     | Disable the tracking (default).                            |

## Example

Minimal input that records the peak memory used by the calculation:

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
check_memory_usage
   -yes
```

Run it; at the end of the run the record `check_memory_usage_result` holds the
peak memory in GB and the line `Peak memory usage: ... GB.` is printed.
