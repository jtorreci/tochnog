# control_system_call

## Description

`control_system_call index integer_value` (manual Professional 6.382): in
control index `index` Tochnog runs a command that you write yourself, the
executable file `tochnog_system_call.sh` of the working directory, with
`integer_value` as its first argument (`$1` in the shell script). Use it to
log the progress of a long run, to copy intermediate results at a
consistent point, or to send yourself a mail when the calculation reaches a
certain index.

## Usage

```
control_system_call index integer_value
```

- `index`: control index.
- `integer_value`: integer (0 or larger) passed as first argument.

## Behaviour

- The script runs at the end of each converged timestep of the index
  (before that step's `control_print_*` output) and once more at the end of
  the index. An index without timesteps gives one call.
- The calculation waits until the script ends. Its standard output goes to
  the Tochnog output; its error output is discarded and its exit status is
  ignored.
- A script without a `#!` line is run by `/bin/sh`.
- If `tochnog_system_call.sh` is missing or not executable, Tochnog prints
  one warning and the calculation continues.
- Linux only (the MS Windows `tochnog_system_call.bat` of the manual is not
  supported).

## Example

```
control_timestep 0 0.1 1.
control_system_call 0 7
```

with the executable `tochnog_system_call.sh`:

```
#!/bin/sh
date >> system_call.out
echo "index value $1" >> system_call.out
```

Test: `validation-suite/test-2026/control_system_call1` (run by
`scripts/build_safe.sh` in an isolated directory with a logging script).

## Differences with the Professional version

- The Professional runs `tochnog_system_call.sh` through the search `PATH`:
  without `.` in the `PATH` the script of the working directory is not run,
  and no message is given. GNU runs the file of the working directory, as
  the manual says (DIVERGENCES `SYSCALL-PATH`, class D).
- The GNU prints a warning when the script is missing; the Professional is
  silent.
- Calls, their order and argument, stdout/stderr and exit status handling
  are equal to the Professional 02-08-2026 (with `.` in its `PATH`).
