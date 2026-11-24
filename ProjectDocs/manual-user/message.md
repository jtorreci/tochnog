# message

## Description

`message switch` (manual Professional 6.856): with `-yes` (default) the
introduction message is printed at the start of the calculation: how to
suppress it, the executable and its build date, and the GPL warranty
notice. `-no` suppresses it. The environment symbol `MESSAGE` (or
`message`) set to `-no` suppresses it for every calculation.

## Usage

```
message switch
```

```
export MESSAGE=-no      # all calculations in this shell
```

## Example

```
message -no
```

Tests: `validation-suite/test-2026/message1` (default, and with
`MESSAGE=-no`) and `message2` (`message -no`), checked by
`scripts/build_safe.sh` on the stdout.

## Differences with the Professional version

- The text names the GNU executable and the GPL warranty notice instead of
  the Professional disclaimer (DIVERGENCES `OUTPUT-TEXT`, class C).
- Before F5 (2026-11-02) the GNU printed no introduction message.
