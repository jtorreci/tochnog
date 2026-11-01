# ranges

## Description

A range `-ra ... -ra` gives several integer numbers at once. It can be
used for the index of a record (`element_group -ra ... -ra 1` puts every
listed element in group 1) and inside the values of a record that takes
node or element numbers (`bounda_dof 0 -ra ... -ra -velx`). The
syntax follows the Professional user manual, section "Ranges".

Between the two `-ra` keywords you can write, in any order and number:

| form | meaning |
|---|---|
| `n` | the single number `n` |
| `-from a -to b` | every number from `a` to `b` (`a` <= `b`; `a` = `b` is the single number `a`) |
| `-from a -to b -step s` | `a`, `a+s`, `a+2s`, ... up to `b` (`s` > 0; `b` is included only when it is reached) |

Important: **`-ra a b -ra` is a LIST of the two numbers `a` and `b`,
never the range `a..b`.** `-ra 1 4 -ra` is nodes 1 and 4 only; the
range 1, 2, 3, 4 is `-ra -from 1 -to 4 -ra`. Recommendation: always
spell a range of consecutive numbers with `-from`/`-to`, and keep plain
lists for scattered numbers. (Some Professional manual pages call a
list such as `-ra 5 4 8 -ra` "a number range"; it is still a list.)

Every word of an index range is read like any other input word: a
`start_define` word can stand for a number, and comments and include
files work inside the range. Arithmetic words are not numbers there (as
in the Professional).

Errors (the run stops):

- `-from` not followed by a number, `-to` and a number.
- a descending block `-from 7 -to 5`: `Error in range, a -from is
  larger than a -to`.
- `-step 0` or a negative step.
- a number with a decimal point (`9.0`) inside the range.
- an empty range `-ra -ra`.
- `-all` instead of a range for a record index (`-all` is only allowed
  for values).

## Usage

```
<record> -ra <items> -ra <values>          index range
<record> <index> ... -ra <items> -ra ...   value range
<items>: n | -from a -to b [-step s]  (repeated, any order)
```

## Example

```
element_group -ra -from 1 -to 8 -step 2 -ra 1      ( elements 1 3 5 7 )
element_group -ra -from 1 -to 2 -from 7 -to 8 -ra 1 ( elements 1 2 7 8 )
bounda_dof 0 -ra 1 -from 5 -to 5 9 -ra -velx      ( nodes 1 5 9 )
```

Tests: `validation-suite/test-2026/range_*.dat`.

## Differences with the Professional version

- Ranges expand to the same numbers and fail on the same inputs as in
  the Professional 02-08-2026 (since F4 T6f, 2026-10-12). The wording of
  some error messages differs (for example `Error in range reading: 9.0`
  for an index, `Error in range, a -step should be larger than 0`), and
  the GNU stops with exit code 1 where the Professional prints the error
  and exits with 0.
- Outside a range the GNU still reads an integer value written with a
  decimal point (`9.0` read as 9); the Professional rejects it
  (DIVERGENCES `INT-DOT`).

Before F4 T6f the GNU expanded an index `-step` range to wrong indices,
rejected `-from a -to a` and a second `-from` block in an index range,
read index range numbers without defines or include files, silently
returned nothing for a descending value block and read `9.0` as 9 in a
value range.
