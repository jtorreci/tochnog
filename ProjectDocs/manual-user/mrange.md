# mrange and mstring

## Description

Two initialization-part limits (manual Professional 4.55 / 4.56):

- `mrange maximum_range_length`: the maximum length of a range
  `-ra ... -ra` in the data part. A longer range stops the program with
  `Error: mrange exceeded.` / `Increase it in the initialisation part.`
  (the Professional message).
- `mstring maximum_number_of_strings`: the maximum number of strings in a
  `start_define ... end_define` block. More strings stop the program with
  `Error: maximal number of strings is N` / `Else, increase mstring in
  the initialisation part.`

They are useful to catch a wrong range or define early (a small limit), or
to state the size an input needs.

## Usage

```
echo -no
number_of_space_dimensions 2
...
mrange 1000
mstring 50
end_initia
```

Both go anywhere in the initialization part after
`number_of_space_dimensions` (`echo` first and `number_of_space_dimensions`
second stay mandatory). The value is an integer larger than 0.

Counting (measured on the Professional 02-08-2026, the same in the GNU):

- a value range (a list in the values of a record, such as the nodes of
  `bounda_dof 0 -ra 1 2 3 4 -ra -velx`) counts its words: 4 here;
  `-ra -from 1 -to 4 -ra` also counts 4;
- an index range (`group_materi_elasti_young -ra 0 1 2 -ra 1000`) counts
  its words plus the two `-ra` markers (5 here), and its expanded list
  counts too (`node -ra -from 5 -to 20 -ra ...` needs 16).

Without `mrange` / `mstring` the GNU limits are the compiled ones: `MRANGE`
500000 (`tochnog.h`) and 799 strings per define (`MSTRING` 800 in
`input.cc`). A larger value is an error that names the constant to raise
and recompile.

Example (`validation-suite/test-2026/mrange1.dat`):

```
mrange 4
end_initia
...
bounda_dof 0 -ra 1 2 3 4 -ra -velx
```

Tests: `mrange1`, `mrange_error1` (list of 4 with `mrange 3`),
`mrange_error2` (`-from 1 -to 4` with `mrange 3`), `mrange_error4` (3-entry
index range with `mrange 4`), `mrange_error3` (above
`MRANGE`), `mstring1`, `mstring_error1`, `mstring_error2` (above
`MSTRING`).

## Differences with the Professional version

- `mstring`: the Professional 02-08-2026 rejects every `mstring` record,
  whatever its value or position ("Error in initialisation part, double
  declared: mstring"); the GNU follows the manual (DIVERGENCES
  `MSTRING-PRO`, class D).
- `mrange` above the compiled `MRANGE` (500000) is an error in the GNU; the
  Professional accepts 100000000 (DIVERGENCES `MRANGE-COMPILED`, class C).
- The expanded length of a value range (`-from a -to b` in a node list)
  is checked in the GNU only against the compiled `MRANGE`; whether the
  Professional checks it against `mrange` was not measured (BACKLOG C95).
