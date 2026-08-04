# check_data

## Description

`check_data -yes` verifies the integrity of the database: it reports data items
that are defined but require another item (`data_required`) that is not present
and active at the same index. It is useful to detect incomplete or inconsistent
configurations before or during a run.

The check is complementary to the strong validation already done by `check()`
(e.g. `CHECK_COMBINATION`); it only reports items whose required companion is
missing.

## Usage

Place it in the data part, as a keyword line:

```
check_data -yes
```

## Parameters

| Parameter | Meaning                                                        |
|-----------|----------------------------------------------------------------|
| `-yes`    | Enable the database-integrity check.                           |
| `-no`     | Disable the check (default).                                   |

## Example

Minimal input that enables the integrity check:

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
check_data
   -yes
```

Run it; for every item that requires a missing companion the log will show a
line like
`Warning: data item <item> (index <n>) requires <other> which is not specified.`
and finish with a summary `check_data: <n> data item(s) have a missing
required item.` or `check_data: no missing required data items.`
