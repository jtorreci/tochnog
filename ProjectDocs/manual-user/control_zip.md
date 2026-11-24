# control_zip

## Description

`control_zip index switch` (manual Professional 6.394): with `-yes` all
`*flavia*`, `*msh`, `vtk*`, `*.plt` and `*.dbs` files of the working
directory are gzipped at the end of control index `index` (the gzip
program must be installed). Later indices write new output files; the
global [zip](zip.md) gzips at the end of the calculation.

## Usage

```
control_zip index switch
```

Test: `validation-suite/test-2026/control_zip1` (run in an isolated
directory: every `*.dbs` of the working directory is gzipped).

## Differences with the Professional version

None measured (Professional 02-08-2026: the GiD files of index 0 are
gzipped, index 1 writes new ones).
