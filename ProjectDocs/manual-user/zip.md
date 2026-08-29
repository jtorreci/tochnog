# zip

## Description

With `zip -yes` all output files - `*flavia*`, `*msh`, `vtk*`,
`*.plt` and `*.dbs` - are compressed with the gzip program at the END
of the calculation (manual Professional 6.1095). Convenient in large
calculations with lots of output where the results are used later:
the plain files are replaced by their `.gz` counterparts and roughly
an order of magnitude of disk space is saved.

The gzip program must be installed; if the compression fails a warning
is printed and the calculation still finishes normally. Only files
that exist are compressed.

## Syntax

```
zip -yes
```

## Example

The `tslv_zip` test of the validation suite runs with `zip -yes` in an
isolated directory and verifies that `tslv_zip.dbs.gz` exists after
the run.

## Note

The per-control `control_zip` record of the Professional (zipping at
each control step) is not implemented; the GNU zips once at the end.
