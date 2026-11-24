# zip — developer notes

Sprint 12 lot 3 (manual Professional 6.1095). Enum `ZIP`, INTEGER,
`no_index = 1`.

## Implementation

In `miscel.cc::exit_tn()`, AFTER the target checks (so target
evaluation reads the plain .dbs) and BEFORE `db_close()` (the record
lives in the database):

```c
db( ZIP, 0, &zip, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
if ( zip==-YES )
  system( "for f in *flavia* *msh vtk* *.plt *.dbs ; do "
          "[ -f \"$f\" ] && gzip -f \"$f\" ; done >/dev/null 2>&1" );
```

The per-file existence test is REQUIRED: a bare
`gzip -f *flavia* ...` returns non-zero when any glob has no match
(the shell passes the unexpanded pattern and gzip errors on it) - the
warning fired even when the compression succeeded (measured). With the
loop, the exit status reflects real failures only.

## Testing

`tslv_zip` runs in an ISOLATED directory (`/tmp/tslv_zip_isolated`,
created and removed by `build_safe.sh`): the gzip glob would compress
EVERY *.dbs of the shared suite directory. The check verifies
`tslv_zip.dbs.gz` exists and no gzip warning was printed.

## Notes

- The Professional's per-control `control_zip` (zipping at each control
step, overruling `zip` per index) is not implemented; the GNU zips
once at the end.
