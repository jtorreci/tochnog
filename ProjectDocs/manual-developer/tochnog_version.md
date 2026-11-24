# tochnog_version — developer notes

Sprint 12 lot 3 (manual Professional 6.1091). Enum `TOCHNOG_VERSION`,
INTEGER, data_length 3 (day, month, year), INDEXED (the record keeps
the standard record index, unlike the solver_* family).

## Implementation

Written ONCE at the start of `top()`, right after `TIME_AT_START`,
parsed from the compiler's `__DATE__` string ("Mmm dd yyyy") with a
month-name table, so the record always matches the binary:

```c
if ( sscanf( __DATE__, "%15s %ld %ld", month_str, &day, &year )==3 ) ...
db( TOCHNOG_VERSION, 0, version, ddum, 3, VERSION_NORMAL, PUT );
```

The user cannot meaningfully write it (the program PUT overwrites at
start); its purpose is being queryable - `target_item -tochnog_version
0 <0|1|2>` (day/month/year) or the `.dbs` dump. Test `tslv_ver`
targets the year.

## Notes

- The historical GNU banner in `date.cc` ("Wed Jan 31 14:36:59 GMT
  2001") stays as is; the RECORD reflects the CURRENT build date.
