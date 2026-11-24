# message — developer

- Record `MESSAGE` (INTEGER, `no_index`; `database.cc` after `ZIP`, enum
  after `ZIP` / `CONTROL_ZIP`).
- `top()` (`top.cc`), right after the `TOCHNOG_VERSION` record is set
  (after `input()`, so the input record is known): default `-YES`; the
  record, or `getenv("MESSAGE")` / `getenv("message")` equal to `-no` /
  `-NO`, gives `-NO`. The executable path is `readlink("/proc/self/exe")`
  (fallback "tochnog"), the build date the `TOCHNOG_VERSION` record.
- Professional 02-08-2026 probe: `MESSAGE=-NO`, `MESSAGE=-no`,
  `message=-no` and the record suppress it; `TOCHNOG_MESSAGE` does not.
- Pending: the Professional "suppress messages" may also cover other
  messages; only the introduction is gated.
