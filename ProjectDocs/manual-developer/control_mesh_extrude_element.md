# control_mesh_extrude_element (F5 T4 U9b)

## Files and functions

- `tochnog.h` / `tochnog-mod.h` — enum values `CONTROL_MESH_EXTRUDE_ELEMENT`
  and `PRISM12` appended before `LAST_DUMMY`.
- `database.cc` — `control_mesh_extrude_element`: INTEGER, length 1,
  CONTROL, `data_required` = `CONTROL_MESH_EXTRUDE`. The word `prism12` is
  registered (name only) so that the record can refuse it with a clear
  message; the element is F5 T4 U9c.
- `mesh.cc` — `mesh_extrude(icontrol)`:
  - `mesh_extrude_source(name, nnol)`: the 2D elements the routine turns
    into 3D elements (bar2, tria3, quad4, quad9 and now tria6); it replaces
    the three copies of that test (attached nodes, contact-spring groups,
    element loop).
  - The record is read after `_direction`: default `-PRISM18`; `-PRISM12`
    stops with the pending message, any other value with `must be -prism18
    or -prism12`.
  - A tria6 sets the mid-plane flag `has_mid` (formerly `has_quad9`), so
    every layer has planes `b = l*ppl*nnod`, `m = (l*ppl+1)*nnod`,
    `t = (l+1)*ppl*nnod`; the prism18 is `nodes+b`, `nodes+m`, `nodes+t` in
    the tria6 order (the hex27 pattern).

## Professional probes (02-08-2026)

- Two tria6, z 0 1 3, n 2 1 1, no record: `-prism18` elements numbered
  `l*nelem+e`, nodes `p*nnod+n` with mid planes (7 planes of 9 nodes):
  the GNU mesh is identical (test `control_mesh_extrude_element_default1`).
- `-prism12`: 12-node prisms on the base and top planes (no mid plane).
- `-hex8`, `-prism15`, `-tria6`, `-yes`: silently `-prism18`.
- A quad4 + tria3 mesh with `-prism18`: hex8 and prism6, as without the
  record.
- tria6 + quad4 with `-prism18`: `-prism12` elements, no mid planes, and
  the calculation stops (`element_intpnt_iso_coord`).
- The record at an index without `control_mesh_extrude`: input error (the
  GNU `data_required` does the same).

## Divergences

- EXTRUDE-ELEMENT-NAME (class C): other names are an input error.
- EXTRUDE-ELEMENT-PRISM12 (class B, pending U9c): `-prism12` refused.
- EXTRUDE-ELEMENT-MIXED (class D): tria6 with a linear element gives
  `-prism18` (manual) instead of the Professional `-prism12`.

## Tests

`control_mesh_extrude_element_default1`, `_prism18_1`, `_quad4_1`
(build_safe `ext_check` / `pv_case`), `_error1`, `_error2`
(`expected_error_ok`). RED: the binary before U9b rejects the record
("I do not know") and leaves the tria6 unextruded.
