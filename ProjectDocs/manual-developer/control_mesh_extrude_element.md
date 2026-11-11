# control_mesh_extrude_element (F5 T4 U9b, U9c)

## Files and functions

- `tochnog.h` / `tochnog-mod.h` — enum values `CONTROL_MESH_EXTRUDE_ELEMENT`
  and `PRISM12` appended before `LAST_DUMMY`.
- `database.cc` — `control_mesh_extrude_element`: INTEGER, length 1,
  CONTROL, `data_required` = `CONTROL_MESH_EXTRUDE`. The word `prism12`
  names the 12-node prism ([prism12](prism12.md), F5 T4 U9c).
- `mesh.cc` — `mesh_extrude(icontrol)`:
  - `mesh_extrude_source(name, nnol)`: the 2D elements the routine turns
    into 3D elements (bar2, tria3, quad4, quad9 and now tria6); it replaces
    the three copies of that test (attached nodes, contact-spring groups,
    element loop).
  - The record is read after `_direction` into `tria6_to`: default
    `-PRISM18`; any value other than `-PRISM18` / `-PRISM12` stops with
    `must be -prism18 or -prism12`.
  - `has_quad9` / `has_tria6` mark the source elements. With both and
    `tria6_to == -PRISM12`, `tria6_to` becomes `-PRISM18` and a warning is
    printed (owner decision 2026-11-10, = Professional: conforming with
    the hex27). `has_mid = has_quad9 || (has_tria6 && tria6_to ==
    -PRISM18)`, so
    every layer has planes `b = l*ppl*nnod`, `m = (l*ppl+1)*nnod`,
    `t = (l+1)*ppl*nnod`; the prism18 is `nodes+b`, `nodes+m`, `nodes+t` in
    the tria6 order (the hex27 pattern).
  - `-PRISM12` (U9c): `nodes+b`, `nodes+t` in the tria6 order (the
    prism6 pattern); it only happens without a quad9.

## Professional probes (02-08-2026)

- Two tria6, z 0 1 3, n 2 1 1, no record: `-prism18` elements numbered
  `l*nelem+e`, nodes `p*nnod+n` with mid planes (7 planes of 9 nodes):
  the GNU mesh is identical (test `control_mesh_extrude_element_default1`).
- `-prism12`: 12-node prisms on the base and top planes (no mid plane),
  numbered as the GNU (`control_mesh_extrude_element_prism12_1`), also
  with new groups and a dropped layer (`_group12_1`) and with a quad4
  (`_mixed12_1`: the tria6 + quad4 mesh is the Professional one with or
  without the record, U9b probe).
- tria6 + quad9 with `-prism12`: the Professional writes `-prism18`
  (mid planes for all), also for a tria6 that only shares a corner node
  with the quad9; the GNU does the same since 2026-11-10.
- `-hex8`, `-prism15`, `-tria6`, `-yes`: silently `-prism18`.
- A quad4 + tria3 mesh with `-prism18`: hex8 and prism6, as without the
  record.
- tria6 + quad4 with `-prism18`: `-prism12` elements, no mid planes, and
  the calculation stops (`element_intpnt_iso_coord`).
- The record at an index without `control_mesh_extrude`: input error (the
  GNU `data_required` does the same).

## Divergences

- EXTRUDE-ELEMENT-NAME (class C): other names are an input error.
- `-prism12` with a quad9: no divergence. The U9c class D row
  (EXTRUDE-ELEMENT-PRISM12-QUAD9, GNU prism12 next to a hex27, not
  conforming) was removed by the owner decision of 2026-11-10: the mesh
  must be conforming, so every tria6 becomes a prism18, with a warning.
- EXTRUDE-ELEMENT-PRISM12 (class B until U9c): resolved, `-prism12`
  generates the Professional mesh.
- EXTRUDE-ELEMENT-MIXED (class D): tria6 with a linear element gives
  `-prism18` (manual) instead of the Professional `-prism12`.

## Tests

`control_mesh_extrude_element_default1`, `_prism18_1`, `_quad4_1`,
`_prism12_1`, `_group12_1`, `_mixed12_1`, `_quad9_1` (build_safe
`ext_check`), `_error1` (`expected_error_ok`). RED: the binary before U9b
rejects the record ("I do not know") and leaves the tria6 unextruded; the
binary before U9c stops on `-prism12` with "is not available yet" (the
U9b `_error2` test, removed).
