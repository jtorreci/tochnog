# group_materi_stokes (developer)

## Where it is implemented

- `database.cc`: record `GROUP_MATERI_STOKES` (`group_materi_stokes`),
  INTEGER, length 1, class MATERI, requires `GROUP_TYPE`.
- `check.cc`: the record requires the `materi_velocity` unknown.
- `general.cc` `general()`: the only consumer.
  - The local `stokes` starts as `-YES` (the manual 6.752 default).
  - For `type==-MATERI` it is read with
    `db( GROUP_MATERI_STOKES, element_group, &stokes, ..., GET_IF_EXISTS )`.
  - In the convection loop, the `MATERI_VELOCITY` branch sets
    `conv_part = dens` only when `stokes!=-YES`. With Stokes flow
    `conv_part` stays 0, so the momentum convection term, its matrix
    entries and the artificial diffusion of the velocity all vanish
    (the Peclet number uses `conv_part`).
  - The other convected dofs (stress, strains, history variables,
    density) keep their `conv_part` whatever the switch.

## History (F4 T6f K21 WU2, 2026-10-12)

- Before this work unit, `stokes` started as `-NO` and was read at index
  0 for every group.
- The default is now `-YES`, as in manual 6.752 and the Professional
  02-08-2026 (`validation_1` without the record: velx 0.1465713).
- The record is now read per element group (owner decision), consistent
  with the other group records (`GROUP_AXISYMMETRIC` just above it).
- Tests: `test-2026/stokes_default1`, `stokes_group1` and
  `stokes_group2`, wired with `run_ok` in `scripts/build_safe.sh`.
- `test-2014/cmat_gate` was rewritten as a well-posed simple shear. The
  old ill-posed model failed its target with the new default; the old
  target was met by chance.

## Investigation notes (Professional 02-08-2026 probes)

- The density (0, 1 or 2) does not change the Professional Stokes result,
  so its `stokes -yes` against `convection_apply -no` gap is not a
  momentum term.
- In the GNU Stokes run the nodal stress differs from the convection
  `-no` run (sigxx -0.200808 against -0.200826), but the velocity is
  identical. In the Professional run the stress feeds back into the
  velocity.
- Conclusion (corrected 2026-10-13, F4 T6f C72, Engram 2775): the
  Professional does NOT convect the integration-point stress either (a
  stress pulse stays bit-frozen), and its Stokes value is unchanged with
  `global_element_dof_apply -no`, so the stress storage is not the cause.
  The residual is a convection-gated Professional effect whose mechanism
  is not identified (BACKLOG C72); it does not belong to this switch.
- Index probes: the result is always the index-0 value, to 1e-12. Probed
  cases: groups `-no/-yes` and `-yes/-no`; the record only on group 1;
  element 0 moved to group 1. The probes are in the session scratch;
  Engram 2765.

## External features / hardcoded parameters

- None. The switch values are the generic `-YES` / `-NO` enums.

## Pending refactorings

- `general()` declares `stokes` together with unrelated locals. A small
  `materi_stokes_switch(element_group)` helper would make the default
  explicit if a second consumer appears.
- The Stokes residual stays open under C72 (Professional-side effect,
  mechanism not identified).

## Differences with the Professional version

- **index of the record** (deliberate divergence, class A; DIVERGENCES
  `STOKES-IDX`, owner 2026-10-12). GNU: per element group (manual 6.752).
  Professional: index 0 for every group. Guard tests: `stokes_group1`
  (0.124105) and `stokes_group2` (-0.0570386). The Professional gives
  0.1465713 for both and fails them.
- **Stokes flow on an Eulerian mesh** (class B; BACKLOG C72).
  Professional `validation_1` without the record: 0.1465713; GNU:
  0.144159, equal to `convection_apply -no`. Neither program convects the
  integration-point stress; the Professional mechanism is not
  identified. The tolerance of
  `stokes_default1` (0.003) covers this gap.
- See also [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
