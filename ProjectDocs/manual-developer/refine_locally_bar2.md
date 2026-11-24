# refine_locally_bar2 (developer)

Bug fix, not a new keyword. `control_mesh_refine_locally` aborted with
`rc=134` (SIGABRT) on a model whose elements are `-bar2`.

## Symptom

`refine5` (sfnet corpus, 1D heat + convection, 2 `-bar2` elements,
`control_mesh_refine_locally 40 10.`) aborts:

```
tochnog: refine_l.cc:231: void refine_locally(): Assertion `name==-TET10' failed.
```

The abort is deterministic and pre-existing: the reference GNU build from the
jan-2014 sfnet sources (F. Lorenzo, `/tmp/opencode/tn2014/build/tochnog`)
aborts at the same assertion, and `refine_l.cc` has not been touched since
`a2407e0` ("Initial Version"). The Bi-CG breakdown message that precedes the
abort on the current GNU is a red herring: the direct-LU retry succeeds and the
calculation would have continued.

## Root cause

`refine_l.cc`, `refine_locally()`, the element-shape dispatch that sets
`nside`/`nnol_side`:

```c
if ( name==-BAR2 ) {          // standalone `if` ...
  nside = 1;
  nnol_side = 2;
}
if ( name==-BAR3 ) {          // ... but the chain starts here
  ...
}
else if ( name==-TRIA3 ) { ... }
else if ( name==-TRIA6 ) { ... }
else if ( name==-TET4  ) { ... }
else {
  assert( name==-TET10 );     // BAR2 falls through to here and aborts
  ...
}
```

For `name == -BAR2` the first `if` runs, then the second `if` is evaluated
independently and is false, so control reaches the final `else` and the
`assert(name==-TET10)` fails. The `else` chain must start at the BAR3 branch.
The element-shape dispatch later in the same function
(`refine_l.cc:677-690`) handles `-BAR2` and `-BAR3` with a proper
`else if` chain, which confirms that BAR2 is a supported shape and the
missing `else` is the defect.

## Fix

`refine_l.cc:203`: `if ( name==-BAR3 )` -> `else if ( name==-BAR3 )`.

## Verification

- `refine5` rc=0, `post_point_dof 0 -temp = 1` (target `1.0` +/- `5e-2`).
- Reference GNU 2014 aborts at the same assert (rc=134); the fixed GNU passes.
- sfnet corpus 197: `184 PASS / 13 RUNFAIL / 0 PARSE` ->
  `186 PASS / 11 RUNFAIL / 0 PARSE`; per-test diff = exactly
  `refine5 RUNFAIL -> PASS` (plus `incnav3`, a separate corpus-input
  alignment, see `CORPUS-SOLVER-TIMEOUT.md`), 0 regressions.
- Internal suite `scripts/build_safe.sh`: 16/16.

## Pending

A second copy-paste defect in the same function, the side-node list inside the
`for (iside=...)` loop (`refine_l.cc:249`): `else if ( name==-BAR2 )` is
unreachable (the previous branch already matched `-BAR2`) and should be
`-BAR3`, so the `-BAR3` local-refinement path would use a stale `ilist`. It is
NOT exercised by any corpus test (`refine5` is the only
`control_mesh_refine_locally` model with bar elements, and it is `-bar2`) and is
left untouched to keep this fix minimal.
