# refine_locally_bar2

## Description

`control_mesh_refine_locally` (local h-refinement driven by an element
error indicator) now works for models built from `-bar2` elements.

Before this fix, a local-refinement step on a `-bar2` model aborted the
calculation:

```
tochnog: refine_l.cc:231: void refine_locally(): Assertion `name==-TET10' failed.
```

The element-shape dispatch in `refine_l.cc` treated `-bar2` with a
standalone `if`, so a bar element fell through to the final branch that
asserts a tetrahedral element. The dispatch is now a proper `else if`
chain. `-bar2` and `-bar3` were already handled by the element-shape
dispatch further down the same routine, so this was a missing `else`, not
a missing feature.

## Usage

No input change is required. A model such as the sfnet `refine5` test:

```
node 1  0
node 2  1
element 1  -bar2 1 2
...
control_mesh_refine_locally 40  10.
```

now refines and completes instead of aborting.

## Verification

The sfnet `refine5` test (1D heat transfer + convection, stationary) runs
to completion with `rc=0` and reaches its target
(`post_point_dof 0 -temp = 1.0`).
