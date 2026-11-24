# group_materi_umat_pardiso_decompose

## Description

In the Professional (manual 6.758) this switch tells the PARDISO solver to
decompose the system matrix every iteration (`-yes`, default) or only
once (`-no`) when a umat routine is present. The GNU has no PARDISO
solver: the record is accepted and has no effect.

## Syntax

```
group_materi_umat_pardiso_decompose index -yes|-no
```

Test: `validation-suite/test-2026/group_materi_umat_pardiso_decompose1`.

## Differences with the Professional version

- No effect in the GNU (DIVERGENCES UMAT-PARDISO-DECOMPOSE, class A):
  the GNU solvers decompose as they always do.
