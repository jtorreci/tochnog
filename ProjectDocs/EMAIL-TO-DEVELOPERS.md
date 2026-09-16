<!--
  DRAFT — intended recipients: the original Tochnog GNU authors (the sf.net /
  GNU fork lineage). Dennis Roddeman is deliberately NOT addressed: he is
  associated with Tochnog Professional, and the Professional-side findings below
  are included as context only, never as questions or requests to the GNU
  authors.
  Before sending: confirm the GNU maintainer address, replace [your name], and
  review the asks. This is a report, not a bug-fix request; the tone is
  evidence-based and non-accusatory. Every claim below is grounded in the
  repository documents (ProjectDocs/CONVERGENCE-CYCLE-REPORT.md and its sources)
  and was measured against the reference binaries named in the text.
-->

# Draft email — Tochnog GNU: fork findings and the corpus recovery

**To:** the original Tochnog GNU authors (sf.net / GNU fork lineage) — [address to be confirmed]
**Subject:** Findings on the GNU Tochnog fork: five defects, a `.dbs` input recovery, and a few questions

---

Dear Tochnog GNU authors,

I have been working on a convergence effort between the Tochnog GNU fork (the
sfnet 2014 sources) and Tochnog Professional. The goal was to understand where
the two versions differ in behaviour and, where the difference is a defect in the
GNU fork, to fix it and prove the fix against a reference binary. I would like to
report what we found on the GNU side, ask you to sanity-check a few conclusions,
and ask whether any of the GNU-fork fixes is worth upstreaming.

## Context

The work ran against two reproducible corpora: the in-repo sfnet 197 test set and
a reconstructed 233-test validation suite. The measurements are documented per
test in the repository (`ProjectDocs/SEGUIMIENTO-CONVERGENCIA.md` and the
`CORPUS-*`/`DIAG-*` notes). A Tochnog Professional binary (build of 02-08-2026,
with the 25-10-2023 build for older paths) was used as the behavioural reference,
together with a build of the reference GNU 2014 sources.

## Method

Each candidate change started from a measured symptom and a measured root cause.
A fix was accepted only after a full corpus A/B showing that the per-test
classification change was exactly the intended set with zero regressions. Where a
hypothesis could not be confirmed against a reference binary, it was recorded as
refuted rather than implemented; several scheduled work units were dropped this
way.

## Defects in the GNU fork (fixed and validated)

These are the changes that look like genuine defects in the fork rather than
deliberate design choices. Each was validated against the corpora.

- **Bi-CG direct-LU retry band width** (`so.cc`, `solve()`, the `bicg_solver`
  band computation). The band was sized `diff*nprinc` instead of
  `(1+diff)*nprinc`; two nodes separated by `diff` node numbers span
  `diff*nprinc + (nprinc-1)` consecutive local indices, so the first-to-last dof
  coupling of adjacent nodes wrote outside `mat`/`ldmat` and silently corrupted
  the `dgbsv` retry. This produced wrong results on coupled consolidation
  systems (e.g. `examp2` gave `-0.7206` instead of the exact `-0.634519741323`,
  which matches the Professional's `-0.6345197413352`). The Cuthill-McKee branch
  already used the correct bound; only the Bi-CG branch omitted the `+1`.

- **`refine_l.cc` element-shape dispatch** (`refine_locally()`). A standalone
  `if (name==-BAR2)` followed by a chain starting at `if (name==-BAR3)` meant a
  `-bar2` element fell through to `assert(name==-TET10)` and aborted. Adding the
  missing `else` fixed it. The same function handles `-BAR2`/`-BAR3` correctly
  further down, which is how we identified it as a missing `else`.

- **Spurious `node_dof` requests** (two sites). A leftover debug
  `db_dbl(NODE_DOF, 1, ...)` in `top.cc` aborted every mesh-only model, and a
  redundant `|| band_solver` in `so.cc solve()` made the eigen task read/write
  `NODE_DOF[VERSION_NEW]` after that version had been deleted.

- **Groundflow velocity sign** (`groundfl.cc`, `groundflow()`). The equation was
  `gvel = -k grad(h)`; the Professional manual (8.2) defines
  `v_i^g = k_ij dh/dx_j`, and the reference binary on `condif5` returns
  `gvelx = +10` for a unit head gradient. The inverted sign reversed the
  convection direction in coupled `-groundflow -condif` models.

- **`group_materi_plasti_camclay` record length** (`database.cc`). The maximum
  stride was 3, which rejected the documented legacy 4-value form at parse time
  and made the `plasti.cc` legacy branch (`if (ldum>=4)`) unreachable. The record
  is variable-length 3..4.

## Divergences we explained rather than fixed

- Several archived targets were met only by the old Bi-CG loop, which reported
  success on three non-convergence exits. With the honest residual criteria the
  same models correctly report divergence. We treated this as an improvement and
  recalibrated the stale targets rather than restoring the dishonest exit.

- A number of legacy test inputs relied on the historic GNU default
  `options_inertia = -yes`. The Professional default is `-no` (manual 6.785), and
  the fork aligned with it. The affected models needed an explicit
  `inertia_apply -yes`; in several cases the Professional's own copies of the
  same tests already carry that record. The `-no` default is a deliberate
  alignment and was not reverted.

## Context from the Professional side (offered for information, not as asks)

We also measured two behaviours on the reference binary. We record them here only
as the context that motivated the alignment; they are not questions to the GNU
authors.

- The reference uses per-operation quadrature: interior Gauss for the element
  integral (2×2/2×2×2 quad4/hex8, 3×3/3×3×3 quad9/hex27) with a row-sum lumped
  mass, whereas the fork used Lobatto. The manual clause "in case inertia is
  applied then lobatto" (6.623) appears to describe the mass, which is row-sum
  lumped and quadrature-invariant; the reference binary contradicts it for the
  stiffness.
- The reference's nodal stress recovery is a separate per-element least-squares
  fit plus an arithmetic mean, not tied to the quadrature.

## The `.dbs` recovery method

Part of our validation suite had been deleted before the effort, and the
directory is not under version control. We recovered the inputs from the archived
`.dbs` database dumps: the dumps mix the parsed input records with run-produced
records, so we kept the input records and dropped the derived ones, and recovered
the one lost scalar (`number_of_integration_points`) from the archived record
lengths:

```
nuknwn = len(dof_label) + 1
number_of_integration_points = len(element_dof) / nuknwn
```

The derivation reproduces the known value in 11/11 cases where the true value is
known, and a fresh K-probe agrees 72/72. One useful corollary: `inertia_apply` is
never echoed in a `.dbs` (0 of 430 dumps), so its absence from a dump is not
evidence that the input lacked it.

## Current status

- sfnet 197: **187 PASS / 9 RUNFAIL / 1 PARSE**.
- validation-suite 233: **229 PASS / 4 RUNFAIL**.
- Combined: **416/430 (96.7%)**.
- Internal build suite: **16/16**.

The remaining failures are stale or malformed inputs, out-of-scope model or
infrastructure limits (PETSc, a long mesh-mapping test), two reconstruction
artifacts, and a small number of post-processing divergences. The full triage is
in `ProjectDocs/BACKLOG-OPEN-ISSUES.md`.

## Asks

1. **Upstream-worthiness of the unambiguous GNU defects.** The Bi-CG direct-LU
   retry band width and the `refine_l.cc` missing `else` look like unambiguous
   defects. Would a patch for either be welcome in the GNU sources? We are happy
   to prepare them in the project's style.

2. **Semantic confirmation on the GNU side.** Three of the fixes are semantic or
   documentation questions rather than clear-cut bugs, and we would rather
   confirm the intended behaviour with you before assuming ours is the right
   reading:
   - `groundflow_velocity`: is `v = +k grad(h)` the intended convention for the
     GNU fork as well, or was the inverted sign deliberate?
   - `group_materi_plasti_camclay`: is the legacy 4-value form intended to remain
     reachable, or is the 3-value form the canonical arity?
   - The spurious `node_dof` requests in `top.cc` and `so.cc`: can you confirm
     these are leftover debug/guard artifacts rather than load-bearing?

3. **Open questions.**
   - `options_inertia`: was the historic GNU default `-yes` intentional, or a
     divergence from the documented Professional `-no`? Several legacy test
     inputs depend on the old default.
   - `control_eigen` / `node_eigen`: we found no Professional equivalent (neither
     reference binary contains the strings). Is the eigen path intended as a
     GNU-only extension? We restored the original GNU semantics and have no
     reference to converge to.
   - Is there a canonical archived input set for the sf.net test suite, or a
     preferred way to keep the test inputs reproducible? We lost part of our
     validation inputs and had to reconstruct them from `.dbs` dumps; any
     archived originals would help us validate the reconstruction.

Thank you for the original code and for keeping the GNU fork available; the
convergence work would not have been possible without a well-defined behavioural
target.

Best regards,
[your name]
