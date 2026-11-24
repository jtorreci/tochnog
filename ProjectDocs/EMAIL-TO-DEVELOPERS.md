<!--
  DRAFT — a personal call from a former Tochnog developer (2003) back to the
  original Tochnog GNU team.

  Recipients: the original Tochnog GNU authors — the 2003-era development team
  and the sf.net / GNU fork lineage. Dennis Roddeman is deliberately NOT
  addressed: he is associated with Tochnog Professional, which appears here only
  as the behavioural reference the GNU line is measured against, never as a
  recipient and never as someone being asked to fix anything.

  Before sending:
    - confirm the recipients' current addresses;
    - state your own 2003 role and period explicitly (the body leaves a
      placeholder for it);
    - replace [your name] and [your email];
    - attach the convergence report (ProjectDocs/report/convergence-report.pdf).

  Every technical claim is grounded in the repository documents
  (ProjectDocs/CONVERGENCE-CYCLE-REPORT.md, SEGUIMIENTO-CONVERGENCIA.md and the
  CORPUS-*/DIAG-* notes) and was measured against the reference binaries named
  in the attached report. This email is a human call; the report carries the
  detail.
-->

# Coming back to Tochnog — and a call to work together again

**To:** the original Tochnog GNU authors (the 2003-era team / sf.net lineage) — [addresses to be confirmed]
**From:** [your name] — [your email]
**Subject:** It has been a long time — coming back to Tochnog, with something to show

---

Dear friends and colleagues,

It is a long time since we enrolled in the Tochnog development, and it is a
genuine pleasure to write to you again. I was part of the team back in 2003 —
[state your role and period here] — and the project never really left my mind.
Recently I picked the GNU line back up, and I would love to tell you what has
happened since, and to ask whether you would like to be part of it again.

**What I have been doing.** I have been working on the GNU fork — the last
sf.net sources — with a simple goal: to understand where it behaves differently
from Tochnog Professional, and to make the same input give the same answer on
both. The two versions had drifted apart over the years, and some of that drift
turned out to be real defects in the GNU line. We measured each one against the
Professional binary before touching anything, fixed the GNU side where it was
wrong, and checked the result by running the full test corpora before and
after. Where a difference was deliberate, or simply an old test, we documented
it instead of "fixing" it.

**What we found (briefly).** A handful of genuine defects, all fixed and
validated: a band-width bug in the Bi-CG direct-LU retry that silently
corrupted coupled consolidation results; a missing `else` in the local
refinement dispatch; spurious `node_dof` requests that aborted mesh-only
models; an inverted sign in the groundflow velocity; and a record-length limit
that made a documented camclay form unreachable. We also recovered a validation
suite whose inputs had been deleted, by reconstructing them from the archived
`.dbs` dumps. The current state is healthy: the in-repo sfnet suite is at
187/197 and the reconstructed validation suite at 233/233 — 420 of 430 combined
(97.7%) — with the internal build suite 16/16. What remains is stale or
malformed test inputs, not source defects. The full account is in the attached
report, and the per-test evidence is in the repository.

**The call.** This is the real reason I am writing. I would love for you to:

1. **Try it.** If you are curious, run your old inputs on the GNU fork with the
   fixes. I can send the repository, a tarball, or the patches — whatever is
   easiest. I would genuinely value your reaction.
2. **Keep building it with us.** The project has life and a direction. There is
   a post-convergence plan — solver robustness, reproducible corpora, and new
   capabilities — and it would mean a great deal to develop it together again.
   It is your code; we are only its current caretakers.
3. **Help us find others.** Do you know colleagues who might be interested —
   from the original days, or people who would enjoy this? I would like to form
   a group and add new capabilities to Tochnog, and a few more hands (or just a
   few more memories) would make a real difference.
4. **Collaborate at whatever level you wish.** Testing, reviewing, coding, or
   simply telling us where we are wrong — all of it is welcome. There is no
   commitment implied; a quick reply is already a gift.

I know a lot of time has passed, and I do not want to impose. But if any part of
this resonates, I would be very happy to talk. Thank you for the original code,
and for having kept Tochnog alive.

With warm regards,
[your name]
