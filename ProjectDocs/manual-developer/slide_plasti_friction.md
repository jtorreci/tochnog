# slide_plasti_friction

## Where

- Registration: `database.cc` (slide family block): SLIDE_PLASTI_FRICTION
  (DOUBLE x2, phi c), SLIDE_PLASTI_TENSION (x1), SLIDE_STIFFNESS (x2),
  SLIDE_PLASTI_RESIDUAL_STIFFNESS (x2), CONTROL_SLIDE_PLASTI_APPLY,
  CONTROL_SLIDE_STIFFNESS_APPLY (INTEGER x1, class CONTROL). Enums in
  `tochnog.h` + `tochnog-mod.h` (in sync).
- Legacy machinery: `slide.cc` (`slide()`, called from `top.cc` in the
  equilibrium loop after the element assembly): penalty constraint of
  the normal velocity + tangential friction `slide_friction * Fn`
  (Fn = assembled node_rhside · normal).

## Diagnosis (blockers of slide1/slide4)

With the Professional slide1.dat inputs (bounda on `-disy`/`-disx` at
the top, `slide_geometry 10 -bottom` on the bottom edge, E=1):
- the run completes but the block is NOT supported by the slide plane:
  the bottom nodes move down with the top (−1e-2), the element does
  not compress and the top reaction is ≈ 0 (target +0.01 / −0.5773e-2
  fail);
- even with explicit `node_slide 1 10`/`node_slide 2 10` membership
  records the constraint has no effect, so the failure is NOT the
  geometry membership test;
- the legacy penalty is added to `node_lhside[vel_indx+...]` of the
  velocity dofs; with the displacement-driven Professional inputs the
  term never produces a normal reaction — the mechanism needs to be
  re-derived against the Professional model (elastic kn/kt springs +
  Mohr-Coulomb cap c+Fn tan(phi), per node, on the TOTAL
  displacements), following the pattern of support_edge_normal
  (area.cc).

## Pending

- Rewrite `slide()` with the Professional elastic-plastic slide law
  (slide_stiffness kn/kt + friction cap + tension cut-off), the
  `node_slide_force`/`node_slide_f`/`node_slide_direction` output
  records and the two control gates.
