# DEV-A HANDOFF — incremental_driver

Sprint: 4-dev parallel convergence (GNU tochnog vs Tochnog Professional 25-10-2023).
Branch: `dev/incdrv` (from 19e995a). Worktree: /tmp/opencode/tn-dev-incdrv.

## State (updated: ISA parse-only registration)
- [x] Worktree verified clean @ 19e995a; HANDOFF created first (c12ad22).
- [x] Registered `materi_strain_isa_c` and `materi_strain_isa_eacc` for
      Professional-compatible parsing and DOF storage only. This does **not**
      implement ISA/Niemunis constitutive behavior; the incremental driver
      deliberately leaves both records unused.
- [x] DECISIONES.md + SEGUIMIENTO notes read (updated_area context, family diagnosis).
- [x] Pro algorithm enumerated from binary symbols (NOT just the old asm dumps —
      full re-disassembly with the Pro name-table extracted from initialise():
      /tmp/incdrv_dis_v2/pro_name_map.json maps Pro db index -> name).
- [x] Pro A/B probes run fresh: oedo, triax drained, triax undrained, syntax,
      + 9 parameter-fit probes (fitA-C, d1-d5). All rc=0, md5-stable.
- [x] Precise spec written: ProjectDocs/manual-developer/incremental_driver.md.
- [ ] Implement: incremental.cc (driver core) + top.cc (hook) + db records.
- [ ] Suite 16/16 green; corpus targets vs Pro values.

## Key findings (all verified, details in the spec doc)
1. Pro hook: if incremental_driver active -> run driver, write result file, exit_tn()
   (normal FEM loop never runs). GNU must hook top.cc accordingly.
2. Driver = single-material-point simulator over combined experiments; NO FEM mesh.
3. Displacement-controlled drained tests are EXACTLY: small-strain eps from
   displacement/length0, elastic trial (E,nu), MC-direct cut (p-constant,
   Lode-preserving at psi=0; measured mean-shift direction at psi>0), 1000 rows.
4. Undrained: p = evol/capacity EXACT; lateral equilibrium on sigma_law + p;
   MC surface on the law (effective) stress. Measured invariants to 1e-9.
5. Result file: 2 header lines + per-combined-experiment blocks, 15 scientific
   columns, absolute time, combined experiments re-run from the initial state.
6. Numerical match strategy: replicate the converged semantics; Pro Newton noise
   ~1e-7..1e-8 relative is part of its output (corpus tolerance 1e-2).

## Evidence paths
- Spec: ProjectDocs/manual-developer/incremental_driver.md
- ISA parse-only boundary: ProjectDocs/manual-developer/incremental_driver_isa_parse_only.md
- Disassembly: /tmp/incdrv_dis_v2/ (entry/orch/deformation/mcd/gflow/calcdef/deften
  + *_annot.txt annotated; pro_name_map.json = Pro record index->name)
- Pro reference results: /tmp/incdrv_pro/{oedo,syntax}/, /tmp/incdrv_probes/
  {orig_triax,rep_undr,fitA,fitB,fitC,d1..d5,triax_s0,triax_n2,triax_nt*}
- Prior killed-attempt artifacts superseded by the above: /tmp/incdrv_pro/* (triax
  and undr dirs contain EMPTY results from a transient Pro glitch — use my probes).

## Commits
- c12ad22 docs(incdrv): create DEV-A HANDOFF (milestone 0: sprint start)
- (next) docs(incdrv): measured Pro algorithm spec
