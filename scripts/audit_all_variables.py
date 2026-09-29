#!/usr/bin/env python3
"""Tochnog GNU vs Tochnog Professional: all-variable audit.

Runs every test of a corpus on both binaries (one isolated work directory
per test and binary), compares the FINAL state variable by variable (dofs
are mapped by name through the `dof_label` record of each .dbs, never by
column; nodes are paired by number, or by coordinates when the meshes are
numbered differently), evaluates each test's own targets on both dumps,
and assigns every mismatching test to a triage cluster.

Subcommands (see scripts/audit_all_variables.sh for the full pipeline):
  run        --bin B --src DIR --out DIR [--convert pro|iter|pro+iter]
             [--timeout 45] [--jobs 7] [--only a,b] [--exclude a,b] [--ld LIB]
  compare    GNU_DIR PRO_DIR NAME            (one test, JSON on stdout)
  aggregate  RUNS                            (-> RUNS/summary.json)
  iter-effect RUNS RUNS_ITER OUT.json         (5-iteration diagnostic)
  clusters   RUNS [ITER_EFFECT.json]         (-> RUNS/clusters.json, cluster_table.json)
  stats      RUNS                            (markdown tables for the report)
  mism-list  RUNS CORPUS                     (mismatching tests, comma list)

Layout expected by aggregate: RUNS/<corpus>/{gnu,pro[,gnuconv]}/<test>/
with a classification.tsv per side. For GNU-dialect corpora `gnu` holds the
original input, `gnuconv` and `pro` the mechanically converted input, and the
comparison is gnuconv vs pro (both binaries run the SAME text).

Report: ProjectDocs/AUDIT-ALL-VARIABLES-2026.md.
"""
import argparse, collections, json, math, os, re, shutil, subprocess, sys, time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

# ===CONVERT===
import re

# (regex, replacement, tag). Anchored on the record name at line start or on
# a whole switch word, so data values are never touched.
RULES = [
    (r'^(\s*)bounda_unknown\b', r'\1bounda_dof', 'bounda_unknown'),
    (r'^(\s*)options_mesh\b', r'\1mesh', 'options_mesh'),
    (r'^(\s*)control_options_inertia\b', r'\1control_inertia_apply', 'control_options_inertia'),
    (r'^(\s*)options_inertia\b', r'\1inertia_apply', 'options_inertia'),
    (r'^(\s*)control_options_convection\b', r'\1control_convection_apply', 'control_options_convection'),
    (r'^(\s*)options_convection\b', r'\1convection_apply', 'options_convection'),
    (r'^(\s*)options_stabilization(\s+)-static\b', r'\1convection_stabilization\2-yes', 'options_stabilization'),
    (r'^(\s*)options_stabilization(\s+)-none\b', r'\1convection_stabilization\2-no', 'options_stabilization'),
    (r'^(\s*)options_stabilization\b', r'\1convection_stabilization', 'options_stabilization'),
    (r'^(\s*)options_processors\b', r'\1processors', 'options_processors'),
    (r'^(\s*)control_options_solver_bicg_error\b', r'\1control_solver_bicg_error', 'control_options_solver'),
    (r'^(\s*)control_options_solver\b(?!_)', r'\1control_solver', 'control_options_solver'),
    (r'^(\s*)options_solver\b', r'\1solver', 'options_solver'),
    (r'^(\s*)group_truss_young\b', r'\1group_truss_elasti_young', 'group_truss_young'),
    (r'(\s)-trussbeam\b', r'\1-truss_beam', 'trussbeam'),
    (r'^(\s*)control_mesh_generate_trussbeam\b', r'\1control_mesh_generate_truss_beam', 'trussbeam'),
    (r'^(\s*)groundflow_phreaticlevel\b', r'\1groundflow_phreatic_level', 'phreaticlevel'),
    (r'^(\s*)contact_geometry\b', r'\1contact_target_geometry', 'contact_geometry'),
    (r'^(\s*)group_materi_plasti_druckprag\b', r'\1group_materi_plasti_druck_prag', 'druckprag'),
    (r'^(\s*)group_materi_failure_cruching\b', r'\1group_materi_failure_crunching', 'cruching'),
    (r'^(\s*)group_materi_failure_voidfraction\b', r'\1group_materi_failure_void_fraction', 'voidfraction'),
    (r'^(\s*)control_mesh_refine_locally_unknown\b', r'\1control_mesh_refine_locally_dof', 'refine_locally_unknown'),
    (r'^(\s*)control_mesh_merge_epscoord\b', r'\1control_mesh_merge_eps_coord', 'merge_epscoord'),
    (r'^(\s*)control_mesh_delete_geometry_movenodes\b', r'\1control_mesh_delete_geometry_move_node', 'movenodes'),
    (r'^(\s*)control_print_unknownsrhside\b', r'\1control_print_dof_rhside', 'unknownsrhside'),
    (r'^(\s*)geometry_bounda_factor\b', r'\1geometry_factor', 'geometry_bounda_factor'),
    (r'^(\s*)group_materi_plasti_hypo_wolfersdorff\b', r'\1group_materi_plasti_hypo_wolffersdorff', 'wolfersdorff'),
    (r'^(\s*)group_materi_plasti_hypo_intergranularstrain\b', r'\1group_materi_plasti_hypo_strain_intergranular', 'intergranularstrain'),
    (r'^(\s*)group_materi_plasti_hypo_pressuredependentvoidratio\b', r'\1group_materi_plasti_hypo_void_ratio_linear', 'pressuredependentvoidratio'),
    (r'^(\s*)area_element_group_sequence_elementgroup\b', r'\1area_element_group_sequence_element_group', 'sequence_elementgroup'),
    (r'(\s)-sizetot\b', r'\1-size_tot', 'sizetot'),
    (r'^(\s*)contactspring', r'\1contact_spring', 'contactspring'),
    (r'^(\s*)group_contactspring_friction\b', r'\1group_contact_spring_plasti_friction', 'contactspring'),
    (r'^(\s*)group_contactspring_', r'\1group_contact_spring_', 'contactspring'),
    (r'(\s)-contactspring\b', r'\1-contact_spring', 'contactspring'),
    (r'^(\s*)force_element_volume', r'\1force_volume', 'force_element_volume'),
    (r'^(\s*)force_element_edge', r'\1force_edge', 'force_element_edge'),
    (r'^(\s*)group_groundflow_porosity\b', r'\1group_porosity', 'group_groundflow_porosity'),
    (r'^(\s*)group_materi_plasti_boundary', r'\1group_materi_plasti_bounda', 'plasti_boundary'),
    (r'^(\s*)group_materi_plasti_cap(\s)', r'\1group_materi_plasti_cap2\2', 'plasti_cap->cap2'),
    # targets only drive the GNU exit code; the Professional rejects many
    # GNU dof basenames used in them (-sph, -pc, -hisv10, ...)
    (r'^(\s*)target_(item|value)\b.*$', r'\1', 'targets_dropped'),
    (r'(\s)-sizedev\b', r'\1-size_dev', 'sizedev'),
    # --- approximations (flagged with '!') ---
    # legacy GNU Mohr-Coulomb (same phi c phiflow layout) -> the
    # Professional/GNU group_materi_plasti_mohr_coul: switches the GNU code
    # path (plasti.cc MOHRCOUL vs MOHR_COUL), checked by gnu-vs-gnuconv
    (r'^(\s*)group_materi_plasti_mohrcoul(\s)', r'\1group_materi_plasti_mohr_coul\2', '!mohrcoul->mohr_coul'),
    (r'(\s)-updated_without_rotation\b', r'\1-updated_linear', '!updated_without_rotation'),
    (r'^(\s*)number_of_integration_points\b.*$', r'\1', '!number_of_integration_points'),
]


def convert(text):
    notes = []
    for pat, rep, tag in RULES:
        new, n = re.subn(pat, rep, text, flags=re.M)
        if n:
            if tag not in notes:
                notes.append(tag)
            text = new
    # the Professional refuses materi_velocity_integrated next to
    # materi_displacement (both integrate the velocity)
    ini = text.split('end_initia')[0]
    if re.search(r'^\s*materi_displacement\s*$', ini, flags=re.M) and \
            re.search(r'^\s*materi_velocity_integrated\s*$', ini, flags=re.M):
        # the drop changes the dof count: explicit node_dof input records
        # are only removable when they are all zero (= the default state)
        nd = re.findall(r'^\s*node_dof\s+(.*)$', text, flags=re.M)
        ok = True
        for rec in nd:
            t = rec.split()
            if not t or not re.match(r'^\d+$', t[0]):
                ok = False      # ranged/indexed forms: not handled
                break
            for v in t[1:]:
                try:
                    if float(v) != 0.0:
                        ok = False
                except ValueError:
                    ok = False
        if ok:
            text = re.sub(r'^(\s*)materi_velocity_integrated[ \t]*$', r'\1', text,
                          count=1, flags=re.M)
            text = re.sub(r'^\s*node_dof\s+.*$', '', text, flags=re.M)
            notes.append('!velocity_integrated_dropped')
    ini = text.split('end_initia')[0]
    # the Professional refuses an explicit 'convection_apply -no' next to
    # materi_displacement; -no is its default, and the GNU ignores
    # convection when materi_displacement is on (general.cc), so the
    # record is dropped: neutral in both binaries
    if re.search(r'^\s*materi_displacement\s*$', ini, flags=re.M):
        text, n = re.subn(r'^\s*convection_apply\s+-no\s*$', '', text, flags=re.M)
        if n:
            notes.append('convection_apply_no_dropped')
    # the Professional requires this post-processing request with pore
    # pressures (post-processing only: adds node_dof_calcul entries)
    if re.search(r'^\s*groundflow_pressure\s*$', ini, flags=re.M) and \
            not re.search(r'^\s*post_calcul\b', text, flags=re.M):
        text = re.sub(r'^(\s*)end_data', r'post_calcul -groundflow_pressure -total_pressure\n\1end_data',
                      text, count=1, flags=re.M)
        notes.append('post_calcul_total_pressure_added')
    # the GNU history-variable initia maps onto the per-model Professional
    # name only for di Prisco (documented alias, mdiprisco_hist)
    if 'group_materi_plasti_diprisco' in text and \
            re.search(r'^\s*materi_history_variables\b', text, flags=re.M):
        text = re.sub(r'^(\s*)materi_history_variables\b',
                      r'\1materi_plasti_diprisco_history', text, flags=re.M)
        notes.append('history_variables->diprisco_history')
    return text, ','.join(notes)


def add_iterations(text, n=5):
    """Diagnostic variant: force n equilibrium iterations in every
    control_timestep index (replaces existing control_timestep_iterations
    records). Used only to separate one-iteration lag artefacts from
    physics differences; never part of the main comparison."""
    body = re.sub(r'\([^)]*\)', lambda m: ' ' * len(m.group(0)), text, flags=re.S)
    idx = sorted(set(int(m.group(1)) for m in
                     re.finditer(r'^\s*control_timestep\s+(\d+)\s', body, flags=re.M)))
    text = re.sub(r'^\s*control_timestep_iterations\s+\d+\s+\S+.*$', '', text, flags=re.M)
    add = ''.join('control_timestep_iterations %d %d\n' % (i, n) for i in idx)
    text = re.sub(r'^(\s*)end_data', add + r'\1end_data', text, count=1, flags=re.M)
    return text, 'iterations%d' % n

# ===RUN===




def aux_files(src_dir, dat_text, all_files):
    """Files of the corpus directory referenced by name in the input."""
    words = set(re.findall(r'[A-Za-z0-9_.\-/]+', dat_text))
    out = []
    for f in all_files:
        if f in words and not f.endswith('.dat'):
            out.append(f)
    return out


def classify(rc, log, dbs_ok, timed_out):
    if timed_out:
        return 'TIMEOUT'
    low = log.lower()
    if 'i do not know' in low or "don't know what to do" in low or \
            'error in data part' in low or 'error in initialisation part' in low:
        return 'PARSE'
    if dbs_ok:
        return 'OK'
    return 'RUNFAIL'


def run_one(args, name, all_files):
    src = os.path.join(args.src, name + '.dat')
    wd = os.path.join(args.out, name)
    if os.path.isdir(wd):
        shutil.rmtree(wd)
    os.makedirs(wd)
    with open(src, errors='replace') as fh:
        text = fh.read()
    note = ''
    if args.convert in ('pro', 'pro+iter'):
        text, note = convert(text)
    if args.convert in ('iter', 'pro+iter'):
        text, n2 = add_iterations(text)
        note = (note + ',' + n2).strip(',')
    with open(os.path.join(wd, name + '.dat'), 'w') as fh:
        fh.write(text)
    for f in aux_files(args.src, text, all_files):
        try:
            sp = os.path.join(args.src, f)
            if os.path.isdir(sp):
                shutil.copytree(sp, os.path.join(wd, f))
            else:
                shutil.copy(sp, wd)
        except Exception:
            pass
    env = dict(os.environ, LC_ALL='C')
    if args.ld:
        env['LD_LIBRARY_PATH'] = args.ld
    t0 = time.time()
    timed_out = False
    with open(os.path.join(wd, 'run.log'), 'wb') as lf:
        try:
            p = subprocess.run([args.bin, name + '.dat'], cwd=wd, stdout=lf,
                               stderr=subprocess.STDOUT, env=env,
                               timeout=args.timeout, stdin=subprocess.DEVNULL)
            rc = p.returncode
        except subprocess.TimeoutExpired:
            rc, timed_out = 124, True
    dt = time.time() - t0
    with open(os.path.join(wd, 'run.log'), errors='replace') as fh:
        log = fh.read()
    dbs = os.path.join(wd, name + '.dbs')
    dbs_ok = os.path.isfile(dbs) and os.path.getsize(dbs) > 0 and \
        'end_data' in open(dbs, errors='replace').read()[-4000:]
    c = classify(rc, log, dbs_ok, timed_out)
    unknown = ''
    m = re.search(r"(?:I do not know|I don't know what to do with)\s*:\s*(\S+)", log)
    if m:
        unknown = m.group(1)
    elif 'initialisation part' in log:
        unknown = '<initia>'
    if not unknown:
        for line in log.splitlines():
            if re.search(r'error|abort|assert|segmentation|diverg', line, re.I):
                unknown = line.strip()[:90].replace('\t', ' ')
                break
    # large outputs we never compare: drop them to keep the scratch small
    for f in os.listdir(wd):
        if f.endswith(('.res', '.msh', '.vtk', '.frd')):
            try:
                os.remove(os.path.join(wd, f))
            except OSError:
                pass
    return name, c, rc, round(dt, 1), unknown, note


# ===COMPARE===

RTOL = 1e-4
AFLOOR = 1e-6
ABSFLOOR = 1e-10  # round-off floor for fields that are identically ~0
SIGN_FLOOR = 1e-3
SIGN_FLOOR_NOTE = 'sign flips are only reported above SIGN_FLOOR*scale'
NOINDEX = {'time_current', 'dof_label', 'time_old', 'dtime', 'time_start', 'post_calcul_label'}

NUM = re.compile(r'^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$')


def isnum(s):
    return bool(NUM.match(s)) or s in ('nan', '-nan', 'inf', '-inf')


def fnum(s):
    try:
        return float(s)
    except ValueError:
        return float('nan')


def parse_dbs(path):
    """-> dict name -> {index(int or None): [tokens]}"""
    recs = {}
    cur = None
    with open(path, errors='replace') as fh:
        for line in fh:
            toks = line.split()
            if not toks:
                continue
            if toks[0].startswith('('):
                continue
            if isnum(toks[0]) and cur is not None:
                cur.extend(toks)
                continue
            if isnum(toks[0]) or toks[0].startswith('-'):
                if cur is not None:
                    cur.extend(toks)
                continue
            name = toks[0]
            rest = toks[1:]
            idx = None
            if rest and name not in NOINDEX and re.match(r'^\d+$', rest[0]):
                idx = int(rest[0])
                rest = rest[1:]
            cur = list(rest)
            recs.setdefault(name, {})[idx] = cur
    return recs


def labels(recs):
    lab = recs.get('dof_label', {}).get(None)
    if not lab:
        return []
    return [t.lstrip('-') for t in lab]


def base_set(lab):
    """indices of labels that are not derivative slots of another label"""
    names = set(lab)
    out = []
    for i, l in enumerate(lab):
        if len(l) > 1 and l[0] in 'xyzt' and l[1:] in names:
            continue
        out.append(i)
    return out


def coords(recs):
    for key in ('node_start_refined', 'node'):
        d = recs.get(key)
        if d:
            return {k: [fnum(x) for x in v] for k, v in d.items() if k is not None}
    return {}


def pair_nodes(cg, cp):
    """-> (mapping gnu_node -> pro_node, method)"""
    common = set(cg) & set(cp)
    if not common:
        return {}, 'none'
    size = 0.0
    for c in list(cg.values()):
        for x in c:
            if not math.isnan(x):
                size = max(size, abs(x))
    tol = 1e-6 * max(size, 1.0)
    bad = 0
    for n in common:
        a, b = cg[n], cp[n]
        if len(a) != len(b) or any(abs(x - y) > tol for x, y in zip(a, b)):
            bad += 1
    if bad == 0 and len(cg) == len(cp):
        return {n: n for n in common}, 'number'
    # pair by coordinates
    key = lambda c: tuple(round(x / tol) if not math.isnan(x) else 0 for x in c)
    idx = {}
    for n, c in cp.items():
        idx.setdefault(key(c), n)
    m = {}
    for n, c in cg.items():
        k = key(c)
        if k in idx:
            m[n] = idx[k]
    if len(m) >= 0.9 * max(len(cg), 1):
        return m, 'coord(%d/%d)' % (len(m), len(cg))
    return {n: n for n in common}, 'number-mismatched-coords(%d bad)' % bad


def group_key(label):
    """tensor/vector family of a dof label: sigxy -> sig, velx -> vel"""
    k = label
    if len(label) > 3 and re.search(r'[xyz]{1,2}$', label):
        k = re.sub(r'[xyz]{1,2}$', '', label)
    # kinematic dofs share one magnitude (dis ~ vel*t), strains too
    if k in ('vel', 'veli', 'dis', 'acc'):
        return 'kin'
    if k.startswith('ep'):
        return 'eps'
    return k


def compare_pairs(pairs, fscale=0.0):
    """pairs: list of (where, g, p) -> summary dict. fscale: magnitude of
    the variable's tensor/vector family in this test (absolute floor)"""
    scale = 0.0
    nnan = 0
    for _, g, p in pairs:
        for v in (g, p):
            if math.isnan(v) or math.isinf(v):
                nnan += 1
            else:
                scale = max(scale, abs(v))
    s = {'n': len(pairs), 'scale': scale, 'nnan': nnan, 'nbad': 0,
         'nsign': 0, 'maxdiff': 0.0, 'worst': None, 'gscale': 0.0,
         'pscale': 0.0}
    for where, g, p in pairs:
        if not (math.isnan(g) or math.isinf(g)):
            s['gscale'] = max(s['gscale'], abs(g))
        if not (math.isnan(p) or math.isinf(p)):
            s['pscale'] = max(s['pscale'], abs(p))
        if math.isnan(g) or math.isnan(p) or math.isinf(g) or math.isinf(p):
            if not (math.isnan(g) and math.isnan(p)):
                s['nbad'] += 1
            continue
        d = abs(g - p)
        fl = max(scale, fscale)
        tol = RTOL * max(abs(g), abs(p)) + AFLOOR * fl + ABSFLOOR
        if d > tol:
            s['nbad'] += 1
            if abs(g) > SIGN_FLOOR * fl and abs(p) > SIGN_FLOOR * fl \
                    and (g > 0) != (p > 0):
                s['nsign'] += 1
        if d > s['maxdiff']:
            s['maxdiff'] = d
            s['worst'] = [where, g, p]
    s['rel'] = s['maxdiff'] / scale if scale > 0 else 0.0
    s['fscale'] = max(scale, fscale)
    s['frel'] = s['maxdiff'] / s['fscale'] if s['fscale'] > 0 else 0.0
    return s


def compare_labelled(rg, rp, lg, lp, rec, pairing, include_deriv):
    out = {}
    dg = rg.get(rec, {})
    dp = rp.get(rec, {})
    if not dg or not dp:
        return out
    posg = {l: i for i, l in enumerate(lg)}
    posp = {l: i for i, l in enumerate(lp)}
    keep = set(lg[i] for i in base_set(lg)) if not include_deriv else set(lg)
    common = [l for l in lg if l in posp and l in keep]
    allpairs = {}
    for l in common:
        ig, ip = posg[l], posp[l]
        pairs = []
        for ng, vg in dg.items():
            if ng is None:
                continue
            np_ = pairing.get(ng, ng) if pairing is not None else ng
            vp = dp.get(np_)
            if vp is None or ig >= len(vg) or ip >= len(vp):
                continue
            pairs.append((ng, fnum(vg[ig]), fnum(vp[ip])))
        if pairs:
            allpairs[l] = pairs
    fam = {}
    for l, pairs in allpairs.items():
        m = 0.0
        for _, g, p in pairs:
            for v in (g, p):
                if not (math.isnan(v) or math.isinf(v)):
                    m = max(m, abs(v))
        k = group_key(l)
        fam[k] = max(fam.get(k, 0.0), m)
    for l, pairs in allpairs.items():
        out[l] = compare_pairs(pairs, fam[group_key(l)])
    return out


def compare_raw(rg, rp, rec, pairing=None):
    dg = rg.get(rec, {})
    dp = rp.get(rec, {})
    out = {}
    for k in dg:
        kp = pairing.get(k, k) if pairing is not None else k
        if kp not in dp:
            continue
        vg, vp = dg[k], dp[kp]
        if len(vg) != len(vp) or not vg:
            continue
        if not all(isnum(x) for x in vg + vp):
            continue
        pairs = [((k, i), fnum(a), fnum(b)) for i, (a, b) in enumerate(zip(vg, vp))]
        for i, pr in enumerate(pairs):
            out.setdefault('%s[%d]' % (rec, i), []).append(pr)
    return {k: compare_pairs(v) for k, v in out.items()}


def read_his(path):
    rows = []
    for line in open(path, errors='replace'):
        t = line.split()
        if t and all(isnum(x) for x in t):
            rows.append([fnum(x) for x in t])
    return rows


def compare_his(dg, dpro):
    out = {}
    hist = lambda f: f.endswith('.his') or f == 'incremental_driver_result.txt'
    fg = {f for f in os.listdir(dg) if hist(f)}
    fp = {f for f in os.listdir(dpro) if hist(f)}
    for f in sorted(fg & fp):
        a, b = read_his(os.path.join(dg, f)), read_his(os.path.join(dpro, f))
        if not a or not b:
            continue
        n = min(len(a), len(b))
        ncol = min(len(a[0]), len(b[0]))
        names = {}
        if f == 'incremental_driver_result.txt':
            for line in open(os.path.join(dpro, f), errors='replace'):
                if line.startswith('time '):
                    names = dict(enumerate(line.split()))
                    break
        for c in range(ncol):
            pairs = [((f, r), a[r][c], b[r][c]) for r in range(n)
                     if c < len(a[r]) and c < len(b[r])]
            s = compare_pairs(pairs)
            s['rows'] = [len(a), len(b)]
            out['%s:%s' % (f, names.get(c, 'col%d' % c))] = s
    return out


# ===TARGETS===


def read_targets(dat_text):
    t = re.sub(r'\([^)]*\)', ' ', dat_text, flags=re.S)
    items, values = {}, {}
    for m in re.finditer(r'^\s*target_item\s+(\d+)\s+(.*)$', t, flags=re.M):
        items[int(m.group(1))] = m.group(2).split()
    for m in re.finditer(r'^\s*target_value\s+(\d+)\s+(.*)$', t, flags=re.M):
        v = m.group(2).split()
        try:
            values[int(m.group(1))] = (float(v[0]), float(v[1]) if len(v) > 1 else 0.0)
        except ValueError:
            pass
    return {k: (items[k], values[k]) for k in items if k in values}


OPLAB = {'total_pressure': 'to_pres', 'static_pressure': 'st_pres',
         'dynamic_pressure': 'dy_pres'}


def calcul_labels(dat_text):
    """GNU dumps carry no post_calcul_label: rebuild it from the input's
    post_calcul (dof, operation) pairs with the Professional label names"""
    t = re.sub(r'\([^)]*\)', ' ', dat_text, flags=re.S)
    m = re.search(r'^\s*post_calcul\s+((?:-\S+\s*)+)', t, flags=re.M)
    if not m:
        return []
    toks = [x.lstrip('-') for x in m.group(1).split()]
    return [OPLAB.get(op, op) for op in toks[1::2]]


def value(recs, item, dat_labels=()):
    lab = labels(recs)
    if len(item) < 3:
        return None
    kind, idx, what = item[0], item[1], item[2]
    try:
        idx = int(idx)
    except ValueError:
        return None
    rec = {'-node_dof': 'node_dof', '-post_point_dof': 'post_point_dof',
           '-node_dof_calcul': 'node_dof_calcul',
           '-post_point_dof_calcul': 'post_point_dof_calcul',
           '-node_rhside': 'node_rhside', '-node': 'node'}.get(kind)
    if rec is None:
        return None
    vals = recs.get(rec, {}).get(idx)
    if vals is None:
        return None
    if rec in ('node_dof', 'post_point_dof', 'node_rhside'):
        name = what.lstrip('-')
        if name not in lab:
            return None
        i = lab.index(name)
    else:
        try:
            i = int(what)
        except ValueError:
            # named post_calcul item: position in post_calcul_label
            cl = [t.lstrip('-') for t in recs.get('post_calcul_label', {}).get(None, [])] \
                or list(dat_labels)
            if what.lstrip('-') not in cl:
                return None
            i = cl.index(what.lstrip('-'))
    if i >= len(vals):
        return None
    return fnum(vals[i])


def evaluate(dat_path, dbs_path):
    txt = open(dat_path, errors='replace').read()
    tg = read_targets(txt)
    cl = calcul_labels(txt)
    try:
        recs = parse_dbs(dbs_path)
    except OSError:
        return []
    out = []
    for k in sorted(tg):
        item, (v, tol) = tg[k]
        x = value(recs, item, cl)
        if x is None:
            out.append((k, ' '.join(item), v, tol, None, 'NA'))
        else:
            out.append((k, ' '.join(item), v, tol, x, 'PASS' if abs(x - v) <= tol else 'FAIL'))
    return out


# ===AGG===


FAMILIES = [
    ('displacement/velocity', r'^(vel|veli|dis|acc)[xyz]$'),
    ('beam/rotation', r'^(rot|beam|mom|she|nor)'),
    ('total strain', r'^ept[xyz]{2}$'),
    ('elastic strain', r'^epe'),
    ('plastic strain', r'^epp'),
    ('stress', r'^sig[xyz]{2}$'),
    ('stress pressure history', r'^sph$'),
    ('plastic kappa/f', r'^(kap|kapshear|f|fnl|plf)$'),
    ('damage', r'^dam$'),
    ('history variables', r'^(hisv\d+|pc|isa|ig\w*|hsg\w*)$'),
    ('groundflow', r'^(pres|tpres|gvel[xyz]|gsat|gpres\w*)$'),
    ('temperature', r'^temp$'),
    ('wave', r'^f?scal$'),
    ('density', r'^dens$'),
    ('maxwell stress', r'^mxs'),
    ('void fraction', r'^void$'),
]

FEATURES = [
    ('vonmises', r'group_materi_plasti_vonmises'),
    ('mohr_coulomb', r'group_materi_plasti_mohr_coul'),
    ('druck_prag', r'group_materi_plasti_druck'),
    ('camclay', r'camclay'),
    ('cap', r'group_materi_plasti_cap'),
    ('hardsoil', r'hardsoil'),
    ('hypo', r'group_materi_(plasti_)?hypo|materi_history_variables'),
    ('diprisco', r'diprisco'),
    ('tension_cutoff', r'plasti_tension'),
    ('damage', r'group_materi_damage'),
    ('viscous', r'visco|maxwell'),
    ('groundflow', r'groundflow'),
    ('condif', r'condif'),
    ('interface', r'-interface|group_interface'),
    ('contact', r'contact'),
    ('truss', r'truss'),
    ('beam', r'-beam|group_beam|beam_rotation'),
    ('spring', r'spring'),
    ('memory_updated', r'-updated\b|-updated_jaumann'),
    ('memory_updated_linear', r'-updated_linear|-updated_without_rotation'),
    ('memory_total', r'-total\b|-total_linear|-total_piola|materi_displacement\b'),
    ('follow_material', r'follow_material'),
    ('inertia', r'inertia_apply\s+-yes|options_inertia\s+-yes'),
    ('refine/remesh', r'refine|remesh|new_mesh'),
    ('delete', r'delete'),
    ('generate', r'control_mesh_generate'),
    ('stress_force', r'materi_stress_force|post_calcul'),
    ('axisym', r'axisym'),
    ('3d', r'number_of_space_dimensions\s+3'),
    ('1d', r'number_of_space_dimensions\s+1'),
]


def family(label):
    for fam, rx in FAMILIES:
        if re.search(rx, label):
            return fam
    return 'other'


def features(text):
    t = re.sub(r'\([^)]*\)', ' ', text, flags=re.S)
    return [f for f, rx in FEATURES if re.search(rx, t)]


def target_vars(text):
    t = re.sub(r'\([^)]*\)', ' ', text, flags=re.S)
    out = []
    for m in re.finditer(r'^\s*target_item\s+\d+\s+(.*)$', t, flags=re.M):
        for w in m.group(1).split():
            if w.startswith('-') and not re.match(r'^-?\d', w) and \
                    w not in ('-node_dof', '-post_point_dof', '-post_node',
                              '-node_dof_calcul', '-post_node_result',
                              '-element_dof', '-post_point', '-node'):
                out.append(w.lstrip('-'))
    return out


FATAL = re.compile(r'^(Program error.*|Error(?!: the Bi-CG iteration broke down).*|.*[Aa]ssert.*)$', re.M)


def load_cls(path):
    d = {}
    if not os.path.isfile(path):
        return d
    base = os.path.dirname(path)
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        while len(f) < 6:
            f.append('')
        e = {'cls': f[1], 'rc': f[2], 'dt': f[3], 'msg': f[4], 'conv': f[5]}
        if e['cls'] == 'OK':
            log = os.path.join(base, f[0], 'run.log')
            if os.path.isfile(log):
                m = FATAL.search(open(log, errors='replace').read())
                if m:
                    # a dump was written but the run hit a fatal error:
                    # its state is partial, never compared
                    e['cls'] = 'ERRDUMP'
                    e['msg'] = m.group(1)[:90]
        d[f[0]] = e
    return d


def run_cmp(a, b, name):
    try:
        return compare_dirs(a, b, name)
    except Exception as exc:  # a malformed dump must not stop the audit
        return {'error': repr(exc)[-400:]}


def mism(res, sections=('node_dof', 'post_point_dof', 'node_dof_calcul', 'his')):
    out = []
    for sec in sections:
        for var, s in res.get(sec, {}).items():
            if s['nbad'] > 0 or s.get('nnan', 0) > 0:
                out.append({'sec': sec, 'var': var, 'nbad': s['nbad'], 'n': s['n'],
                            'rel': s['rel'], 'scale': s['scale'], 'nsign': s['nsign'],
                            'frel': s.get('frel', s['rel']), 'fscale': s.get('fscale', s['scale']),
                            'nnan': s.get('nnan', 0), 'worst': s['worst'],
                            'gscale': s.get('gscale'), 'pscale': s.get('pscale')})
    return out


# ===CLUSTERS===

TH = 1e-2  # a variable counts as mismatching when frel >= 1%

CLUSTERS = [
    # id, title, verdict, relevance (1 post-processing .. 3 geotechnical core)
    ('K01', 'GNU -updated_linear memory model', 'GNU-suspect', 3),
    ('K02', 'Iteration lag: non-principal dofs report the pre-solve state', 'GNU-suspect', 2),
    ('K03', 'Default convection_apply (GNU -yes, Pro manual -no)', 'GNU-suspect', 2),
    ('K04', 'Pro-authored analytic targets missed by the GNU (feature bugs)', 'GNU-suspect', 2),
    ('K05', 'Contact (1D follow_material contact; contact_target_geometry alias C5)', 'GNU-suspect', 2),
    ('K06', 'Elastic power laws: Pro adds the base young/poisson', 'GNU-suspect (semantic)', 3),
    ('K07', 'Positional node_dof initial state over different dof layouts', 'INPUT-artifact', 1),
    ('K08', 'bounda_time_increment load tables', 'PRO-suspect', 1),
    ('K09', 'Groundflow: -pres dof datum (Pro stores head h, GNU pore pressure); totals agree', 'BOTH-plausible / INPUT-artifact', 2),
    ('K10', 'Groundflow velocity gvel = 0 in the Professional', 'PRO-suspect', 1),
    ('K11', 'Nodal recovery of total strain (stress recovered with the Pro operator, strain not)', 'GNU-suspect', 2),
    ('K12', 'post_calcul item conventions (prival order, calcul items, section forces)', 'BOTH-plausible', 1),
    ('K13', 'Direct/visco plasticity: GNU stores no plastic strain / f', 'GNU-suspect', 2),
    ('K14', 'Hypoplasticity: history-variable slots and stress path', 'UNKNOWN', 3),
    ('K15', 'Other plasticity return-mapping differences (MC, von Mises, cap, hardsoil)', 'UNKNOWN (pre-existing MC)', 3),
    ('K16', 'GNU extensions / Pro ignores the record (GNU-only features)', 'INPUT-artifact', 1),
    ('K17', 'Truss / beam / support divergence (GNU blows up on Pro inputs)', 'GNU-suspect', 2),
    ('K18', 'Mesh refinement, higher-order and special elements', 'UNKNOWN', 2),
    ('K19', 'Nonlinear kinematics residual (updated / follow_material / large strain)', 'BOTH-plausible / UNKNOWN', 2),
    ('K20', 'Axisymmetric hoop (zz) nodal strain/stress', 'BOTH-plausible', 2),
    ('K21', 'Viscous fluid / Navier-Stokes (incnav, tube, examp1)', 'UNKNOWN', 1),
    ('K22', 'Phreatic level: wet/dry density in the initial gravity stress', 'GNU-suspect', 3),
    ('K23', 'Round-off in identically-zero fields (rigid motion); false positives', 'BOTH-plausible', 1),
    ('K24', '3D coupled consolidation (pore pressure and displacement 4-5%)', 'UNKNOWN', 3),
    ('K99', 'Unclassified', 'UNKNOWN', 1),
]
TITLE = {c[0]: c for c in CLUSTERS}

K04 = {'pro_other/condif7', 'pro_other/force18', 'pro_other/elasti7', 'pro_other/tension_direct1',
       'pro_other/mohr_coul_direct8', 'pro_other/stress_pressure_factor', 'pro_other/incnav3',
       'pro_other/axisym9', 'pro_other/distri2', 'pro_other/ground9', 'pro_other/tet10_2',
       'sfnet/incnav3', 'pro_other/mohr_coul_direct4'}
K05 = {'pro_other/contac1', 'pro_other/contac2', 'pro_other/contac4', 'pro_other/axisym3',
       'sfnet/contac2', 'sfnet/contac3', 'sfnet/examp12', 'sfnet/examp9', 'sfnet/examp17',
       'pro_validation/validation_10', 'pro_other/pressed_rubber_cylinder', 'vs/contact_ctrl_apply'}
K06 = {'vs/myoung6', 'vs/myoung6_e2', 'vs/myoung6_e3', 'vs/mpower'}
K12 = {'pro_other/calcul1', 'pro_other/calcul2', 'sfnet/calcul1', 'sfnet/calcul2', 'pro_other/ctvmis2',
       'sfnet/ctvmis2', 'pro_other/force10', 'pro_other/force13', 'vs/gforce10', 'vs/tpef_force10',
       'pro_validation/validation_5', 'sfnet/examp10', 'vs/msf_axisym', 'vs/msf_tunnel3d'}
K16 = {'sfnet/transv1', 'sfnet/temp4', 'sfnet/temp5', 'sfnet/slide1', 'vs/condif_rad', 'vs/strain_settle',
       'vs/strain_settle_diag', 'vs/mfactor', 'vs/mdirect_gate', 'vs/mvoid', 'vs/mvoid_low',
       'vs/vtk_other1', 'vs/mstrain_compression',
       'vs/mcrunch', 'vs/mcrunch_low', 'vs/cmat_gate'}
K17 = {'pro_other/truss2', 'pro_other/truss3', 'pro_other/truss5', 'pro_other/truss8', 'pro_other/truss13',
       'pro_other/truss14', 'pro_other/trubea1', 'pro_other/trubea2', 'pro_other/trubea3',
       'pro_other/support_edge_normal_2', 'pro_other/support_edge_normal_3', 'pro_other/update2',
       'pro_other/dynamic3'}


def fams(e, th=TH, secs=('node_dof',)):
    return {m['family'] for m in e.get('mism', []) if m['sec'] in secs and m['frel'] >= th}


def vars_(e, th=TH, secs=('node_dof',)):
    return {m['var'] for m in e.get('mism', []) if m['sec'] in secs and m['frel'] >= th}


def input_text(runs, e):
    p = os.path.join(runs, e['corpus'], e['gside'], e['name'], e['name'] + '.dat')
    try:
        return re.sub(r'\([^)]*\)', ' ', open(p, errors='replace').read(), flags=re.S)
    except OSError:
        return ''


def nonzero_node_dof(txt):
    for m in re.finditer(r'^\s*node_dof\s+(.*)$', txt, flags=re.M):
        t = m.group(1).split()
        for v in t[1:]:
            try:
                if float(v) != 0.0:
                    return True
            except ValueError:
                return True  # -ra / -from forms: positional block
    return False


def assign(runs, key, e, vanish):
    f = fams(e)
    v = vars_(e)
    feats = set(e['features'])
    txt = input_text(runs, e)
    conv = (e.get('pro') or {}).get('conv', '')
    meta = e['cmp_meta']
    if key in vanish:
        return 'K02'
    if key in ('sfnet/condif2', 'sfnet/condif3'):
        return 'K03'
    if key in K05 or ('contact_geometry' in conv and 'contact' in feats):
        return 'K05'
    if key in K04:
        return 'K04'
    if key in ('pro_other/ground18', 'pro_other/ground8'):
        return 'K22'
    if key == 'pro_other/large2':
        return 'K24'
    if key in K06:
        return 'K06'
    if key in K17:
        return 'K17'
    if key in K12 or (not f and fams(e, secs=('node_dof_calcul',))):
        return 'K12'
    if nonzero_node_dof(txt) and (meta['labels_only_gnu'] or meta['labels_only_pro']):
        return 'K07'
    if 'bounda_time_increment' in txt and 'displacement/velocity' in f:
        return 'K08'
    if f and f <= {'groundflow'} and v <= {'gvelx', 'gvely', 'gvelz'}:
        return 'K10'
    if 'groundflow' in f:
        return 'K09'
    if key in K16:
        return 'K16'
    if key in ('vs/dpline_move', 'vs/dpoint1', 'vs/dpoint_time'):
        return 'K23'
    if 'hypo' in feats:
        return 'K14'
    ul = bool(re.search(r'-updated_linear|-updated_without_rotation', txt)) or \
        '!updated_without_rotation' in conv
    orig_clean = 'orig_mism' in e and not vars_({'mism': e['orig_mism']})
    if ('!updated_without_rotation' in conv and orig_clean and f - {'total strain'}) or \
            key in ('pro_other/interac1', 'pro_other/matrix7', 'vs/mesh_act_grav'):
        return 'K01'
    if 'axisym' in feats and f and v <= {'eptzz', 'sigzz', 'eptxx', 'sigxx', 'velx', 'disx'}:
        return 'K20'
    if re.search(r'^\s*group_materi_viscosity\b', txt, flags=re.M) and \
            re.search(r'incnav|tube|examp1\b', key):
        return 'K21'
    if f and f <= {'total strain', 'elastic strain'}:
        return 'K11'
    if f & {'plastic strain', 'plastic kappa/f'} and \
            re.search(r'plasti_(tension_direct|mohr_coul_direct|compression_direct)|visco', txt) and \
            not ({'stress', 'displacement/velocity'} & f):
        return 'K13'
    if f & {'plastic strain', 'plastic kappa/f'} or \
            feats & {'mohr_coulomb', 'vonmises', 'cap', 'hardsoil', 'druck_prag', 'camclay'}:
        return 'K15'
    if feats & {'refine/remesh'} or re.search(r'-(quad8|tria6|tet10|hex20|prism15|quad9|hex27)\b', txt):
        return 'K18'
    if feats & {'memory_updated', 'follow_material'} or ul:
        return 'K19'
    return 'K99'



# ===========================================================================
# command entry points
# ===========================================================================

def compare_dirs(dg, dpro, name):
    """Final-state comparison of one test (both run directories)."""
    rg = parse_dbs(os.path.join(dg, name + '.dbs'))
    rp = parse_dbs(os.path.join(dpro, name + '.dbs'))
    lg, lp = labels(rg), labels(rp)
    cg, cp = coords(rg), coords(rp)
    pairing, method = pair_nodes(cg, cp)
    res = {
        'name': name,
        'labels_gnu': len(lg), 'labels_pro': len(lp),
        'labels_only_gnu': [l for l in lg if l not in lp],
        'labels_only_pro': [l for l in lp if l not in lg],
        'nodes_gnu': len(cg), 'nodes_pro': len(cp), 'node_pairing': method,
        'time_gnu': rg.get('time_current', {}).get(None),
        'time_pro': rp.get('time_current', {}).get(None),
        'node_dof': compare_labelled(rg, rp, lg, lp, 'node_dof', pairing, False),
        'post_point_dof': compare_labelled(rg, rp, lg, lp, 'post_point_dof', None, False),
        'node_dof_calcul': compare_raw(rg, rp, 'node_dof_calcul', pairing),
        'his': compare_his(dg, dpro),
    }
    alld = compare_labelled(rg, rp, lg, lp, 'node_dof', pairing, True)
    res['node_dof_deriv'] = {k: v for k, v in alld.items() if k not in res['node_dof']}
    nr = {}
    rgg, rpp = rg.get('node_rhside', {}), rp.get('node_rhside', {})
    if rgg and rpp:
        vg = next(iter(rgg.values()))
        vp = next(iter(rpp.values()))
        if len(vg) == len(lg) and len(vp) == len(lp):
            nr = compare_labelled(rg, rp, lg, lp, 'node_rhside', pairing, False)
            # only the primary (principal) unknowns carry a residual/reaction
            nr = {k: v for k, v in nr.items()
                  if re.match(r'^(vel[xyz]|temp|pres|f?scal|rot[xyz]|dens)$', k)}
    res['node_rhside'] = nr
    return res


def cmd_run(argv):
    ap = argparse.ArgumentParser(prog='audit_all_variables.py run')
    ap.add_argument('--bin', required=True)
    ap.add_argument('--src', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--timeout', type=float, default=45)
    ap.add_argument('--jobs', type=int, default=7)
    ap.add_argument('--convert', default='')
    ap.add_argument('--only', default='')
    ap.add_argument('--exclude', default='')
    ap.add_argument('--ld', default='')
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    all_files = sorted(os.listdir(args.src))
    names = sorted(f[:-4] for f in all_files if f.endswith('.dat'))
    if args.only:
        keep = set(args.only.split(','))
        names = [n for n in names if n in keep]
    if args.exclude:
        drop = set(args.exclude.split(','))
        names = [n for n in names if n not in drop]
    with ThreadPoolExecutor(args.jobs) as ex:
        res = list(ex.map(lambda n: run_one(args, n, all_files), names))
    path = os.path.join(args.out, 'classification.tsv')
    old = {}
    if args.only and os.path.isfile(path):
        for line in open(path):
            old[line.split('\t')[0]] = line
    for r in res:
        old[r[0]] = '\t'.join(str(x) for x in r) + '\n'
    with open(path, 'w') as fh:
        for k in sorted(old):
            fh.write(old[k])
    print(args.out, dict(collections.Counter(r[1] for r in res)))


def _cmp_job(j):
    key, kind, a, b, name = j
    return key, kind, run_cmp(a, b, name)


def cmd_aggregate(argv):
    runs = argv[0]
    summary = {}
    jobs = []
    for corpus in sorted(os.listdir(runs)):
        cdir = os.path.join(runs, corpus)
        if not os.path.isdir(os.path.join(cdir, 'gnu')):
            continue
        g = load_cls(os.path.join(cdir, 'gnu', 'classification.tsv'))
        p = load_cls(os.path.join(cdir, 'pro', 'classification.tsv'))
        gc = load_cls(os.path.join(cdir, 'gnuconv', 'classification.tsv'))
        for name in sorted(set(g) | set(p)):
            e = {'corpus': corpus, 'name': name, 'gnu': g.get(name), 'pro': p.get(name),
                 'gnuconv': gc.get(name)}
            src = os.path.join(cdir, 'gnu', name, name + '.dat')
            txt = open(src, errors='replace').read() if os.path.isfile(src) else ''
            e['features'] = features(txt)
            e['targets'] = target_vars(txt)
            summary[corpus + '/' + name] = e
            # GNU-dialect corpora: both binaries run the SAME converted input
            gside = 'gnuconv' if gc else 'gnu'
            e['gside'] = gside
            ge = e[gside]
            if ge and e['pro'] and ge['cls'] == 'OK' and e['pro']['cls'] == 'OK':
                jobs.append((corpus + '/' + name, 'cmp', os.path.join(cdir, gside, name),
                             os.path.join(cdir, 'pro', name), name))
            if e['gnu'] and e['gnuconv'] and e['gnu']['cls'] == 'OK' and \
                    e['gnuconv']['cls'] == 'OK':
                jobs.append((corpus + '/' + name, 'convcheck', os.path.join(cdir, 'gnu', name),
                             os.path.join(cdir, 'gnuconv', name), name))
            # GNU on its ORIGINAL input vs Pro on the converted one: tells
            # whether a GNU(conv) mismatch comes from the conversion itself
            if gc and e['gnu'] and e['pro'] and e['gnu']['cls'] == 'OK' and \
                    e['pro']['cls'] == 'OK':
                jobs.append((corpus + '/' + name, 'origcmp', os.path.join(cdir, 'gnu', name),
                             os.path.join(cdir, 'pro', name), name))
    with ProcessPoolExecutor(min(14, os.cpu_count() or 4)) as ex:
        for key, kind, res in ex.map(_cmp_job, jobs, chunksize=4):
            e = summary[key]
            if 'error' in res:
                e[kind + '_error'] = res['error']
                continue
            if kind == 'cmp':
                e['cmp_meta'] = {k: res[k] for k in ('labels_only_gnu', 'labels_only_pro',
                                                     'nodes_gnu', 'nodes_pro', 'node_pairing',
                                                     'time_gnu', 'time_pro')}
                e['mism'] = mism(res)
                e['mism_deriv'] = mism(res, ('node_dof_deriv',))
                e['mism_rhside'] = mism(res, ('node_rhside',))
                e['nvars'] = sum(len(res[s]) for s in ('node_dof', 'post_point_dof',
                                                        'node_dof_calcul', 'his'))
                e['nhis'] = len(res['his'])
            elif kind == 'origcmp':
                e['orig_mism'] = mism(res)
            else:
                e['conv_mism'] = mism(res)
    for e in summary.values():
        for m in e.get('mism', []):
            v = m['var']
            m['family'] = 'history-file' if m['sec'] == 'his' else (
                'post_calcul' if m['sec'] == 'node_dof_calcul' else family(v))
            m['targeted'] = v in e['targets']
    # the test's own targets, evaluated identically on both dumps
    for key, e in summary.items():
        if 'mism' not in e:
            continue
        cdir = os.path.join(runs, e['corpus'])
        dat = os.path.join(cdir, 'gnu', e['name'], e['name'] + '.dat')
        if not os.path.isfile(dat):
            continue
        e['tgt_gnu'] = evaluate(dat, os.path.join(cdir, e['gside'], e['name'], e['name'] + '.dbs'))
        e['tgt_pro'] = evaluate(dat, os.path.join(cdir, 'pro', e['name'], e['name'] + '.dbs'))
    json.dump(summary, open(os.path.join(runs, 'summary.json'), 'w'), indent=1)
    print('tests', len(summary), 'compared', sum(1 for e in summary.values() if 'mism' in e))


def cmd_iter_effect(argv):
    """which node_dof mismatches vanish when both binaries are forced to 5
    equilibrium iterations (RUNS_ITER built with --convert iter|pro+iter)"""
    a_dir, b_dir, out = argv[:3]
    A = json.load(open(os.path.join(a_dir, 'summary.json')))
    B = json.load(open(os.path.join(b_dir, 'summary.json')))

    def sig(e):
        return {m['var'] for m in e.get('mism', []) if m['sec'] == 'node_dof' and m['frel'] >= 1e-2}
    rows = []
    for k, e in A.items():
        if k not in B or 'mism' not in B[k] or not sig(e):
            continue
        a, b = sig(e), sig(B[k])
        tag = 'VANISH' if not b else ('PARTIAL' if b < a else 'PERSIST')
        rows.append((tag, k, sorted(a), sorted(b)))
    json.dump(rows, open(out, 'w'))
    print(dict(collections.Counter(r[0] for r in rows)))


def cmd_clusters(argv):
    runs = argv[0]
    S = json.load(open(os.path.join(runs, 'summary.json')))
    vanish = set()
    if len(argv) > 1 and os.path.isfile(argv[1]):
        for tag, k, a, b in json.load(open(argv[1])):
            if tag == 'VANISH':
                vanish.add(k)
    by = collections.defaultdict(list)
    for key, e in sorted(S.items()):
        if 'mism' not in e or not fams(e) and not fams(e, secs=('node_dof_calcul',)):
            continue
        by[assign(runs, key, e, vanish)].append(key)
    json.dump(dict(by), open(os.path.join(runs, 'clusters.json'), 'w'), indent=1)
    rows = []
    for cid, title, verdict, rel in CLUSTERS:
        ks = by.get(cid, [])
        if not ks:
            continue
        mags = sorted(min(max([x['frel'] for x in S[k]['mism']
                               if x['sec'] in ('node_dof', 'node_dof_calcul')] or [0]), 1.0)
                      for k in ks)
        med = mags[len(mags) // 2]
        rows.append((len(ks) * med * rel, cid, title, verdict, rel, len(ks), med, ks))
    rows.sort(key=lambda r: -r[0])
    for r in rows:
        print('%s score=%5.1f n=%3d med=%.2f rel=%d %-28s %s' % (r[1], r[0], r[5], r[6], r[4], r[3], r[2]))
    json.dump(rows, open(os.path.join(runs, 'cluster_table.json'), 'w'), indent=1)


def cmd_stats(argv):
    """markdown tables of the report"""
    runs = argv[0]

    S = json.load(open(os.path.join(runs, 'summary.json')))
    CT = json.load(open(os.path.join(runs, 'cluster_table.json')))
    CORP = collections.OrderedDict([
        ('sfnet', 'sfnet 2014 (GNU dialect, converted)'),
        ('vs', 'validation-suite test-2014 (reconstructed, converted)'),
        ('t26', 'validation-suite test-2026 (own tests, converted)'),
        ('pro_other', 'Pro distribution test/other (Pro dialect, as is)'),
        ('pro_validation', 'Pro distribution test/validation'),
        ('pro_small', 'Pro distribution tutorial / easy_tochnog / incremental_driver'),
    ])


    def corp(e):
        c = e['corpus']
        return c if c in CORP else 'pro_small'


    def cls(x):
        return x['cls'] if x else 'NA'


    cov = collections.defaultdict(collections.Counter)
    why = collections.defaultdict(collections.Counter)
    for k, e in S.items():
        c = corp(e)
        g, p = cls(e[e['gside']]), cls(e['pro'])
        cov[c]['total'] += 1
        if g == 'OK' and p == 'OK':
            cov[c]['both'] += 1
        elif g == 'OK':
            cov[c]['gnu_only'] += 1
            why['Pro fails'][p] += 1
        elif p == 'OK':
            cov[c]['pro_only'] += 1
            why['GNU fails'][g] += 1
        else:
            cov[c]['neither'] += 1

    print('## coverage')
    print('| corpus | tests | run on both | GNU only | Pro only | neither |')
    print('|---|---|---|---|---|---|')
    tot = collections.Counter()
    for c, name in CORP.items():
        x = cov[c]
        tot.update(x)
        print('| %s | %d | %d | %d | %d | %d |' % (name, x['total'], x['both'], x['gnu_only'], x['pro_only'], x['neither']))
    print('| **all** | %d | %d | %d | %d | %d |' % (tot['total'], tot['both'], tot['gnu_only'], tot['pro_only'], tot['neither']))
    print()
    for side, cnt in why.items():
        print('%s (when the other binary ran): %s' % (side, dict(cnt)))
    print()


    def maxfrel(e, secs=('node_dof', 'post_point_dof', 'node_dof_calcul', 'his')):
        return max([m['frel'] for m in e['mism'] if m['sec'] in secs] or [0.0])


    res = collections.defaultdict(collections.Counter)
    for k, e in S.items():
        if 'mism' not in e:
            continue
        c = corp(e)
        r = res[c]
        r['compared'] += 1
        mf = maxfrel(e)
        nd = maxfrel(e, ('node_dof',))
        if mf == 0 and not e['mism']:
            r['identical'] += 1
        elif mf < 1e-2:
            r['minor'] += 1
        else:
            r['major'] += 1
            big = [m for m in e['mism'] if m['frel'] >= 1e-2]
            if not any(m['targeted'] for m in big):
                r['major_untargeted_only'] += 1
            if nd < 1e-2:
                r['major_not_in_node_dof'] += 1
        if e.get('nhis'):
            r['with_history'] += 1
    print('## match summary')
    print('| corpus | compared | match (within tol.) | minor only (<1%) | mismatch >=1% | of which only untargeted vars | of which node_dof ok (post_point/calcul/history only) | history files compared |')
    print('|---|---|---|---|---|---|---|---|')
    tot = collections.Counter()
    for c, name in CORP.items():
        r = res[c]
        tot.update(r)
        print('| %s | %d | %d | %d | %d | %d | %d | %d |' % (name, r['compared'], r['identical'], r['minor'], r['major'],
                                                             r['major_untargeted_only'], r['major_not_in_node_dof'], r['with_history']))
    print('| **all** | %d | %d | %d | %d | %d | %d | %d |' % (tot['compared'], tot['identical'], tot['minor'], tot['major'],
                                                              tot['major_untargeted_only'], tot['major_not_in_node_dof'], tot['with_history']))
    print()

    tg = collections.Counter()
    for k, e in S.items():
        if not e.get('tgt_gnu'):
            continue
        def st(rs):
            x = [r[5] for r in rs if r[5] != 'NA']
            return 'NA' if not x else ('PASS' if all(y == 'PASS' for y in x) else 'FAIL')
        src = 'Pro-authored targets' if e['corpus'].startswith('pro') else 'GNU-authored targets'
        tg[(src, st(e['tgt_gnu']), st(e['tgt_pro']))] += 1
    print('## targets (test\'s own target_value, evaluated on both dumps)')
    print('| targets | GNU | Pro | tests |')
    print('|---|---|---|---|')
    for (src, g, p), n in sorted(tg.items()):
        print('| %s | %s | %s | %d |' % (src, g, p, n))
    print()

    fam = collections.defaultdict(set)
    for k, e in S.items():
        for m in e.get('mism', []):
            if m['sec'] == 'node_dof' and m['frel'] >= 1e-2:
                fam[m['family']].add(k)
    print('## families (tests with a node_dof variable off by >=1% of its family scale)')
    print('| family | tests |')
    print('|---|---|')
    for f, v in sorted(fam.items(), key=lambda x: -len(x[1])):
        print('| %s | %d |' % (f, len(v)))
    print()
    print('## clusters')
    print('| rank | id | cluster | verdict | tests | median max rel. diff | relevance | score |')
    print('|---|---|---|---|---|---|---|---|')
    for i, r in enumerate(CT, 1):
        score, cid, title, verdict, rel, n, med, ks = r
        print('| %d | %s | %s | %s | %d | %.2f | %d | %.1f |' % (i, cid, title, verdict, n, med, rel, score))
    print()
    for r in CT:
        print('- %s: %s' % (r[1], ', '.join(r[7])))



def cmd_mism_list(argv):
    """comma-separated names of the mismatching tests of one corpus"""
    S = json.load(open(os.path.join(argv[0], 'summary.json')))
    print(','.join(sorted(e['name'] for e in S.values()
                          if e['corpus'] == argv[1] and e.get('mism'))))


COMMANDS = {'run': cmd_run, 'mism-list': cmd_mism_list, 'aggregate': cmd_aggregate, 'iter-effect': cmd_iter_effect,
            'clusters': cmd_clusters, 'stats': cmd_stats}


def main():
    if len(sys.argv) < 2 or sys.argv[1] in ('-h', '--help'):
        print(__doc__)
        return 0
    cmd, argv = sys.argv[1], sys.argv[2:]
    if cmd == 'compare':
        json.dump(compare_dirs(*argv[:3]), sys.stdout)
        return 0
    if cmd not in COMMANDS:
        print('unknown subcommand: %s' % cmd, file=sys.stderr)
        return 2
    COMMANDS[cmd](argv)
    return 0


if __name__ == '__main__':
    sys.exit(main())
