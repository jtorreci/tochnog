#!/usr/bin/env python3
"""Audit fixed-size scratch buffers passed to db() / get_group_data() /
array_*() calls (BACKLOG C49, F4 T13).

The database getter db(idat, index, ival, dval, length, version, task)
writes `length` values into ival (INTEGER records) or dval (DOUBLE records)
on GET / GET_IF_EXISTS, where `length` is the stored length of the record,
at most data_length[idat] (database.cc). A caller that passes a buffer
shorter than data_length[idat] overflows it as soon as a record with more
values is stored (input with more values, or a fixed-length record).

What the script does, mechanically:
  1. Reads #define constants from tochnog.h and every
     `type[X] = ...;` / `data_length[X] = expr;` / `fixed_length[X] = ...;`
     from database.cc. Run-time lengths (ndim, nuknwn, ...) are evaluated at
     their upper bound (ndim = MDIM, nuknwn = MUKNWN, ...).
  2. For every db(...) call with a literal record name and a reading task
     (GET, GET_IF_EXISTS, GET_AND_CHECK), it takes the buffer argument the
     record type selects (ival or dval), resolves its size from the nearest
     preceding local declaration in the same function (`name[N]` -> N,
     `&scalar` -> 1, `&a[k]` / `a+k` -> N-k) and reports the call when
     buffer size < data_length.
  3. The same for get_group_data(idat, gr, element, unknowns, values, n,
     task) (values receives up to data_length doubles).
  4. For array_*(...) helpers whose length argument is a constant N, every
     local array argument shorter than N is reported.

Usage:
  scripts/audit_db_buffers.py [--repo DIR] [--all] [--tsv FILE]
    default: print the hits (buffer smaller than the record) and a summary
    --all : also list the REVIEW calls (record chosen at run time with a
            buffer below DATA_ITEM_SIZE: usually the dummy of the other
            type behind a db_type() check) and the calls whose buffer size
            could not be resolved (pointer, parameter, variable offset)
    --tsv : write every scanned call (file, line, function, record, task,
            buffer, buffer size, data_length, fixed, verdict) to FILE
Exit status 0 when there is no hit, 1 otherwise (usable as a CI check).
"""
import argparse
import os
import re
import sys

RUNTIME_BOUNDS = {
    'ndim': 'MDIM', 'nuknwn': 'MUKNWN', 'npuknwn': 'MPUKNWN', 'nprinc': 'MPRINC',
    'npointmax': 'MPOINT', 'materi_maxwell_stress': '(MMAXWELL-1)',  # input.cc asserts < MMAXWELL
    'nder': '3',
}
READ_TASKS = {'GET', 'GET_IF_EXISTS', 'GET_AND_CHECK'}
# Output (or read) array arguments of other helpers with a fixed size
# contract: function -> {argument index: minimum size}. geometry() always
# clears normal[MDIM] and writes projection[ndim]; co is read (ndim) when
# inod < 0.
CALLEE_CONTRACT = {
    'geometry': {1: 'MDIM', 5: 'MDIM', 7: 'MDIM'},
    'msf_element_rule': {5: '3'},  # calcul_force.cc: writes nper[0..2]
}
TYPE_KW = r'(?:static\s+|const\s+|register\s+)*(?:long\s+int|unsigned\s+long|long|int|double|float|char)'


def strip_comments(text):
    """Blank comments and string/char literals, keeping offsets and newlines."""
    out = []
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if text.startswith('//', i):
            j = text.find('\n', i)
            j = n if j < 0 else j
            out.append(' ' * (j - i)); i = j
        elif text.startswith('/*', i):
            j = text.find('*/', i + 2)
            j = n if j < 0 else j + 2
            out.append(re.sub(r'[^\n]', ' ', text[i:j])); i = j
        elif c in '"\'':
            j = i + 1
            while j < n and text[j] != c:
                j += 2 if text[j] == '\\' else 1
            j = min(j + 1, n)
            out.append(c + ' ' * (j - i - 2) + c if j - i >= 2 else text[i:j]); i = j
        else:
            out.append(c); i += 1
    return ''.join(out)


class Const:
    def __init__(self, repo):
        self.defs = {}
        for f in ('tochnog.h',):
            for m in re.finditer(r'^\s*#define\s+(\w+)\s+([^/\n]+)', open(os.path.join(repo, f), errors='replace').read(), re.M):
                self.defs[m.group(1)] = m.group(2).strip()
        for k, v in RUNTIME_BOUNDS.items():
            self.defs.setdefault(k, v)
        self.base = dict(self.defs)
        self.globals = []
        self.globals = parse_decls(file_scope(strip_comments(open(os.path.join(repo, 'tochnog.h'), errors='replace').read())), 0, DECL_START_TOP)

    def file_defines(self, text):
        """Add the #define constants of one source file (reset per file)."""
        self.defs = dict(self.base)
        for m in re.finditer(r'^\s*#define\s+(\w+)\s+([^/\n]+)', text, re.M):
            self.defs[m.group(1)] = m.group(2).strip()

    def eval(self, expr, depth=0):
        expr = expr.strip()
        if not expr or depth > 10:
            return None
        def rep(m):
            w = m.group(0)
            if w in self.defs:
                v = self.eval(self.defs[w], depth + 1)
                return '(%d)' % v if v is not None else 'None'
            return 'None'
        e = re.sub(r'[A-Za-z_]\w*', rep, expr)
        if 'None' in e:
            return None
        try:
            return int(eval(e, {'__builtins__': {}}))
        except Exception:
            return None


def read_database(repo, const):
    text = strip_comments(open(os.path.join(repo, 'database.cc'), errors='replace').read())
    dtype, dlen, fixed = {}, {}, {}
    for m in re.finditer(r'\btype\[(\w+)\]\s*=\s*(\w+)\s*;', text):
        dtype[m.group(1)] = m.group(2)
    for m in re.finditer(r'\bdata_length\[(\w+)\]\s*=\s*([^;]+);', text):
        if m.group(1) == 'idat' or m.group(1) == 'data_number':
            continue
        v = const.eval(m.group(2))
        dlen[m.group(1)] = (max(v, 1) if v is not None else None, ' '.join(m.group(2).split()))
    for m in re.finditer(r'\bfixed_length\[(\w+)\]\s*=\s*(\d+)\s*;', text):
        fixed[m.group(1)] = int(m.group(2))
    return dtype, dlen, fixed


def split_args(text, start):
    """text[start] == '(' -> (list of args, end offset after ')')."""
    depth, i, args, cur = 0, start, [], start + 1
    while i < len(text):
        c = text[i]
        if c in '([{':
            depth += 1
        elif c in ')]}':
            depth -= 1
            if depth == 0:
                args.append(text[cur:i]); return [a.strip() for a in args], i + 1
        elif c == ',' and depth == 1:
            args.append(text[cur:i]); cur = i + 1
        i += 1
    return None, len(text)


def functions(text):
    """Yield (name, body_start, body_end) for top-level brace blocks."""
    depth, i, start = 0, 0, None
    for m in re.finditer(r'[{}]', text):
        if m.group(0) == '{':
            if depth == 0:
                start = m.start()
            depth += 1
        else:
            depth -= 1
            if depth == 0 and start is not None:
                head = text[max(0, start - 1500):start]
                nm = re.findall(r'(\w+)\s*\([^;{}]*\)\s*(?:const\s*)?$', head, re.S)
                yield (nm[-1] if nm else '?'), start, m.end()
                start = None


def buffer_size(arg, text, fstart, pos, const, sig):
    """Size in elements of the buffer expression arg, or (None, why)."""
    a = arg.strip()
    while a.startswith('(') and a.endswith(')'):
        a = a[1:-1].strip()
    off = 0
    m = re.fullmatch(r'&\s*(\w+)\s*\[\s*([^\]]+)\]', a)
    if m:
        name = m.group(1); off = const.eval(m.group(2))
        if off is None:
            return None, 'variable offset'
    else:
        m = re.fullmatch(r'(\w+)\s*\+\s*(.+)', a)
        if m and const.eval(m.group(2)) is not None:
            name, off = m.group(1), const.eval(m.group(2))
        else:
            if re.fullmatch(r'&\s*\w+\s*(\.|->)\s*\w+', a):
                return 1, 'scalar'
            m = re.fullmatch(r'&\s*(\w+)', a)
            if m:
                name = m.group(1)
                n, kind = decl_size(name, text, fstart, pos, const, sig)
                if kind == 'scalar':
                    return 1, 'scalar'
                if kind == 'array':
                    return n, 'array'
                return None, kind
            if re.fullmatch(r'\w+', a):
                name = a
            else:
                return None, 'expression'
    n, kind = decl_size(name, text, fstart, pos, const, sig)
    if kind != 'array':
        return None, kind
    if n is None:
        return None, 'array of unresolved size'
    return n - off, 'array'


DECL_START = re.compile(r'(?:^|(?<=[;{}(]))\s*(?:extern\s+)?' + TYPE_KW + r'\b', re.S)
DECL_START_TOP = re.compile(r'(?:^|(?<=[;{}]))\s*(?:extern\s+)?' + TYPE_KW + r'\b', re.S)


def file_scope(text):
    """text with every function body and every prototype parameter list
    blanked: what remains are globals and struct/class members."""
    t = list(text)
    for name, fs, fe in functions(text):
        head = text[max(0, fs - 1500):fs].rstrip()
        if head.endswith(')') or head.endswith('const'):
            t[fs:fe] = re.sub(r'[^\n]', ' ', text[fs:fe])
    t = ''.join(t)
    return re.sub(r'\([^;{}]*\)', lambda m: re.sub(r'[^\n]', ' ', m.group(0)), t)


def parse_decls(seg, base=0, start=None):
    """Declarations in seg -> list of (offset, name, kind, dims) where kind is
    'array', 'scalar' or 'pointer' and dims the list of bracket expressions."""
    out = []
    for m in (start or DECL_START).finditer(seg):
        i, depth, cur, parts = m.end(), 0, m.end(), []
        # the type keyword must be followed by a declarator, not by '('
        while i < len(seg):
            c = seg[i]
            if c in '([{':
                depth += 1
            elif c in ')]}':
                if depth == 0:
                    break
                depth -= 1
            elif c == ',' and depth == 0:
                parts.append((cur, seg[cur:i])); cur = i + 1
            elif c == ';' and depth == 0:
                break
            i += 1
        parts.append((cur, seg[cur:i]))
        for off, d in parts:
            d = d.split('=', 1)[0]
            dm = re.match(r'\s*(\*|&)?\s*([A-Za-z_]\w*)\s*((?:\[[^\]]*\]\s*)*)\s*$', d)
            if not dm:
                continue
            dims = re.findall(r'\[([^\]]*)\]', dm.group(3))
            kind = 'pointer' if dm.group(1) == '*' else ('scalar' if dm.group(1) == '&' else ('array' if dims else 'scalar'))
            out.append((base + off, dm.group(2), kind, dims))
    return out


_cache = {}


def lookup(decls, name, const):
    n = None
    for _, nm, kind, dims in decls:
        if nm == name:
            n = (kind, dims)
    if n is None:
        return None
    kind, dims = n
    if kind == 'pointer':
        return None, 'pointer'
    if kind == 'scalar':
        return 1, 'scalar'
    if any(x.strip() == '' for x in dims):
        return None, 'parameter'
    size = 1
    for x in dims:
        v = const.eval(x)
        if v is None:
            return None, 'array'
        size *= v
    return size, 'array'


def decl_size(name, text, fstart, pos, const, sig):
    """Size of the nearest declaration of name before pos: function body,
    then the signature, then file scope, then tochnog.h globals."""
    key = (id(text), fstart)
    if key not in _cache:
        _cache[key] = parse_decls(text[fstart:], fstart)
    body = [d for d in _cache[key] if d[0] < pos]
    fkey = (id(text), 'file')
    if fkey not in _cache:
        _cache[fkey] = parse_decls(file_scope(text), 0, DECL_START_TOP)
    params = []
    for part in (sig.split('(', 1)[1].rsplit(')', 1)[0] if '(' in sig else '').split(','):
        params += parse_decls(part.strip())
    for decls in (body, params,
                  [d for d in _cache[fkey] if d[0] < fstart],
                  const.globals):
        r = lookup(decls, name, const)
        if r is not None:
            return r
    return None, 'not found (global/macro)'


def scan(repo):
    const = Const(repo)
    dtype, dlen, fixed = read_database(repo, const)
    rows = []
    files = sorted(f for f in os.listdir(repo) if f.endswith(('.cc', '.c')))
    for f in files:
        raw = open(os.path.join(repo, f), errors='replace').read()
        text = strip_comments(raw)
        const.file_defines(text)
        _cache.clear()
        lines = [0]
        for m in re.finditer('\n', text):
            lines.append(m.end())
        import bisect
        lineno = lambda p: bisect.bisect_right(lines, p)
        for fname, fs, fe in functions(text):
            sig = text[max(0, fs - 1500):fs]
            sig = sig[sig.rfind(fname + '('):] if fname + '(' in sig else ''
            for m in re.finditer(r'\b(db|get_group_data|array_\w+|' + '|'.join(CALLEE_CONTRACT) + r')\s*\(', text[fs:fe]):
                call = m.group(1)
                p = fs + m.start()
                args, _ = split_args(text, fs + m.end() - 1)
                if args is None:
                    continue
                if call == 'db':
                    if len(args) != 7:
                        continue
                    rec = args[0].lstrip('-').strip()
                    task = args[6]
                    if task not in READ_TASKS:
                        rows.append((f, lineno(p), fname, call, args[0], task, '', '', '', '', 'skip:write'))
                        continue
                    if not re.fullmatch(r'[A-Z][A-Z0-9_]*', rec):
                        # record chosen at run time: the type is unknown, so
                        # both buffers may be written; any record can hold
                        # up to DATA_ITEM_SIZE values
                        sz = [buffer_size(x, text, fs, p, const, sig) for x in (args[2], args[3])]
                        small = [(x, n) for x, (n, k) in zip((args[2], args[3]), sz) if n is not None and n < const.eval('DATA_ITEM_SIZE')]
                        if small:
                            rows.append((f, lineno(p), fname, call, args[0], task, ' '.join(x for x, _ in small), min(n for _, n in small), const.eval('DATA_ITEM_SIZE'), 'dynamic record', 'REVIEW'))
                        else:
                            rows.append((f, lineno(p), fname, call, args[0], task, '', '', '', '', 'skip:dynamic'))
                        continue
                    t = dtype.get(rec, 'DOUBLE_PRECISION' if rec not in dtype else dtype[rec])
                    buf = args[2] if t == 'INTEGER' else args[3]
                    need = dlen.get(rec, (1, '1'))
                    if rec not in dlen:
                        need = (1, '1 (default)')
                    if task == 'GET_AND_CHECK':
                        v = const.eval(args[4])
                        need = (v, args[4]) if v is not None else need
                elif call == 'get_group_data':
                    if len(args) != 7:
                        continue
                    rec = args[0].lstrip('-').strip(); task = args[6]
                    if not re.fullmatch(r'[A-Z][A-Z0-9_]*', rec):
                        rows.append((f, lineno(p), fname, call, args[0], task, '', '', '', '', 'skip:dynamic'))
                        continue
                    buf = args[4]
                    need = dlen.get(rec, (1, '1 (default)'))
                elif call in CALLEE_CONTRACT:
                    if text[:fs + m.start()].rstrip().endswith(('void', 'long int', 'double')):
                        continue  # a declaration, not a call
                    for i, nexpr in CALLEE_CONTRACT[call].items():
                        if i >= len(args):
                            continue
                        size, kind = buffer_size(args[i], text, fs, p, const, sig)
                        nval = const.eval(nexpr)
                        if size is None:
                            rows.append((f, lineno(p), fname, call, '', 'arg%d' % i, args[i], '', nval, nexpr, 'unresolved:' + kind))
                            continue
                        rows.append((f, lineno(p), fname, call, '', 'arg%d' % i, args[i], size, nval, nexpr, 'HIT' if size < nval else 'ok'))
                    continue
                else:
                    if call in ('array_outproduct_2D', 'array_outproduct_3D'):
                        continue
                    if call == 'array_member':
                        cand, nexpr = [args[0]], args[2] if len(args) > 2 else ''
                    else:
                        if not args:
                            continue
                        nexpr = args[-1]
                        cand = [a for a in args[:-1]]
                    nval = const.eval(nexpr)
                    if nval is None:
                        continue
                    for a in cand:
                        if not re.fullmatch(r'&?\s*\w+(\s*\[[^\]]*\])?|\w+\s*\+\s*\w+', a.strip()) or re.fullmatch(r'[-+.\d]+|[A-Z_0-9]+', a.strip()):
                            continue
                        size, kind = buffer_size(a, text, fs, p, const, sig)
                        if size is None:
                            continue
                        verdict = 'HIT' if size < nval else 'ok'
                        rows.append((f, lineno(p), fname, call, '', '', a, size, nval, nexpr, verdict))
                    continue
                size, kind = buffer_size(buf, text, fs, p, const, sig)
                fx = fixed.get(rec, 1)
                if need[0] is None:
                    verdict = 'unresolved-length'
                elif size is None:
                    verdict = 'unresolved:' + kind
                else:
                    verdict = 'HIT' if size < need[0] else 'ok'
                rows.append((f, lineno(p), fname, call, rec, task, buf, size, need[0], need[1] + ('' if fx else ' (var)'), verdict))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--repo', default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--tsv')
    a = ap.parse_args()
    rows = scan(a.repo)
    if a.tsv:
        with open(a.tsv, 'w') as fh:
            fh.write('file\tline\tfunction\tcall\trecord\ttask\tbuffer\tbuffer_size\tneeded\tneeded_expr\tverdict\n')
            for r in rows:
                fh.write('\t'.join(str(x) for x in r) + '\n')
    hits = [r for r in rows if r[-1] == 'HIT']
    for r in hits:
        print('%s:%d %s(): %s %s %s -> buffer %s[%s] < %s (%s)' % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9]))
    for r in rows if a.all else []:
        if r[-1] == 'REVIEW':
            print('%s:%d %s(): REVIEW %s %s %s -> buffer(s) %s, min size %s (run-time record)' % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]))
    if a.all:
        for r in rows:
            if r[-1].startswith('unresolved'):
                print('%s:%d %s(): %s %s %s buffer %s: %s' % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[-1]))
    from collections import Counter
    c = Counter(r[-1].split(':')[0] for r in rows)
    calls = Counter(r[3] if not r[3].startswith('array_') else 'array_*' for r in rows)
    calls['geometry'] = len({(r[0], r[1]) for r in rows if r[3] == 'geometry'})
    print('# scanned: ' + ', '.join('%s %d' % kv for kv in sorted(calls.items())) +
          ' | verdicts: ' + ', '.join('%s %d' % kv for kv in sorted(c.items())), file=sys.stderr)
    return 1 if hits else 0


if __name__ == '__main__':
    sys.exit(main())
