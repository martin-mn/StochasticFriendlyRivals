#!/usr/bin/env python3
"""mk_strata_table.py -- the SI table of the thirteen boundary strata of the partner reduction (SI 'Partners').

    python3 mk_strata_table.py [--style box|table] [--label 'SI Box 3'] [--p1 DIR] [--out FILE] [--legend FILE]

Reads the certificate output P1/prove_strata.json (and P1/prove_strata.log for a cross-check) of
check/ref3/eq3_boundary and writes the box (default: SI Box 3, in the style of SI Box 2, so that no SI table is renumbered) or an \\sitable.  Every entry is
exact (sympy over Q).  Nothing is typed by hand: each formula in the table is computed here from the JSON, and
the script asserts every statement that the legend makes.

The title and the legend are part of the paper: they are read from PRFull/Text/si_box3_legend.tex (found from this
script's directory; --legend gives another file), which the code-and-data release 1.0.0 leaves out (it is added with
the paper in release 1.1.0).  The script always writes the box without its title and legend, si_box3_strata_body.tex
(<out>_body.tex, beside --out), and, when the legend file is found, the whole box, si_box3_strata.tex (--out).  The
assertions run in both cases.

What the table shows, per stratum (1,a,b,c) of {p_CC = 1, p_DD > 0} off the segments (a,b) = (0,1), (1,0):
  L      the eps^0 coefficient of the normaliser (the sum of the tree polynomials) of the self-play chain,
         written as a product of factors each positive on the open stratum.  It is recomputed here with the
         prover's own routine (prove_strata.limit_W) and compared with the line of prove_strata.log.
  gains  the deviators whose gain differs from the open-facet display of the SI with the stratum's values
         substituted (or for which that display is undefined), with their gain; all other deviators gain
         exactly what the display gives.
  Bernstein factors   the irreducible factors (of L, of the deviators' normalisers, of the multipliers) that are
         not a free coordinate x, 1-x or a constant, each written with the sign that makes it positive; each
         was certified by the coefficients of one strict sign after x = t/(1+t).

Asserted (the run stops otherwise): the JSON reports overall = True and, per stratum, self-play limit (1,0,0,0),
limit valid, sufficiency and necessity; mu = 0 for all 13 x 16 rows; every listed coefficient non-negative and
every denominator positive at the centre of the stratum (signs certified in the JSON); ALLD and 0111 gain
I/(1-a+c) and II/(1-a+2b) on every stratum; where a = 1, W_CD = W_DD = 0 for all 16 deviators; the self-play
eps-order is 0 and L agrees with the log; the two excluded segments are zeros of L on the closure of (1,a,b,1).
"""
import argparse
import json
import os
import re
import sys

import sympy as sp

sys.dont_write_bytecode = True     # import the prover without writing __pycache__ into the certificate directory

# the prover's directory: ../P1 when this script sits in check/ref3/eq3_boundary/<dir>/, else the project path
_P1_CAND = [os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'P1'),
            '/Users/martin/Documents/CL1/PRFull/check/ref3/eq3_boundary/P1']
DEFAULT_P1 = next(d for d in _P1_CAND if os.path.isfile(os.path.join(d, 'prove_strata.json')))
# the title and legend, part of the paper: PRFull/Text/si_box3_legend.tex, four levels up from check/ref3/eq3_boundary/box
DEFAULT_LEGEND = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), *(['..'] * 4), 'Text',
                                               'si_box3_legend.tex'))

a, b, c, u, v = sp.symbols('a b c u v')
LOC = {'a': a, 'b': b, 'c': c, 'u': u, 'v': v}
I_ = c * v - (1 - a)
II_ = (1 - a) * (u - 1) + b * v
# the open-facet display of the SI 'Partners' paragraph: gain = (k1 I + k2 II)/den
DISP = {
    '1001': (0, 0, 1), '1010': (0, 0, 1), '1011': (0, 0, 1), '1101': (0, 0, 1),
    '1110': (0, 0, 1), '1111': (0, 0, 1),
    '0000': (1, 0, 1 - a + c), '0100': (1, 0, 1 - a + c),
    '1000': (2, 0, 2 - 2 * a + 3 * c),
    '1100': (1, 0, 2 - 2 * a + 2 * c),
    '0110': (0, 1, 1 - a + 2 * b), '0111': (0, 1, 1 - a + 2 * b),
    '0001': (1, 1 - c, 2 - 2 * a + b + c - b * c),
    '0010': (1 - b, c, (1 - a) * (1 - b) + 2 * c),
    '0011': (1 - b, 1, 2 * (1 - a + c) + 2 * b * (a - c)),
    '0101': (b, 1 - c, (1 - a) * (1 - c) + b * (3 - 2 * a)),
}
SIG_TEX = {'0000': r'\ALLD', '1111': r'\ALLC', '1010': r'\TFT', '1001': r'\WSLS', '1000': r'\Grim'}


def S(x):
    return sp.sympify(x, locals=LOC)


# ----------------------------------------------------------------- printing, with every factor made positive

def poly_tex(f):
    """A polynomial in a, b, c: terms by degree, constant first; e.g. 1-a+2b, 2b+c-bc."""
    f = sp.expand(f)
    if f == 0:
        return '0'
    P = sp.Poly(f, a, b, c)
    terms = sorted(P.terms(), key=lambda t: (sum(t[0]), [-e for e in t[0]]))
    out = ''
    for k, (mon, cf) in enumerate(terms):
        cf = sp.Rational(cf)
        var = ''.join(('%s^%d' % (n, e) if e > 1 else n) for n, e in zip('abc', mon) if e > 0)
        sgn = '-' if cf < 0 else '+'
        mag = abs(cf)
        if var:
            cs = '' if mag == 1 else sp.latex(mag)
            body = cs + var
        else:
            body = sp.latex(mag)
        if k == 0:
            out = ('-' if sgn == '-' else '') + body
        else:
            out += sgn + body
    return out


def nterms(f):
    return len(sp.Add.make_args(sp.expand(f)))


def positive_parts(expr, centre):
    """expr = K * prod f_i^e_i with every f_i > 0 at `centre` (signs certified on the stratum in the JSON).
    Returns (K, [(f_i, e_i)])."""
    expr = sp.factor(sp.cancel(expr))
    if expr == 0:
        return 0, []
    K, facs = sp.factor_list(expr, a, b, c)
    K = sp.Rational(K)
    out = []
    for f, e in facs:
        val = f.subs(centre)
        assert val != 0, ('factor vanishes at the centre', f)
        if val < 0:
            f = -f
            if e % 2:
                K = -K
        out.append((sp.expand(f), e))
    out.sort(key=lambda t: fkey(t[0]))
    return K, out


def fkey(f):
    """Order of factors: by degree, the variables they contain, the number of terms, then the coefficients."""
    P = sp.Poly(f, a, b, c)
    terms = sorted(P.terms(), key=lambda t: (sum(t[0]), [-e for e in t[0]]))
    return (P.total_degree(), ''.join(sorted(str(x) for x in f.free_symbols)), len(terms),
            [(sum(m), [-e for e in m], abs(sp.Rational(cf))) for m, cf in terms])


def prod_tex(K, facs):
    """Product of positive factors with a positive rational constant K in front, no outer brackets."""
    assert K > 0
    parts = []
    for f, e in facs:
        s = poly_tex(f)
        if nterms(f) > 1:
            s = '(' + s + ')'
        if e > 1:
            s += '^%d' % e
        parts.append(s)
    body = ''.join(parts)
    if K == 1:
        return body or '1'
    ks = sp.latex(K) if K.q == 1 else r'\tfrac{%d}{%d}' % (K.p, K.q)
    return ks + body


def den_tex(Kd, dfacs):
    """'' (no denominator), or '/x' with the brackets of the SI display."""
    Kd = sp.Rational(Kd)
    assert Kd > 0 and Kd.q == 1
    if Kd == 1 and not dfacs:
        return ''
    d = prod_tex(Kd, dfacs)
    atom = (not dfacs) or (Kd == 1 and len(dfacs) == 1 and dfacs[0][1] == 1)
    return '/' + d if atom else '/[' + d + ']'


def gain_tex(l1, l2, centre, a_is_one, cval):
    """gain = l1 I + l2 II (mu = 0), or lam v where a = 1 (then I = c v, II = b v and l2 = 0)."""
    if a_is_one:
        lam = sp.cancel(l1 * cval)
        n, d = lam.as_numer_denom()
        Kn, nf = positive_parts(n, centre)
        Kd, df = positive_parts(d, centre)
        K = sp.Rational(Kn) / sp.Rational(Kd)
        assert K > 0
        num = prod_tex(sp.Integer(K.p), nf)
        num = 'v' if num == '1' else num + 'v'
        return num + den_tex(sp.Integer(K.q), df)
    l1 = sp.cancel(l1)
    l2 = sp.cancel(l2)
    n1, d1 = l1.as_numer_denom()
    n2, d2 = l2.as_numer_denom()
    if l1 != 0 and l2 != 0:
        D = sp.lcm(sp.expand(d1), sp.expand(d2))
    else:
        D = d1 if l1 != 0 else d2
    k1 = sp.expand(sp.cancel(l1 * D))
    k2 = sp.expand(sp.cancel(l2 * D))
    Kd, df = positive_parts(D, centre)
    if Kd < 0:
        Kd, k1, k2 = -Kd, -k1, -k2
    # a common rational content of k1, k2 moves into the denominator constant
    pieces = []
    parts = []
    for k, sym in ((k1, 'I'), (k2, r'I\!I')):
        if k == 0:
            continue
        Kk, kf = positive_parts(k, centre)
        assert Kk > 0, 'negative multiplier'
        parts.append((sp.Rational(Kk), kf, sym))
    g = sp.gcd_list([P[0] for P in parts]) if len(parts) > 1 else parts[0][0]
    Ktot = sp.Rational(g) / sp.Rational(Kd)
    for Kk, kf, sym in parts:
        Kk = Kk / g
        assert Kk.q == 1
        body = prod_tex(Kk, kf)
        pieces.append(sym if body == '1' else body + sym)
    num = '+'.join(pieces)
    if len(pieces) > 1:
        num = '[' + num + ']'
    if Ktot.p != 1:
        num = sp.latex(Ktot.p) + num
    return num + den_tex(sp.Integer(Ktot.q), df)


# ----------------------------------------------------------------- the strata

def stratum_values(label):
    code = label[3:-1].split(',')
    sub, free = {}, []
    for name, cd in zip('abc', code):
        if cd == '*':
            free.append(LOC[name])
        else:
            sub[LOC[name]] = sp.Integer(int(cd))
    return sub, free


def stratum_tex(label):
    code = label[3:-1].split(',')
    return '(1,' + ','.join(n if cd == '*' else cd for n, cd in zip('abc', code)) + ')'


def selfplay_lead(p1dir, label):
    """L by the prover's own routine; returns (k0, L, certs)."""
    sys.path.insert(0, p1dir)
    import prove_strata as ps                                   # noqa: E402  (the certificate program itself)
    code = tuple(label[3:-1].split(','))
    vals_r, _, free = ps.setup(code)
    p = (ps.R(1),) + tuple(vals_r)
    W, inf = ps.limit_W(p, p, free)
    assert [sp.simplify(w) for w in W] == [1, 0, 0, 0] and inf['lead_ok']
    assert inf['k0'] == 0     # the legend: L is the term of order eps^0, from the trees rooted at CC alone
    return inf['k0'], sp.sympify(str(inf['lead']), locals=LOC), inf['lead_certs']


def log_leads(p1dir):
    out = {}
    cur = None
    for line in open(os.path.join(p1dir, 'prove_strata.log')):
        m = re.match(r'STRATUM (\(\S+\))', line)
        if m:
            cur = m.group(1)
        m = re.search(r'self-play W\(p,p\) = \[1, 0, 0, 0\] ; eps-order (\d+), normaliser lead (.*), certified>0: True', line)
        if m and cur:
            out[cur] = (int(m.group(1)), S(m.group(2).replace('**', '^')))
            cur = None
    return out


def bernstein_factors(rows, lead_certs):
    facs = {}
    certs = list(lead_certs)
    for r in rows.values():
        certs += r['lead_certs'] + r['certs']['l1'] + r['certs']['l2'] + r['certs']['mu']
    for cs in certs:
        m = re.match(r'(num|den)\[(.*)\]\^(\d+):(.*)', cs)
        if not m:
            continue
        meth = m.group(4)
        assert meth.startswith('subst') or meth.startswith('manifest') or meth == 'const', cs
        if meth.startswith('subst'):
            facs[m.group(2)] = meth
    return facs


def read_legend(path):
    """(title, legend): the two lines of the legend file that are not comments; None if the file is absent."""
    if not os.path.isfile(path):
        return None
    lines = [l.rstrip('\n') for l in open(path) if l.strip() and not l.startswith('%')]
    assert len(lines) == 2, (path, 'expected a title line and a legend line')
    return lines[0], lines[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--p1', default=DEFAULT_P1)
    ap.add_argument('--style', choices=('box', 'table'), default='box',
                    help="'box': an SI Box in the style of SI Box 2 (default, no renumbering); 'table': an \\sitable")
    ap.add_argument('--label', default=None, help="default 'SI Box 3' (box) or 'SI Table 1' (table)")
    ap.add_argument('--out', default=None)
    ap.add_argument('--legend', default=None,
                    help='the title and legend (default: PRFull/Text/si_box3_legend.tex, found from this script)')
    args = ap.parse_args()
    if args.legend is not None and not os.path.isfile(args.legend):
        ap.error('--legend %s: no such file' % args.legend)
    if args.label is None:
        args.label = 'SI Box 3' if args.style == 'box' else 'SI Table 1'
    if args.out is None:
        args.out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                'si_box3_strata.tex' if args.style == 'box' else 'si_strata_table.tex')
    d = json.load(open(os.path.join(args.p1, 'prove_strata.json')))
    assert d['overall'] is True and d['open_facet_ok'] is True
    loglead = log_leads(args.p1)
    rows_tex = []
    n_diff = 0
    n_bern = 0
    seg_checked = False
    for e in d['strata']:
        lab = e['stratum']
        assert e['selfplay_cc'] and e['limit_valid'] and e['sufficiency'] and e['necessity'], lab
        sub, free = stratum_values(lab)
        centre = dict(sub)
        centre.update({x: sp.Rational(1, 2) for x in free})
        # centre must avoid zeros of every factor; 1/2 is used, and checked by positive_parts' assertion
        a_is_one = sub.get(a) == 1
        cval = sub.get(c, c)
        Is, IIs = I_.subs(sub), II_.subs(sub)
        # self-play normaliser
        k0, L, Lcerts = selfplay_lead(args.p1, lab)
        assert k0 == 0
        assert lab in loglead and loglead[lab][0] == 0 and sp.expand(loglead[lab][1] - L) == 0, lab
        KL, Lf = positive_parts(L, centre)
        assert KL > 0
        L_tex = prod_tex(KL, Lf)
        # in the table: one unbreakable box per factor, a break allowed between factors
        L_cell = r'\mbox{$L=' + (prod_tex(KL, Lf[:1]) if Lf else '1') + '$}'
        for f, ee in Lf[1:]:
            L_cell += r'\allowbreak\mbox{$' + prod_tex(sp.Integer(1), [(f, ee)]) + '$}'
        if KL != 1:
            L_cell = r'\mbox{$L=' + L_tex + '$}'
        if lab == '(1,*,*,1)':
            # the two excluded segments are where L vanishes on the closure
            assert sp.expand(L - (1 + b - a) * (1 - a - b + 2 * a * b)) == 0
            assert L.subs({a: 0, b: 1}) == 0 and L.subs({a: 1, b: 0}) == 0
            seg_checked = True
        # gains
        diffs = []
        for s in sorted(e['rows']):
            r = e['rows'][s]
            l1, l2, mu = S(r['l1']), S(r['l2']), S(r['mu'])
            assert mu == 0, (lab, s)
            assert r['ok'] and r['identity'] and r['lead_ok'], (lab, s)
            if a_is_one:
                assert r['W_CD_zero(a=1)'] is True and S(r['W'][1]) == 0 and S(r['W'][3]) == 0, (lab, s)
                assert l2 == 0
            g = sp.cancel(l1 * Is + l2 * IIs)
            k1, k2, den = DISP[s]
            dens = sp.sympify(den).subs(sub)
            mark = ''
            if sp.simplify(dens) == 0:
                same = False
                mark = r'$^\dagger$'
            else:
                go = (sp.sympify(k1).subs(sub) * Is + sp.sympify(k2).subs(sub) * IIs) / dens
                same = sp.simplify(g - go) == 0
            if s == '0000':
                assert sp.simplify(g - Is / (1 - a + c).subs(sub)) == 0, lab
            if s == '0111':
                assert sp.simplify(g - IIs / (1 - a + 2 * b).subs(sub)) == 0, lab
            if not same:
                assert g != 0
                name = (SIG_TEX[s] + '{}') if s in SIG_TEX else '$%s$' % s
                diffs.append(r'%s%s:~$%s$' % (name, mark, gain_tex(l1, l2, centre, a_is_one, cval)))
        n_diff += len(diffs)
        bf = bernstein_factors(e['rows'], Lcerts)
        bf_tex = []
        for f in bf:
            fe = S(f.replace('^', '**'))
            if fe.subs(centre) < 0:
                fe = -fe
            bf_tex.append((fkey(fe), poly_tex(fe)))
        bf_tex.sort()
        n_bern += len(bf_tex)
        rows_tex.append((lab, stratum_tex(lab), L_cell, diffs, [t[1] for t in bf_tex]))
    assert seg_checked
    # order: a = 0, a = 1, a free
    # JSON order (a, b in 0, 1, free; c free, then 1), which already groups a = 0, a = 1, a free
    # the open facet: L = c(1-a+b)(2ab+c-ac-bc); with c > 0 its zeros on the closed facet are the two segments
    k0, Lo, _ = selfplay_lead(args.p1, '(1,*,*,*)')
    assert k0 == 0 and sp.expand(Lo - c * (1 - a + b) * (2 * a * b + c - a * c - b * c)) == 0
    # 2ab + c(1-a-b) >= 0 on [0,1]^3 with c > 0 vanishes only at (a,b) = (0,1), (1,0): for a+b <= 1 both terms are
    # non-negative; for a+b > 1 it is >= 2ab + 1-a-b = ab + (1-a)(1-b) > 0.  Spot check on a grid:
    for i in range(11):
        for j in range(11):
            for k in range(1, 11):
                aa, bb, cc = sp.Rational(i, 10), sp.Rational(j, 10), sp.Rational(k, 10)
                val = Lo.subs({a: aa, b: bb, c: cc})
                assert (val == 0) == ((aa, bb) in ((0, 1), (1, 0))) and val >= 0
    # the title and legend are the paper's (Text/si_box3_legend.tex); every statement of the legend is asserted above,
    # whether or not the file is here.  head: the lines of L that carry them, which the body leaves out
    tl = read_legend(args.legend or DEFAULT_LEGEND)
    title, legend = tl if tl else ('', '')
    L, head = [], []
    L.append('%% %s -- written by mk_strata_table.py from the certificate output prove_strata.json (P1); do not edit by hand.' % args.label)
    if args.style == 'box':
        L.append(r'\begin{table}[!tp]')
        L.append(r'\noindent\rule{\textwidth}{0.4pt}\par\smallskip')
        head.append(len(L))
        L.append(r'\noindent{\bf %s: %s}\par\smallskip' % (args.label, title))
        L.append(r'{\footnotesize\renewcommand{\baselinestretch}{1.1}\selectfont')
        head.append(len(L))
        L.append(r'\noindent ' + legend + r'\par\smallskip')
        L.append(r'\noindent\setlength{\tabcolsep}{4pt}%')
    else:
        L.append(r'\begin{table}[!tp]\footnotesize\let\small\footnotesize\renewcommand{\baselinestretch}{1.1}\selectfont')
        head.append(len(L))
        L.append(r'\sitable{%s: %s}{%s}' % (args.label, title, legend))
        L.append(r'\begin{center}\setlength{\tabcolsep}{4pt}')
    L.append(r'\begin{tabular}{@{}>{\raggedright\arraybackslash}p{0.2\textwidth}>{\raggedright\arraybackslash}p{0.48\textwidth}>{\raggedright\arraybackslash}p{0.28\textwidth}@{}}')
    L.append(r'\toprule')
    L.append(r'stratum; $L$ & gains that differ from the open-facet display & Bernstein factors\\')
    prev = None
    for lab, st, Lt, diffs, bfs in rows_tex:
        grp = lab[3]
        if prev is None or grp != prev:
            L.append(r'\midrule')
        prev = grp
        g = '; '.join(diffs) if diffs else 'none'
        bt = ', '.join('$%s$' % x for x in bfs) if bfs else 'none'
        L.append(r'$%s$\newline %s & %s & %s\\' % (st, Lt, g, bt))
    L.append(r'\bottomrule')
    L.append(r'\end{tabular}')
    if args.style == 'box':
        L.append(r'\par}')
        L.append(r'\noindent\rule{\textwidth}{0.4pt}')
    else:
        L.append(r'\end{center}')
    L.append(r'\end{table}')
    # the body, written in every case: the box without its title and legend, beside --out
    body = ['%% %s without its title and legend (they are part of the paper) -- written by mk_strata_table.py from the '
            'certificate output prove_strata.json (P1); do not edit by hand.' % args.label]
    body += [x for k, x in enumerate(L) if k > 0 and k not in head]
    body_out = '%s_body%s' % os.path.splitext(args.out)
    with open(body_out, 'w') as fh:
        fh.write('\n'.join(body) + '\n')
    # file names without their directories, so that the output does not depend on where the files are written
    print('wrote %s: 13 strata, %d gains that differ from the display, %d Bernstein factors'
          % (os.path.basename(body_out), n_diff, n_bern))
    if tl is None:
        print('%s not written: its title and legend (Text/si_box3_legend.tex) are part of the paper, not of this release'
              % os.path.basename(args.out))
        return
    with open(args.out, 'w') as fh:
        fh.write('\n'.join(L) + '\n')
    print('wrote %s: the same with its title and legend' % os.path.basename(args.out))


if __name__ == '__main__':
    main()
