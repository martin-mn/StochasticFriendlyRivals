"""Prover 2 -- symbolic proof on the 13 boundary strata (and re-derivation of the open facet).

For each stratum (codes for a, b, c in {0, 1, *}, * = open (0,1)) and each of the 16 binary
deviators sigma:
  1. W(sigma, p) by the limit lemma of p2lib (closed classes of M0, stationary laws on them,
     absorption probabilities, first-order reduced generator Q) -- no cofactors.
  2. Every object the lemma needs positive is certified positive on the WHOLE open stratum
     (mu_i entries, denominators of absorption probabilities, positive Q entries, w entries,
     denominators of W), by factoring into irreducibles and the x = t/(1+t) coefficient test.
  3. Gain g = (u-1) W_CD + v W_DC - W_DD; decomposition g = l1 I + l2 II - mu with
     l1, l2, mu certified >= 0 (and the identity checked symbolically in u, v).
  4. Self-play: W(p,p) = (1,0,0,0) by the same lemma (unique closed class {CC}).
  5. Necessity witnesses: a deviator with g = k1 * I, k1 > 0 and one with g = k2 * II, k2 > 0
     (for a = 1: g = k v, k > 0, since then I = c v, II = b v with b, c > 0).
Run:  python3 p2_proof.py      (about 10-20 s)
"""
import itertools, json, sys, time
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from p2lib import *

U, V = sp.symbols('u v', real=True)

STRATA = []
for ca, cb, cc in itertools.product('01*', '01*', '1*'):
    if (ca, cb, cc) == ('*', '*', '*'):
        continue
    if (ca, cb) in (('0', '1'), ('1', '0')):
        continue
    STRATA.append(ca + cb + cc)
assert len(STRATA) == 13

# the open-facet display of the SI 'Partners' paragraph, as gain = expression in I, II
Isym, IIsym = sp.symbols('I II')
DISPLAY = {}
for s in ('1001', '1010', '1011', '1101', '1110', '1111'):
    DISPLAY[s] = sp.Integer(0)
for s in ('0000', '0100'):
    DISPLAY[s] = Isym / (1 - A + C)
DISPLAY['1000'] = 2 * Isym / (2 - 2 * A + 3 * C)
DISPLAY['1100'] = Isym / (2 - 2 * A + 2 * C)
for s in ('0110', '0111'):
    DISPLAY[s] = IIsym / (1 - A + 2 * B)
DISPLAY['0001'] = (Isym + (1 - C) * IIsym) / (2 - 2 * A + B + C - B * C)
DISPLAY['0010'] = ((1 - B) * Isym + C * IIsym) / ((1 - A) * (1 - B) + 2 * C)
DISPLAY['0011'] = ((1 - B) * Isym + IIsym) / (2 * (1 - A + C) + 2 * B * (A - C))
DISPLAY['0101'] = (B * Isym + (1 - C) * IIsym) / ((1 - A) * (1 - C) + B * (3 - 2 * A))


def cert(expr, free, want):
    """want '+' (>0), '0+' (>=0: identically 0 or >0) or 'nz' (sign-definite, never 0) on the
    whole open stratum; returns the certificate string or raises."""
    s, rep = factor_sign_cert(expr, free)
    ok = {'+': s == '+', '0+': s in ('0', '+'), 'nz': s in ('+', '-')}[want]
    if not ok:
        raise AssertionError('certificate failed for %s (want %s): %s %s' % (expr, want, s, rep))
    if s == '0':
        return '0'
    return ('sign %s = ' % s) + (' * '.join('(%s)^%d:%s[%s]' % (f, m, sg, how) for f, m, sg, how in rep)
                                 or 'nonzero constant')


def run():
    t0 = time.time()
    report, summary = [], {}
    all_polys = []  # (stratum, role, polynomial) -- for the independent root search
    for code in STRATA + ['***']:
        p_, free = stratum_values(code)
        a, b, c = p_
        p = (sp.Integer(1), a, b, c)
        I = c * V - (1 - a)
        II = (1 - a) * (U - 1) + b * V
        report.append('=' * 100)
        report.append('STRATUM (a,b,c) = (%s,%s,%s)   free = %s' % (code[0], code[1], code[2], free))
        # ---- self-play
        Wpp, ipp = limit_W(p, p, free)
        selfplay = (ipp['classes'] == [[0]] and Wpp == [1, 0, 0, 0])
        report.append('  self-play: closed classes of M0(p,p) = %s ; W(p,p) = %s ; permanent CC: %s'
                      % ([[NAMES[s] for s in cl] for cl in ipp['classes']], Wpp, selfplay))
        assert selfplay
        rows = {}
        nmax_cls = 1
        for sg in SIGMAS:
            sig = tuple(sp.Integer(int(x)) for x in sg)
            W, info = limit_W(sig, p, free)
            nmax_cls = max(nmax_cls, info['k'])
            # --- certify positivity of every lemma object on the open stratum
            certs = []
            for role, e in info['pos']:
                certs.append((role, cert(e, free, 'nz' if role in ('h-den', 'stat-det', 'trans-det') else '+')))
                all_polys += [(code, sg, role + '-num', f) for f in [sp.fraction(sp.together(e))[0]]]
                all_polys += [(code, sg, role + '-den', sp.fraction(sp.together(e))[1])]
            for role, e in info['zero_or_pos']:
                certs.append((role, cert(e, free, '0+')))
                all_polys += [(code, sg, role + '-den', sp.fraction(sp.together(e))[1])]
            # --- consistency of the limit: W is a law, stationary for M0, and solves the
            #     first-order equations  W M1 h_j = 0
            M0, M1, M2 = expansion(sig, p)
            Wv = sp.Matrix([W])
            assert sp.simplify(sum(W) - 1) == 0
            assert (Wv * M0 - Wv).applyfunc(sp.simplify) == sp.zeros(1, 4)
            if info['k'] > 1:
                for hj in info['h']:
                    hv = sp.Matrix([hj[s] for s in range(4)])
                    assert sp.simplify((Wv * M1 * hv)[0]) == 0
            for Wx in W:
                if Wx != 0:
                    cert(Wx, free, '+')
                    all_polys.append((code, sg, 'W-den', sp.fraction(sp.together(Wx))[1]))
            # --- gain and decomposition
            g = sp.expand((U - 1) * W[1] + V * W[2] - W[3])
            if a != 1:
                l2 = sp.factor(sp.cancel(W[1] / (1 - a)))
                l1 = sp.factor(sp.cancel((W[2] - b * l2) / c))
                mu = sp.factor(sp.cancel(W[3] - (1 - a) * l1))
            else:
                assert sp.simplify(W[1]) == 0, ('a=1 but W_CD != 0', code, sg, W)
                l2 = sp.Integer(0)
                l1 = sp.factor(sp.cancel(W[2] / c))
                mu = sp.factor(sp.cancel(W[3]))
            assert sp.simplify(g - (l1 * I + l2 * II - mu)) == 0
            cl1, cl2, cmu = cert(l1, free, '0+'), cert(l2, free, '0+'), cert(mu, free, '0+')
            for nm, e in (('l1', l1), ('l2', l2), ('mu', mu)):
                if e != 0:
                    all_polys.append((code, sg, nm + '-num', sp.fraction(sp.together(e))[0]))
                    all_polys.append((code, sg, nm + '-den', sp.fraction(sp.together(e))[1]))
            rows[sg] = dict(W=W, l1=l1, l2=l2, mu=mu, classes=info['classes'], k=info['k'],
                            Q=info['Q'], w=info['w'])
            if code == '***':
                d = DISPLAY[sg].subs({Isym: I, IIsym: II})
                assert sp.simplify(g - d) == 0, ('open facet display mismatch', sg)
                # the display AS WRITTEN: denominator > 0, coefficients of I, II >= 0, no constant
                dn = sp.fraction(sp.together(DISPLAY[sg]))
                if DISPLAY[sg] != 0:
                    dden = sp.fraction(sp.together(DISPLAY[sg]))[1]
                    cI = sp.diff(DISPLAY[sg], Isym); cII = sp.diff(DISPLAY[sg], IIsym)
                    assert sp.simplify(DISPLAY[sg] - cI * Isym - cII * IIsym) == 0
                    report.append('        display: den {%s}  coef(I) {%s}  coef(II) {%s}' % (
                        cert(dden, free, '+'), cert(cI, free, '0+'), cert(cII, free, '0+')))
            report.append('  %s  classes=%-30s W=%s' % (
                sg, [[NAMES[s] for s in cl] for cl in info['classes']], [sp.factor(x) for x in W]))
            if info['k'] > 1:
                report.append('        Q=%s  w=%s' % (info['Q'].tolist(), info['w']))
            report.append('        g = l1*I + l2*II - mu :  l1=%s  l2=%s  mu=%s' % (l1, l2, mu))
            report.append('        certs: l1 {%s}  l2 {%s}  mu {%s}' % (cl1, cl2, cmu))
            report.append('        lemma objects: %s' % ('; '.join('%s{%s}' % rc for rc in certs)))
        # ---- necessity witnesses
        wit = {}
        if a != 1:
            wI = [s for s, r in rows.items() if r['l2'] == 0 and r['mu'] == 0 and r['l1'] != 0]
            wII = [s for s, r in rows.items() if r['l1'] == 0 and r['mu'] == 0 and r['l2'] != 0]
            assert wI and wII, (code, wI, wII)
            for s in wI:
                cert(rows[s]['l1'], free, '+')
            for s in wII:
                cert(rows[s]['l2'], free, '+')
            wit = {'I': {s: str(rows[s]['l1']) for s in wI}, 'II': {s: str(rows[s]['l2']) for s in wII}}
        else:
            # I = c v, II = b v, both > 0 iff v > 0; witness: g = k v, k > 0
            wv = [s for s, r in rows.items() if r['mu'] == 0 and r['l1'] != 0]
            for s in wv:
                cert(rows[s]['l1'], free, '+')
            assert wv
            wit = {'v (I=cv, II=bv)': {s: str(sp.factor(rows[s]['l1'] * c)) for s in wv}}
        report.append('  NECESSITY witnesses (gain = positive multiple of I, resp. II): %s' % wit)
        mus_pos = [s for s, r in rows.items() if r['mu'] != 0]
        summary[code] = dict(selfplay=selfplay, max_classes=nmax_cls, witnesses=wit,
                             mu_nonzero=mus_pos,
                             rows={s: dict(W=[str(x) for x in r['W']], l1=str(r['l1']), l2=str(r['l2']),
                                           mu=str(r['mu']), classes=r['classes'], k=r['k'])
                                   for s, r in rows.items()})
        print('stratum %s done (%.1fs); max #closed classes %d; witnesses %s; mu != 0 for %s'
              % (code, time.time() - t0, nmax_cls, wit, mus_pos))
        sys.stdout.flush()
    return report, summary, all_polys


if __name__ == '__main__':
    import os, pickle
    here = os.path.dirname(os.path.abspath(__file__))
    report, summary, all_polys = run()
    open(os.path.join(here, 'p2_proof_report.txt'), 'w').write('\n'.join(report) + '\n')
    json.dump(summary, open(os.path.join(here, 'p2_summary.json'), 'w'), indent=1, default=str)
    pickle.dump([(c, s, r, sp.srepr(f)) for c, s, r, f in all_polys],
                open(os.path.join(here, 'p2_polys.pkl'), 'wb'))
    print('ALL CERTIFICATES PASSED; %d polynomials recorded for root search' % len(all_polys))
