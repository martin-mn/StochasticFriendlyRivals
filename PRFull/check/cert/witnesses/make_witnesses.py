#!/usr/bin/env python3
r"""Machine-readable archive of every exact memory-two witness of the paper (ms.tex), recomputed from
scratch in exact rational arithmetic.  One JSON per witness in this directory, the exact margins
against all 65536 binary co-players in margins/<name>.tsv.gz, and index.json.

    cd PRFull/check && nice -n 10 python3 cert/witnesses/make_witnesses.py            # all (about 5 min, 3 processes)
    cd PRFull/check && nice -n 10 python3 cert/witnesses/make_witnesses.py W_SI2 N_SI2  # some

Per witness (strategy p given as exact rationals, wedge -> sign of T-S and the efficiency class):
  selfplay     exact limit (w_CC, w_CD, w_DC, w_DD) of p against itself (cert/exact.py, mpq), the L1
               efficiency class of p's stratum (integer exponents only), and whether p is efficient;
  margins      the exact margin sign*(w_DC - w_CD) against EVERY binary co-player g (genome, bit s =
               q_s), by a full exact scan -- no floating point anywhere;
  categories   of the co-players:
                 structural_zero   w_CD = w_DC = 0 by exponents alone (on the whole open stratum, L1)
                 identity          exact 0, not structural, and N_q = D* margin is identically 0 on the
                                   open stratum (every stratum corner with E_c = E has corner margin 0,
                                   L3; decided by the interval corner scan with integrality);
                                   for a binary point the stratum is the point, so every exact
                                   non-structural zero is an identity
                 pointwise_zero    exact 0 at p but N_q not identically 0 on the stratum
                 identity_undetermined   exact 0, corner scan inconclusive (not expected)
                 positive / negative
  min_margin   the exact minimum, the co-players attaining it (argmin, by category when it is 0),
               and the smallest positive margin with its co-players;
  cross-checks the interval route (cert.certify.certify_point: outward-rounded enclosures, structural
               zeros, exact rationals on the unresolved) must give the same counts; a sample of up to
               16 identities is re-verified exactly on the 2-level grid {1/3,2/3}^d of the stratum
               (a multilinear polynomial vanishing there is 0), when 2^d <= 2048;
  friendly_rival   efficient and no negative margin (competitiveness over all co-players then follows
               from the box-vertex lemma).
"""
import os, sys, json, time, gzip, hashlib, itertools, datetime
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
CHECK = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, CHECK); sys.path.insert(0, os.path.join(CHECK, 'families'))
from fractions import Fraction as F
import numpy as np
import famlib as fl
from cert import engine as en, exact as ex, certify as ce

H, Q, T3 = F(1, 2), F(1, 4), F(1, 3)
MS = 'PRFull/Text/ms.tex'

def Pt(t): return [0, 0, 0, t, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0]
def Pab(a, b): return [1, 1, 1, 0, 0, 0, 1, 0, 1, a, 1, b, 0, 0, 0, 0]
def centre_of(pat): return [H if c == '*' else F(int(c)) for c in pat]

WITNESSES = [
    # name, strategy, wedge, source, paper claim (stratum, dim, selfplay, min margin, #co-players at min)
    ('W_SI2', [1, H, 1, H, 0, 0, 1, H, 1, Q, 1, Q, Q, 0, Q, 0], 'W',
     MS + ' eq:certW (main text, "The Prisoner\'s Dilemma wedge") and SI Table 2 row W; M2f check/certify7.py',
     dict(stratum='1*1*001*1*1**0*0', dim=7, selfplay=['1', '0', '0', '0'], min_margin='0', n_at_min=37996)),
    ('S_SI2', [1, H, H, H, 1, H, H, H, 1, H, H, 1, 0, H, H, H], 'S',
     MS + ' SI Table 2 row S; M2f check/certify11S.py',
     dict(stratum='1***1***1**10***', dim=11, selfplay=['1', '0', '0', '0'], min_margin='0', n_at_min=30850)),
    ('N_SI2', [F(2, 3), 1, 0, T3, T3, T3, 0, 0, F(2, 3), 1, F(2, 3), T3, T3, 0, 1, 0], 'N',
     MS + ' SI Table 2 row N; M2f check/certify8N.py (corrected copy PRFull/check/certify8N.py)',
     dict(stratum='*10***00*1***010', dim=8, selfplay=['0', '1/2', '1/2', '0'], min_margin='0', n_at_min=31024)),
    ('E_SI2', [1, 0, 1, F(2, 3), F(2, 3), T3, 0, T3, 1, 1, F(2, 3), F(2, 3), F(2, 3), 1, 0, T3], 'E',
     MS + ' SI Table 2 row E (the C/D mirror image of row N)',
     dict(stratum='101***0*11***10*', dim=8, selfplay=['0', '1/2', '1/2', '0'], min_margin='0', n_at_min=31024)),
]
for tt in (F(1, 10), F(1, 4), F(1, 2), F(3, 4), F(9, 10)):
    WITNESSES.append(('pt_t%s' % str(tt).replace('/', '_'), Pt(tt), 'N',
                      MS + ' main text "binary-invisible" p(t), Figure 3c, SI "The other wedges"; M2f check/certify_segment.py '
                      '(the whole open segment is certified by PRFull/check/families/cert_pt.py)',
                      dict(stratum='000*000001011010', dim=1, selfplay=['0', '1/2', '1/2', '0'],
                           **({'min_margin': '0'} if tt == H else {}))))
WITNESSES.append(('pt_t0', Pt(0), 'N', MS + ' p(t) at t = 0: "not even efficient" (Figure 3c)',
                  dict(stratum='0000000001011010', dim=0, selfplay=['1/6', '1/3', '1/3', '1/6'], efficient=False)))
WITNESSES.append(('pt_t1', Pt(1), 'N', MS + ' p(t) at t = 1: "not even efficient" (Figure 3c)',
                  dict(stratum='0001000001011010', dim=0, selfplay=['2/9', '1/3', '1/3', '1/9'], efficient=False)))
WITNESSES.append(('fig3a_left_centre', centre_of('1*1*0*0101*1*010'), 'W',
                  MS + ' Figure 3a left: the maximal stratum of the binary-anchored component (the 5-cube), its centre',
                  dict(stratum='1*1*0*0101*1*010', dim=5)))
WITNESSES.append(('fig3a_right_centre', centre_of('1*1*001*1*1**0*0'), 'W',
                  MS + ' Figure 3a right: the maximal 7-dimensional stratum of the binary-invisible component, its centre '
                  '(the paper\'s witness in this stratum is W_SI2, not the centre)',
                  dict(stratum='1*1*001*1*1**0*0', dim=7)))
WITNESSES.append(('ab_a1_2_b1_2', Pab(H, H), 'W', MS + ' SI "The wedge W": the family p(a,b), extreme point a = b = 1/2 '
                  '(the whole family is certified by PRFull/check/families/cert_ab.py)', dict(stratum='111000101*1*0000', dim=2)))
WITNESSES.append(('ab_a1_4_b1_4', Pab(Q, Q), 'W', MS + ' SI "The wedge W": the family p(a,b) at a = b = 1/4',
                  dict(stratum='111000101*1*0000', dim=2)))
WITNESSES.append(('ab_sliver', Pab(H + F(1, 2**40), Q), 'W',
                  'negative control: p(a,b) at a = 1/2 + 2^-40, b = 1/4, accepted by the float census rule (>= -1e-9) '
                  'but not competitive (triage report 3)', dict(stratum='111000101*1*0000', dim=2, friendly_rival=False)))
for bits in itertools.product((0, 1), repeat=3):      # 111*000101*1*010: p3, p10, p12
    p = [1, 1, 1, bits[0], 0, 0, 0, 1, 0, 1, bits[1], 1, bits[2], 0, 1, 0]
    WITNESSES.append(('MB8_%d%d%d' % bits, p, 'W',
                      MS + ' SI "The wedge W": the eight binary friendly rivals 111*000101*1*010 of Murase and Baek '
                      '(p3, p10, p12 = %d, %d, %d)' % bits, dict(dim=0, selfplay=['1', '0', '0', '0'], friendly_rival=True)))
_S15 = [32, 27, 15, 25, 29, 28, 9, 6, 21, 14, 27, 15, 6, 12, 30, 22]
WITNESSES.append(('S15_box_centre', [F(n, 32) for n in _S15], 'S',
                  'NEW (not yet in the paper): centre of the certified 15-dimensional box of wedge-S friendly rivals '
                  '(cert.certify.certify_ratio_box, half-width 1/1024; triage scratch cert/certS_ratio.py)',
                  dict(stratum='1' + '*'*15, dim=15)))

CAT = ('structural_zero', 'identity', 'pointwise_zero', 'identity_undetermined', 'positive', 'negative')
CODE = dict(structural_zero='S', identity='I', pointwise_zero='Z', identity_undetermined='U', positive='P', negative='N')

def dumps(obj):
    """JSON with indentation, but every list of integers on one line"""
    import re
    txt = json.dumps(obj, indent=1)
    return re.sub(r'\[\s*(-?\d+(?:,\s*-?\d+)*)\s*\]', lambda mo: '[' + ','.join(mo.group(1).split(',')).replace(' ', '').replace('\n', '') + ']', txt)

def mirror_CD(p):
    """the C/D exchange: p'_s = 1 - p_{15-s}"""
    return [1 - F(p[15 - s]) for s in range(16)]

def build(name, p, wedge, source, claim):
    t0 = time.time()
    p = [F(x) for x in p]; sign, eff = en.WEDGES[wedge]
    pat = en.pattern_str(en.pattern_of(p)); fr = en.free_coords(pat); d = len(fr)
    # ---- self-play and efficiency
    Wsp = ex.selfplay(p); cls = en.efficiency_class(pat)
    efficient = en.meets(cls, eff)
    exact_eff = {'C': Wsp[0] == 1, 'A': Wsp[1] + Wsp[2] == 1, 'Z': Wsp[3] == 0}[eff]
    assert efficient == exact_eff, (name, cls, Wsp)
    # ---- full exact scan
    R = fl.exact_scan([p], sign)
    m = {g: R[g][0]['m'] for g in range(fl.NQ)}
    supp = {g: R[g][0]['supp'] for g in range(fl.NQ)}
    struct = [g for g in range(fl.NQ) if fl.structural_zero(supp[g])]
    sset = set(struct)
    for g in struct: assert m[g] == 0
    zns = [g for g in range(fl.NQ) if m[g] == 0 and g not in sset]
    pos = [g for g in range(fl.NQ) if m[g] > 0]; neg = [g for g in range(fl.NQ) if m[g] < 0]
    # ---- identities among the non-structural zeros (L3 corner scan over the stratum cube)
    ident, pwz, und = [], [], []
    if zns:
        if d == 0:
            ident = list(zns)
        elif ce.corner_cost_seconds(pat, nq=len(zns)) < 300:
            C = en.corner_scan(pat, sign, np.array(zns, dtype=np.int64), nthreads=en.MAX_THREADS)
            cnt = C['cnt']
            for i, g in enumerate(zns):
                c = cnt[i]
                if c[en.C_POS] > 0 or c[en.C_NEG] > 0: pwz.append(g)
                elif c[en.C_LOW] == 0 and c[en.C_UNR] == 0: ident.append(g)
                else: und.append(g)
        else:
            und = list(zns)
    cats = dict(structural_zero=struct, identity=ident, pointwise_zero=pwz, identity_undetermined=und, positive=pos, negative=neg)
    assert sum(len(v) for v in cats.values()) == fl.NQ
    # exact grid spot check of identities (and of pointwise zeros being non-identities)
    grid_check = None
    if d and (1 << d) <= 2048 and (ident or pwz):
        grid = []
        for bits in itertools.product((T3, 2*T3), repeat=d):
            x = list(p)
            for j, v in zip(fr, bits): x[j] = v
            grid.append(x)
        samp_i = ident[::max(1, len(ident)//16)][:16]; samp_z = pwz[::max(1, len(pwz)//8)][:8]
        Rg = fl.exact_scan(grid, sign, qlist=samp_i + samp_z)
        bad_i = sum(1 for g in samp_i if any(v['N'] != 0 for v in Rg[g]))
        bad_z = sum(1 for g in samp_z if all(v['N'] == 0 for v in Rg[g]))
        grid_check = dict(identities_checked=len(samp_i), identities_failing=bad_i, pointwise_zeros_checked=len(samp_z),
                          pointwise_zeros_that_vanish_on_grid=bad_z, grid_points=len(grid), ok=(bad_i == 0 and bad_z == 0))
    # ---- minimum and arg-min
    mn = min(m.values())
    if mn == 0:
        argmin = {k: cats[k] for k in ('structural_zero', 'identity', 'pointwise_zero', 'identity_undetermined')}
    else:
        argmin = {'at_min': sorted(g for g in range(fl.NQ) if m[g] == mn)}
    n_at_min = sum(1 for v in m.values() if v == mn)
    mpos = min((m[g] for g in pos), default=None)
    argpos = sorted(g for g in pos if m[g] == mpos) if mpos is not None else []
    # ---- interval route
    cp = ce.certify_point(p, wedge=wedge, nproc=en.MAX_THREADS)
    interval_agrees = (cp['n_pos'] == len(pos) and cp['n_neg'] == len(neg) and cp['n_struct_zero'] == len(struct)
                       and cp['n_exact_zero'] == len(zns) and tuple(cp['selfplay']) == tuple(Wsp))
    fr_ok = bool(efficient and not neg)
    # ---- margins file
    os.makedirs(os.path.join(HERE, 'margins'), exist_ok=True)
    cat_of = {}
    for k, v in cats.items():
        for g in v: cat_of[g] = CODE[k]
    lines = ''.join('%d\t%s\t%s\n' % (g, m[g], cat_of[g]) for g in range(fl.NQ)).encode()
    sha = hashlib.sha256(lines).hexdigest()
    mfile = os.path.join('margins', name + '.tsv.gz')
    with gzip.GzipFile(os.path.join(HERE, mfile), 'wb', mtime=0) as f: f.write(lines)
    # ---- paper claim comparison
    got = dict(stratum=pat, dim=d, selfplay=fl.vec_str(Wsp), min_margin=fl.fstr(mn), n_at_min=n_at_min,
               efficient=bool(efficient), friendly_rival=fr_ok)
    mism = {k: dict(paper=v, recomputed=got[k]) for k, v in claim.items() if got.get(k) != v}
    out = dict(
        name=name, source=source, strategy=fl.vec_str(p), stratum=pat, dimension=d, free_coordinates=fr,
        wedge=wedge, sign_T_minus_S=sign, required_efficiency_class=eff,
        conventions="p_s = Pr(cooperate) at state s = 4x+y, x the last outcome, y the one before, CC=0 CD=1 DC=2 DD=3, "
                    "focal action first; errors p -> eps+(1-2eps)p on both players, eps -> 0 before any inequality; "
                    "co-player code g in 0..65535 with bit s = q_s (1 = C); margin = sign*(w_DC - w_CD)",
        selfplay=dict(W=fl.vec_str(Wsp), efficiency_class_L1=cls, efficient=bool(efficient)),
        competitiveness=dict(
            n_coplayers=fl.NQ, min_margin=fl.fstr(mn), min_margin_float=float(mn), n_at_min=n_at_min,
            counts={k: len(v) for k, v in cats.items()}, argmin=argmin,
            min_positive_margin=fl.fstr(mpos) if mpos is not None else None,
            min_positive_margin_float=float(mpos) if mpos is not None else None, argmin_positive=argpos,
            negative_coplayers=neg if len(neg) <= 5000 else ('%d co-players, listed in the margins file' % len(neg))),
        identity_grid_check=grid_check,
        interval_route=dict(n_pos=cp['n_pos'], n_neg=cp['n_neg'], n_struct_zero=cp['n_struct_zero'],
                            n_exact_zero=cp['n_exact_zero'], n_unresolved_by_interval=cp['n_unresolved_by_interval'],
                            agrees_with_exact_scan=bool(interval_agrees)),
        friendly_rival=fr_ok,
        paper_claim=claim, recomputed=got, matches_paper=(not mism), mismatches=mism,
        margins_file=mfile, margins_sha256=sha,
        regenerate="cd /Users/martin/Documents/CL1/PRFull/check && nice -n 10 python3 cert/witnesses/make_witnesses.py " + name,
        engine="PRFull/check/cert (exact.py: GTH in leading-order exact rationals, gmpy2.mpq; engine.py/certeng.c: "
               "interval engine with outward rounding)",
        generated=datetime.date.today().isoformat(), seconds=round(time.time() - t0, 1))
    if name == 'E_SI2':
        out['mirror_of_N_SI2'] = fl.vec_str(mirror_CD([F(x) for x in dict((w[0], w[1]) for w in WITNESSES)['N_SI2']])) == out['strategy']
    with open(os.path.join(HERE, name + '.json'), 'w') as f: f.write(dumps(out))
    print("%-20s %-16s d=%-2d %s  self-play %-26s class %s eff %-5s  min %-14s at %6d  %s  FR %-5s  interval agrees %s  paper %s  grid %s  (%.1f s)"
          % (name, pat, d, wedge, '(' + ','.join(fl.vec_str(Wsp)) + ')', cls, efficient, fl.fstr(mn)[:14], n_at_min,
             {CODE[k]: len(v) for k, v in cats.items() if v}, fr_ok, interval_agrees, 'ok' if not mism else mism,
             None if grid_check is None else grid_check['ok'], time.time() - t0), flush=True)
    return out

def main(names):
    T = fl.timer(); idx = []
    todo = [w for w in WITNESSES if not names or w[0] in names]
    for w in todo:
        o = build(*w)
        idx.append(dict(name=o['name'], file=o['name'] + '.json', stratum=o['stratum'], dimension=o['dimension'], wedge=o['wedge'],
                        selfplay=o['selfplay']['W'], efficient=o['selfplay']['efficient'], min_margin=o['competitiveness']['min_margin'],
                        n_at_min=o['competitiveness']['n_at_min'], counts=o['competitiveness']['counts'],
                        friendly_rival=o['friendly_rival'], matches_paper=o['matches_paper'],
                        interval_agrees=o['interval_route']['agrees_with_exact_scan'], source=o['source']))
    if not names:
        with open(os.path.join(HERE, 'index.json'), 'w') as f:
            json.dump(dict(generated=datetime.date.today().isoformat(), n=len(idx), seconds=round(T(), 1),
                           regenerate="cd /Users/martin/Documents/CL1/PRFull/check && nice -n 10 python3 cert/witnesses/make_witnesses.py",
                           families=dict(p_t="PRFull/check/families/cert_pt.py -> cert_pt.json",
                                         p_ab="PRFull/check/families/cert_ab.py -> cert_ab.json"),
                           witnesses=idx), f, indent=1)
    print("done: %d witnesses in %.1f s" % (len(idx), T()))

if __name__ == '__main__':
    main(sys.argv[1:])
