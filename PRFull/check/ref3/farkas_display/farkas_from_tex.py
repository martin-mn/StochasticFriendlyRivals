import re, random, sys
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application
sys.path.insert(0, sys.argv[1])
from m1rows import stationary_limit
s = open(sys.argv[2]).read()
seg = s[s.index(r'the gain $\pi(\sigma,p)-\pi(p,p)$'):s.index('Every denominator is positive')]
seg = seg.split('is ',1)[1]
names = {r'\WSLS': '1001', r'\TFT': '1010', r'\Grim': '1000', r'\ALLC': '1111', r'\ALLD': '0000'}
a, b, c, u, v, I, J = sp.symbols('a b c u v I J')
loc = {'a': a, 'b': b, 'c': c, 'I': I, 'J': J}
T = standard_transformations + (implicit_multiplication_application,)
claims = []
for part in re.split(r';\s*(?:and\s+)?', seg.strip().rstrip('.')):
    m = re.match(r'\$(.*?)\$ for (.*)', part.strip())
    expr, who = m.group(1), m.group(2)
    expr = expr.replace(r'I\!I', 'J').replace('[', '(').replace(']', ')')
    if expr.startswith('0'): f = sp.Integer(0)
    else:
        num, den = re.match(r'(.*)/(\(.*\)|\d)$', expr).groups() if not expr.startswith('(') or ')/' in expr else (None, None)
        f = parse_expr(expr.replace('/', ')/(').join(['(', ')']), local_dict=loc, transformations=T)
    for w in re.split(r',\s*|\s+and\s+', who):
        w = w.strip().replace('{}', '')
        if w: claims.append((names.get(w, w.strip('$')), f))
print(len(claims), 'deviators:', sorted(x for x, _ in claims))
assert len({x for x, _ in claims}) == 16
Iex = c*v - (1 - a); Jex = (1 - a)*(u - 1) + b*v
random.seed(7); bad = 0
for t in range(8):
    pt = {a: sp.Rational(random.randint(1, 99), 100), b: sp.Rational(random.randint(1, 99), 100), c: sp.Rational(random.randint(1, 99), 100)}
    p = (sp.Integer(1), pt[a], pt[b], pt[c])
    assert stationary_limit(p, p) == [1, 0, 0, 0]
    for bits, f in claims:
        W = stationary_limit(tuple(sp.Integer(int(x)) for x in bits), p)
        gain = -W[3] + (u - 1)*W[1] + v*W[2]
        pred = f.subs({I: Iex, J: Jex}).subs(pt)
        if sp.simplify(gain - pred) != 0:
            bad += 1; print('MISMATCH', bits, pt, sp.factor(gain), sp.factor(pred))
print('checked', 8*len(claims), 'gains; mismatches', bad)
