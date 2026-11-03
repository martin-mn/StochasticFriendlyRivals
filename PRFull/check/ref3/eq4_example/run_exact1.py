import sympy as sp, time
from exact import *
p = list(EX); p[9] = x
fam = [52620, 36236, 56716, 40332, 52612, 36228, 56708, 40324]
t=time.time()
for c in fam:
    m, w = margin_limit(p, code2q(c))
    print(c, ''.join(map(str,code2q(c))), 'q6=',code2q(c)[6], ' margin(x)=', m, ' w=', w)
print('time', time.time()-t)
target = (5-13*x)/(2*(25-13*x))
print('matches (5-13x)/(2(25-13x)):', sp.simplify(m-target)==0)
