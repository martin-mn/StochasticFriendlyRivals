from ref_check import occ, payoffdiff, selfpay, Emax, COP, is_eff, is_rival
from fractions import Fraction as F
print("=== CX1: delta=0, wedge S (u,v)=(2/5,-1).  p=(p0=1; 0,0,0,0)  i.e. p1=0 ===")
u,v,d=F(2,5),F(-1),F(0); p0=F(1); p=[F(0)]*4
print(" self-play payoff",selfpay(p0,p,d,u,v),"  Emax",Emax(u,v),"  efficient:",is_eff(p0,p,d,u,v))
print(" min margin over 32 det co-players:",min(payoffdiff(p0,p,q0,q,d,u,v) for q0,q in COP),
      "  unbeatable:",is_rival(p0,p,d,u,v))
print(" => FRIENDLY RIVAL with p1=0, contradicting 'efficient iff strictly nice' for delta in [0,1)")
print("    (at delta=0 the efficient set is all of {p0=1}: 5^4=625 grid pts vs 5^3=125 for delta>0)")

print("\n=== CX2: T=S boundary ray of the CLOSED wedge S, (u,v)=(1/2,-1/2), delta=9/10, Grim ===")
u,v,d=F(1,2),F(-1,2),F(9,10); p0=F(1); p=[F(1),F(0),F(0),F(0)]
print(" p2+p3 =",p[1]+p[2],"   threshold 2-1/delta =",2-F(1)/d)
print(" self-play",selfpay(p0,p,d,u,v),"Emax",Emax(u,v),"efficient:",is_eff(p0,p,d,u,v))
print(" min margin:",min(payoffdiff(p0,p,q0,q,d,u,v) for q0,q in COP),"unbeatable:",is_rival(p0,p,d,u,v))
print(" => FRIENDLY RIVAL below the stated threshold, because T=S makes rivalry vacuous.")
print("    On the closed wedge S the set is {p0=p1=1, p2+p3>=2-1/d} UNION ({u-v=1} x {p0=p1=1}).")
