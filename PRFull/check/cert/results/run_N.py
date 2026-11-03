"""Referee P2 (b): positive-volume certificate in wedge N (T>S, above the switch line).
Stratum *10***00*1***010 (dimension 8), box = SI Table 2 witness +- 1/12, vertex route, every
unresolved vertex margin computed exactly (4.5M mpq evaluations, about 9 min on 3 processes).
Run: cd PRFull/check && nice -n 10 python3 cert/results/run_N.py"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p2lib import *
PAT = '*10***00*1***010'
CEN = [F(2,3), 1, 0, F(1,3), F(1,3), F(1,3), 0, 0, F(2,3), 1, F(2,3), F(1,3), F(1,3), 0, 1, 0]
H = F(1, 12)
if __name__ == '__main__':
    exact_all = '--corner-identities' not in sys.argv
    fcs = en.free_coords(PAT)
    out = box_certificate('N_box', 'N', PAT, CEN, [CEN[j] - H for j in fcs], [CEN[j] + H for j in fcs],
                          exact_all=exact_all, max_exact_seconds=840)
    st = ("STATEMENT.  Every p in the closed 8-dimensional box below (pinned coordinates as in the stratum) is a\n"
          "friendly rival of wedge N in the eps -> 0 limit: its self-play limit is (0,1/2,1/2,0), and against\n"
          "every binary memory-two co-player (hence every memory-two co-player) w_DC >= w_CD.  The box is the\n"
          "SI Table 2 witness +- 1/12 in every free coordinate; the witness is an interior point of the\n"
          "friendly-rival set.  (Trial boxes: +-1/9, +-1/8, +-1/6 fail, a certified-negative vertex against co-player 28672.)")
    print(write_results(out, os.path.join(HERE, 'N_box'), st))
