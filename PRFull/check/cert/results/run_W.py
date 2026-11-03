"""Referee P2 (a): positive-volume certificate in wedge W (T>S, below the switch line).
Stratum 1*1*001*1*1**0*0 (dimension 7), box = SI Table 2 witness +- 1/16, vertex route, every
unresolved vertex margin computed exactly.  Run: cd PRFull/check && nice -n 10 python3 cert/results/run_W.py"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p2lib import *
if __name__ == '__main__':
    pat = '1*1*001*1*1**0*0'
    c = [1, F(1,2), 1, F(1,2), 0, 0, 1, F(1,2), 1, F(1,4), 1, F(1,4), F(1,4), 0, F(1,4), 0]
    h = F(1, 16); fcs = en.free_coords(pat)
    out = box_certificate('W_box', 'W', pat, c, [c[j] - h for j in fcs], [c[j] + h for j in fcs], exact_all=True)
    st = ("STATEMENT.  Every p in the closed 7-dimensional box below (pinned coordinates as in the stratum) is a\n"
          "friendly rival of wedge W in the eps -> 0 limit: its self-play limit is (1,0,0,0), and against every\n"
          "binary memory-two co-player (hence, by the box-vertex lemma, every memory-two co-player)\n"
          "w_DC >= w_CD.  The box is the SI Table 2 witness +- 1/16 in every free coordinate; the witness is\n"
          "therefore an interior point of the friendly-rival set.")
    print(write_results(out, os.path.join(HERE, 'W_box'), st))
