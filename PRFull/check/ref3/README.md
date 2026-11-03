# check/ref3: checks behind the Ref3 revisions of 2026-10-24

Each script is exact (sympy over Q, or Fractions) unless noted; each runs in under a minute and a half on 3 processes of the laptop (times below). The scripts
with a process pool (`P2/p2_numcheck.py`, `stable_sets/certify.py`, `normlead.py`, `audit.py`, `audit3.py`,
`eq4_example/offer/screen_all.py`) start `CERT_MAX_THREADS` processes when it is set (otherwise 12, or 10 for
`screen_all.py`). The printed output of each script is deposited beside it as `<script>.log` (except as noted).
`PRFull/release/reproduce_certificates.sh` reruns them (steps `m1_*`, `offer_screen`, `offer`; the scripts of
`eq4_example/` itself and `farkas_display/` are not among them) and compares every file and log they write with the
deposited one.

## eq3_boundary: the reduction to Eq. (3) on the 13 boundary strata (SI "Partners")

Claim: on {p_CC = 1, p_DD > 0} off the segments (p_CD, p_DC) = (0,1), (1,0), for every game (u,v), p = (1,a,b,c) is
stable against the 16 binary deviators exactly when I = cv - (1-a) <= 0 and II = (1-a)(u-1) + bv <= 0. On the open
facet this is the display in the SI. These scripts prove it on the 13 lower-dimensional strata as well (a, b in {0, 1,
open}, c in {1, open}), and re-derive the open-facet display.

- `P1/prove_strata.py` (5 s): the eps -> 0 limit by Markov-chain-tree cofactors in Q[eps,a,b,c]. It checks that the
  lower eps-orders of every numerator vanish identically and certifies the leading coefficient of the normaliser
  positive on the whole open stratum, so the limit holds at every point, not only generically. It then writes each gain
  as l1 I + l2 II - mu and certifies l1, l2 >= 0 and mu = 0 by factoring into irreducibles and checking each factor
  (manifestly positive, or all coefficients of one sign after x = t/(1+t)). `P1/crosscheck.py` compares the limits
  with `PartnersRivalsM1f/check/exact_partners.py` (832 exact comparisons). Output: `prove_strata.log`, `.json`;
  `crosscheck.log`.
- `P2/run_all.sh` (47 s on 12 processes, 81 s on 3): an independent proof by regular perturbation of the eps = 0
  chain (closed classes, absorption probabilities, first-order generator), with no cofactors; exact numerical limits at
  245 points near the stratum edges; and a root search for zeros of the 106 factors. Prints "ALL CERTIFICATES PASSED".
  Output: `p2_summary.json`, `p2_proof_report.txt`, `p2_polys.pkl`; `run_all.log` is the printed output of one run,
  whose lists of sympy factors come in an order that depends on Python's hash seed.
- `referee/referee_check.py` (7 s): P1 and P2 limits agree on all 832 entries (13 strata x 16 deviators x 4 states);
  1568 pointwise limits by exact elimination at eps = 1e-60 with extrapolation; 38 514 (u,v) checks of the
  equivalence, including points on and 1e-9 off I = 0 and II = 0. 0 problems (`referee_check.log`).

Result: proved on all 13 strata; mu = 0 everywhere; ALLD and 0111 gain I/(1-a+c) and II/(1-a+2b) on every stratum;
where a = 1, W_CD = 0 for all 16 deviators and both conditions read v <= 0. The two excluded segments are exactly
where the self-play normaliser lead (1+b-a)(1-a-b+2ab) vanishes.

Added for Ref4 (2026-10-24):
- `box/mk_strata_table.py` (3 s) writes SI Box 3 from `P1/prove_strata.json` and asserts its legend. It always writes the
  box without its title and legend (`box/si_box3_strata_body.tex`); the title and legend are part of the paper
  (`Text/si_box3_legend.tex`, not in release 1.0.0), and where they are present it also writes the whole box
  (`box/si_box3_strata.tex`, copied to `Text/`). `box/check_box_vs_exact_partners.py si_box3_strata.tex` (or
  `si_box3_strata_body.tex`) reads the typeset box back against exact_partners.py (204 listed gains, 1044 display gains,
  78 L > 0 checks, 0 mismatches). Logs: `mk_strata_table.log` (with the legend present), `check_box.log`.
- `stable_dim/facets_cert.py` (9 s): gains on the open facet p_DD = 0 and on the alternating face, certified; on the six
  other open facets every outcome recurs in self-play. `stable_dim/zd_check.py`: the zero-determinant plane of the SI
  (eq. si-zd), 422 points tie for all 16 deviators, 200 nearby points do not. `stable_dim/f2_display_check.py`: the
  p_DD = 0 display as typeset, 960 triples. `stable_dim/zd_derivation.py`: the plane from the optimality equations.

## stable_sets: the other memory-one stability statements (Ref4 audit)

`certify.py` (6 s on 12 processes, 12 s on 3): symbolic certificates, stratum by stratum, for the facet p_DD = 0 off its
edges, the alternators, the frozen face (stable only on T = S) and (1,0,1,y) on the switch line; 24 jobs, one left open
(D1, below). `certify_edges.py`: Eq. (3) on the two efficient edges with p_DD = 0, (1,a,1,0) and (1,1,b,0).
`d1_switchline.py`: (1,1,0,y) on the switch line is stable only at (1,0) (every gain a multiple of u - 1; ALLD gains
1 - u, ALLC (u-1)/(y+2)). `audit.py`, `audit3.py`, `normlead.py` and `audit_summary.txt`: the exact sweep over all 81
strata. Logs: `certify.log` (with a FAIL line for D1, the job it leaves open), `certify_edges.log`, `d1_switchline.log`,
`normlead.log`, `audit_summary.txt` (the output of `audit.py`, which also writes `audit_out.json`, deposited as well) and
`audit3.log`.

## cycle_family: the strategies (0, p_CD, p_DC, 1) (Ref4 item 1)

`prove.py` (3 s, run from this directory): certified gains on all nine strata. Off the corners (0,0,1,1) and (0,1,0,1)
stable exactly where p_CD <= -2v and p_DC <= 1 - 2u; (0,0,1,1) only at (0,0), (0,1,0,1) only at (1/2,-1/2).
`crosscheck.py` (70 s): 405,769 exact tests against exact_partners.py, 0 disagreements for the corrected statement,
1212 for the old one, all at the two corners. `counterexample.py`: the evaluation's example at finite eps.

## farkas_display: the 16 gains as typeset

`python3 farkas_from_tex.py . ../../../Text/ms.tex` parses the gain formulas from the SI sentence in ms.tex and compares
them with the exact limits (`m1rows.py`) at 8 random points of the open facet: 128 gains, 0 mismatches
(`farkas_from_tex.log`). It reads the paper's source, `Text/ms.tex`, which release 1.0.0 leaves out, so it runs only
in release 1.1.0 (with the paper).

## eq4_example: the worked example of the Eq. (4) strategy (main text, after Eq. (4))

- `workedex.py`: self-play after one error (each player exploited 8 times on average counting the error, exploitations
  balanced on every path); the competitiveness threshold in p_<DC,CD> with the other coordinates fixed (bisection
  against all 65536 binary co-players at eps = 1e-10: 0.3846157, against 5/13 = 0.3846154).
- `exact.py`, `renewal_check.py`: exact limits over Q(x, eps). The eight co-players 001*00011011*0*1 bind at 5/13; all
  defect against the offer; their margin is (5-13x)/(2(25-13x)). Renewal at the offer state: with an offer the
  co-player nets +1 exploitation over 3 rounds, without one it nets -5/8 over 25/4 rounds on average. At x = 1/4 their
  margin is 7/174 > 0, so they are not the strategy's worst co-players there.
- `screen.py`, `screen2.py`, `worst_at_quarter.py` (with `we_core.py`): batched screens of all 65536 co-players.
- The printed output of each script is in `<script>.log`. `workedex.py` and the three screens compute in floating
  point (numpy), and the last line of `screen.log` is its run time.
- `offer/` (Ref4 item 3): `screen_all.py N D` writes the exact limit margins of all 65536 co-players at x = N/D (about
  22 s on 3 processes; outputs `margins_*.txt` for 1/100, 1/4, 2/5, 5/13, and `screen_all.log`, the output of the four
  runs); `interval.py` reconstructs each margin as a ratio of affine functions of x (structural lemma), checks it at
  5/13 and 7/10, and reads off the range: no numerator negative at 0, smallest root 5/13, shared by 400 co-players, so
  the point is a friendly rival exactly for 0 < x <= 5/13. `renewal.py`, `lengths.py`, `binders.py`, `two_class.py`: all
  400 defect against the offer and net +1 per offer and -5/8 per non-offer; the eight 001*00011011*0*1 have the shortest
  excursions (3 and 25/4 rounds). `selfplay.py`: recovery after one error (2/x exploitations each, 8 at x = 1/4). Each
  script's output is in `<script>.log`.
