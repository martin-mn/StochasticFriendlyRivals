Scripts of the adversarial review of 2026-09-27 (outputs in ../../logs/review/).  Local only; not part of the Cannon job.
  indep.py lib|full <sid>...   independent recomputation of a probe strategy's atom with the cert library (not audit.py):
                               E from exact self-play W, C by ce.certify_point, S in the q-focal chain (i_eval(q,q,p) +
                               exact lead_exact(q,p)); 'full' adds exact rationals against all 65536 q in both chains
  compare.py                   indep results against the Cannon probe records
  repcheck.py                  dumped atoms per replicate against the REP lines of the paper's runs (atom bit order)
  bounds.py                    pi_i recomputed, the draw reproduced, the 95% bounds, the W friendly rivals of the sample
  rerun.py <sid>...            current src/audit.py on probe records against Cannon's output
