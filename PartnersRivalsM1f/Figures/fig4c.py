#!/usr/bin/env python3
"""
fig4c.py -- Figure 4 of ../../PartnersRivals1 with the CONTINUOUS memory-one run added.

    python3 fig4c.py               ->  Fig4c.pdf / .png            (needs ../cannon/cw/packed/{cw1_e4,kappa_e4}.tsv)
    python3 fig4c.py --synthetic   ->  tests/Fig4c_synthetic.pdf / .png  (fabricated stand-in, watermarked SYNTHETIC)
    python3 fig4c.py --short       ->  tests/Fig4c_localshort.pdf / .png (a SHORT LOCAL RUN of the real chain,
                                       ../cannon/cw/packed_local, watermarked SHORT LOCAL RUN)
    python3 fig4c.py --rows m1b,m1c,m2b   (the default; any subset/order of the keys of ROWSPEC)
    python3 fig4c.py --run F2 [--synthetic | --short]
                                   ->  SIFig4c, tests/SIFig4c_{synthetic,localshort}.pdf / .png: SI Figure 4, the N = 100 run F2
                                       (N = 100, beta = 3, mu = 1e-4, eps = 1e-4, 10 x 1e8 generations): rows
                                       wfdata.load("m1", "F2"), cwdata.load(run="F2") (packed/cw2_e4.tsv, short mode
                                       packed_local2/), wfdata.load("m2", "F2").  --run F1 is the default and is
                                       exactly the invocation it always was.

Rows, top to bottom (all N = 1000, beta = 100, mu = 1e-2, eps = 1e-4, 10 x 1e7 generations, second half sampled):
    m1b  Memory-1, binary       the 16 binary memory-one strategies   (DiskM1WF run F1, wfdata.load("m1", "F1"))
    m1c  Memory-1, continuous   the whole cube [0,1]^4 under the exact-alpha kernel (../cannon/cw, cwdata.load())
    m2b  Memory-2, binary       the 65536 binary memory-two strategies (DiskM2WF run F1, wfdata.load("m2", "F1"))
Columns and colour keys are the original's exactly -- the drawing is fig4.draw itself, imported from FinalFigures
(read-only; nothing there is modified): (a) efficiency = mean per-round payoff / Emax, viridis 0..1, black below 0;
(b) the atom with the largest share; (c) the most enriched atom, argmax_a share(a) / kappa(a), where kappa(a) is the
share of the strategy space in atom a at that game: |a|/16 or |a|/65536 for the binary rows (fig4.winners: the
argmax of S/|a|), the mutation kernel's mass on a (kappa_e4.tsv, M Monte Carlo draws per game) for the continuous
row.  Atoms with kappa = 0 are excluded, as the original excludes |a| = 0.  Column (c) has a categorical key, so no
enrichment colour range enters the figure; the enrichment VALUES of the continuous row are printed, and are not
bounded by 16 as the binary-M1 ones are but by 1/kappa(a) (see cwdata.enrichment for the kappa_hat = 0 rule).
"""
import argparse
import os
import sys

import numpy as np

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cwdata                                                      # noqa: E402  puts FinalFigures on sys.path
import wfdata                                                      # noqa: E402
import fig4                                                        # noqa: E402  the drawing (read-only import)

ORDER = wfdata.ORDER
# the originals label a row 'Memory-1, $N = 1000$'; the continuum sheets also say which strategy space
NLAB = {"F1": "$N = 1000$", "F2": "$N = 100$"}
ROWSPECS = {run: {"m1b": ("Memory-1 binary, " + NLAB[run], ("m1", run)),
                  "m1c": ("Memory-1 continuous, " + NLAB[run], None),
                  "m2b": ("Memory-2 binary, " + NLAB[run], ("m2", run))} for run in ("F1", "F2")}
ROWSPEC = ROWSPECS["F1"]
DEFAULT = "m1b,m1c,m2b"
STEMS = {"F1": "Fig4c", "F2": "SIFig4c"}


def continuous(synthetic, short=False, run="F1"):
    """(E, wa, we) of the continuous run, with the printed summary fig4.winners gives the binary ones"""
    uv_ref = wfdata.load("m1", run)[3]
    d = cwdata.load(synthetic, uv_ref=uv_ref, short=short, run=run)
    S, K, E, M = d["S"], d["K"], d["E"], d["M"]
    wa = np.where((K > 0) | (S > 0), S, -np.inf).argmax(1)
    enr, nfl = cwdata.enrichment(S, K, M)
    we = enr.argmax(1)
    tag = "m1 continuous%s%s" % ("" if run == "F1" else " " + run, " [SYNTHETIC]" if synthetic else " [SHORT LOCAL RUN]" if short else "")
    print("%s (N = %d, beta = %d, mu = %s): efficiency %.3f .. %.3f, median %.3f, negative at %d games, below 0.9 at %d"
          % (tag, d["D"]["N"], d["D"]["beta"], d["D"]["muplain"], E.min(), E.max(), np.median(E), (E < 0).sum(), (E < 0.9).sum()))
    print("   most abundant: " + ", ".join("%s %d" % (ORDER[i], (wa == i).sum()) for i in range(7) if (wa == i).any()))
    print("   most enriched: " + ", ".join("%s %d" % (ORDER[i], (we == i).sum()) for i in range(7) if (we == i).any()))
    ew = enr[np.arange(512), we]
    kw = K[np.arange(512), we]
    ceil = np.where(K > 0, 1.0 / np.where(K > 0, K, 1.0), 0.0).max(1)
    print("   enrichment of the winner: %.3g .. %.3g, median %.3g; above 16 (the binary-M1 ceiling) at %d games; "
          "ceiling max_a 1/kappa(a): %.3g .. %.3g" % (ew.min(), ew.max(), np.median(ew), (ew > 16).sum(), ceil.min(), ceil.max()))
    # how well kappa_hat resolves the winner: relative s.e. of kappa_hat is sqrt((1-k)/(k M)); the share's own
    # error (replicate spread) is not in the pack and is not included
    rse = np.sqrt((1 - kw) / np.maximum(kw * M, 1))
    e2 = np.sort(enr, 1)[:, -2]
    close = (e2 > -np.inf) & (ew / np.where(e2 > 0, e2, 1) - 1 < 3 * np.sqrt(2) * rse)
    print("   winner's kappa: min %.3g (%d kernel draws); kappa_hat < 1e-4 at %d games; winner not resolved from the "
          "runner-up by kappa_hat's own noise (3 sigma) at %d games; kappa_hat = 0 with share > 0 (floored at 3/M): %d entries"
          % (kw.min(), int(round(kw.min() * M)), (kw < 1e-4).sum(), close.sum(), nfl))
    # the enrichment floor: an atom the population only visits as fresh mutants holds a share of order mu * kappa, so
    # share / kappa bottoms out near mu (1e-2 in F1, 1e-4 in F2); and a winner held at a small share is an atom of
    # tiny kernel mass that the population enters rarely, not the atom it mostly plays
    pos = (S > 0) & (K > 0)
    sw = S[np.arange(512), we]
    print("   enrichment floor: min share/kappa over held atoms %.3g, 1st percentile %.3g (mu = %s); the winner's share "
          "is below 1e-2 at %d games, below 1e-1 at %d; the winner is not the most abundant atom at %d games"
          % ((S[pos] / K[pos]).min(), np.percentile(S[pos] / K[pos], 1), d["D"]["muplain"], (sw < 1e-2).sum(),
             (sw < 1e-1).sum(), (we != wa).sum()))
    return E, wa, we


def rows(keys, synthetic, short=False, run="F1"):
    out = []
    for k in keys:
        label, src = ROWSPECS[run][k]
        if src is None:
            out.append((label,) + continuous(synthetic, short, run))
        else:
            out.append((label,) + fig4.run(*src)[:3])
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--synthetic", action="store_true", help="draw the continuous row from fabricated data, watermarked")
    ap.add_argument("--short", action="store_true", help="draw the continuous row from the kit's SHORT LOCAL RUN "
                    "(../cannon/cw/packed_local), watermarked")
    ap.add_argument("--rows", default=DEFAULT, help="comma-separated keys of ROWSPEC, top to bottom (default %s)" % DEFAULT)
    ap.add_argument("--run", default="F1", choices=["F1", "F2"], help="F1: N = 1000, Figure 4 (default); "
                    "F2: N = 100, SI Figure 4 -> SIFig4c*")
    a = ap.parse_args(argv)
    keys = a.rows.split(",")
    assert all(k in ROWSPEC for k in keys), "unknown row key in %r; known: %s" % (a.rows, ", ".join(ROWSPEC))
    if a.synthetic and not all(os.path.exists(p) for p in cwdata.paths(True, run=a.run)):
        cwdata.make_synthetic(run=a.run)
    assert not (a.synthetic and a.short), "--synthetic and --short are exclusive"
    stem = STEMS[a.run] + ("_synthetic" if a.synthetic else "_localshort" if a.short else "")
    if a.synthetic or a.short:                                     # test renders live in tests/, not beside the results
        os.makedirs(os.path.join(HERE, "tests"), exist_ok=True)
        stem = os.path.join("tests", stem)
    fig4.HERE = HERE                                               # in memory only: fig4.draw writes into this directory
    R = rows(keys, a.synthetic, a.short, a.run)
    with cwdata.watermark(a.synthetic), cwdata.watermark(a.short, "SHORT LOCAL RUN", cwdata.short_footer(a.run) if a.short
                                                         else "", 80):
        fig4.draw(R, stem)
    print("(the file above is Figures/%s.pdf / .png of PartnersRivalsM1f, not FinalFigures)" % stem)


if __name__ == "__main__":
    main()
