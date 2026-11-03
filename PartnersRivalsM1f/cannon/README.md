# `riv1` — exhaustive verification of the rival closed form over the memory-one continuum

Verifies, over the full four-parameter space p = (p1,p2,p3,p4) in [0,1]^4, the eps -> 0
characterisation of the RIVAL sets

    Riv+  (T > S)  =  { p4 = 0,  p2 + p3 <= 1 }  minus the edge { p2 = 1, p3 = 0 }
    Riv-  (T < S)  =  { p1 = 1,  p2 + p3 >= 1 }  minus the same edge

against direct computation of the exact eps -> 0 limit of the pair chain.  See
`../Text/rivals.md` for the theorem, its proof, and its attribution — the closed form itself
is Hilbe, Traulsen & Sigmund (GEB 2015) Prop. 2/3 in the error-free model; what is verified
here is the error-limit version, which differs on exactly one edge.

## The engine

`m1f_pairs.c` is `../../PartnersRivals/check/pairs.c` generalised from a binary resident to a
continuous one.  The leading-order (coefficient, exponent) GTH reduction is unchanged; the only
new ingredient is the alphabet of leading terms,

    p_i = 0   ->  C branch (1, eps^1),  D branch (1, eps^0)
    p_i = 1   ->  C branch (1, eps^0),  D branch (1, eps^1)
    0<p_i<1   ->  C branch (p_i, eps^0), D branch (1-p_i, eps^0)

which is the eps-clamp p^eps = eps + (1-2 eps) p read to leading order.  GTH never subtracts,
so every leading coefficient is exact; there is no threshold and no extrapolation.  The limit
of a margin is therefore either exactly 0 or bounded away from it, and the classification uses
1e-12 only to separate those two cases.

## Segments

| segment | mode | what it sweeps | co-players |
|---|---|---|---|
| 0 | `m1` | exhaustive grid of spacing 1/400 over the cube, 2.585e10 points | the 16 memory-one vertices |
| 1 | `rand` | 1e9 random points drawn from all 81 strata (each coordinate 0, 1, or uniform), off-lattice | the 16 memory-one vertices |
| 2 | `m2` | grid of spacing 1/20, 194481 points | **all 65536 memory-two** |
| 3 | `rand2` | 1e6 random points from all 81 strata | **all 65536 memory-two** |

Segments 2 and 3 test the quantifier: they check that no memory-two co-player outperforms a
memory-one resident more than the best memory-one co-player does.  That is a corollary of the
payoff-completeness theorem of LaPorte, Hilbe, Glynatsi & Nowak (PNAS 2026); these segments
confirm it numerically over the continuum, as `../../PartnersRivalsM1` did over the 16.

Each segment is cut into `SHARDS` = 128 shards, so the array is 0-511 with
`segment = id / 128`, `shard = id % 128`.

## Running

```bash
./compile.sh
sbatch --parsable riv1_probe.sbatch          # measure the true per-unit rate first
sbatch --parsable riv1.sbatch > JOBID        # then the full array
```

Record the array job id in `JOBID` and append one line to `~/jobs.tsv` in the same step.
Completion is established from `sacct -j $(cat JOBID) --format=JobID,State,ExitCode -X`, never
from an empty `squeue`.

## What a clean run looks like

Every task prints one `SUMMARY` line with `bad_plus=0 bad_minus=0 largest_violation=0` and no
`DISAGREE` lines.  Any `DISAGREE` line is a point where the closed form and the computed limit
part company and must be investigated, not averaged away.

## Harvested and removed from Cannon, 2026-09-26

The Cannon working directory of this kit was `~/PartnersRivalsM1f/riv1` (from `sacct --format=WorkDir` of its own
job ids). Its 1548 output files (`out/`: riv1 45636422, fr1 45658567, and the two probes) had never been copied back;
they are now in `out/` here, md5-identical to Cannon's. The Linux binary that ran the campaigns is kept as
`m1f_pairs.cannon-linux`. The Cannon directory was then removed (Cannon home-directory clean-up requested by Martin
via the General session).
