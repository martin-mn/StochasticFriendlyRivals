# Task-2 verification kit (PartnersRivalsM2f)

Independent re-verification of the C/D-mirror detection-loss claim of `Text/friendly.md` §5.6.
Nothing here reuses `cannon/census.c`'s driver; only `m2f.c` (the validated engine) and
`witnesses.h` (the sound screen) are shared.

* `t2.c ncensus <sign> <lo> <hi> <K1> <K2> <seed>` — exhaustive re-census of wedge N (sign +1,
  candidates d15 = 0, i.e. g < 3^15, by the *proved* Riv+ ⊆ {p15 = 0}) or wedge E (sign −1,
  candidates d0 = 2 by the mirror image of that condition).  K1 = 64 reproduces the census's
  sampling depth; any stratum with 0 hits in K1 gets K2 more draws (early stop at 4 hits).
  This measures the K = 64 detection loss directly, with no extrapolation.
* `t2.c list <file> <shard> <nshard> <seed>` — heavy per-stratum sampling of a stratified
  sample of census strata, for the Horvitz–Thompson correction.
* `t2.c seedtest <g> <sign> <s0> <ns> <K>` — does census.c's per-stratum seeding
  `rs = g*6364136223846793005 + 1442695040888963407` produce K draws that behave like K iid
  draws?  Histogram of hits over many seeds vs Binomial(K, f).

Build: `./compile.sh`.  Submit with `--parsable`, record the id in `JOBID`, append to `~/jobs.tsv`.
