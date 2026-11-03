#! /bin/bash
# cw / kappa: the mutation kernel's mass on each atom (E,N,R) and each of the 81 strata, per game -- the
# neutral reference of Fig 4(c)'s enrichment, share / kappa.  kappa.x = kappa.c + m1cls.c.
#
#     one task, looping over the 512 cells isl = 3, 6, ..., 1536;  M = 1e7 kernel draws per cell
#     RNG: kappa.x's own MT19937, init_genrand(41000000 + isl), independent of the dynamics
#
# ONE task, not 512: a cell takes ~2 s (2.07 s locally at M = 1e7), so 512 array elements would spend more
# time in the scheduler than computing.  Idempotent per cell: kappa/k<isl> is written through a temporary
# file and a rename, and a cell whose file has its single 98-field line is skipped, so a rerun after a
# kill finishes the rest.  Local estimate 512 x 2.1 s = 18 min; -t 120.
#SBATCH -p sapphire
#SBATCH -J PartnersRivalsM1f.cwkappa
#SBATCH --comment=PartnersRivalsM1f
#SBATCH -t 120
#SBATCH -c 1
#SBATCH --mem=100
#SBATCH -o log/o_kappa.txt
#SBATCH -e log/e_kappa.txt

module load intel/25.2.1-fasrc01
M=10000000
cd kappa || { echo "kappa/ missing -- use submit.sh, which stages games.dat"; exit 1; }
[ -f games.dat ] || { echo "kappa/games.dat missing -- use submit.sh"; exit 1; }
nbad=0; nrun=0; nskip=0
for isl in $(seq 3 3 1536); do
  if [ -f "k$isl" ] && [ "$(wc -l < "k$isl")" -eq 1 ] && [ "$(awk '{print NF}' "k$isl")" -eq 98 ]; then
    nskip=$((nskip+1)); continue
  fi
  if ../kappa.x "$isl" "$M" > ".k$isl.new" && [ "$(awk '{print NF}' ".k$isl.new")" -eq 98 ]; then
    mv -f ".k$isl.new" "k$isl"; nrun=$((nrun+1))
  else
    echo "kappa.x failed at isl $isl" >&2; rm -f ".k$isl.new"; nbad=$((nbad+1))
  fi
done
echo "kappa: $nrun computed, $nskip already present, $nbad failed"
[ "$nbad" -eq 0 ]
