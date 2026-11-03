#!/bin/bash
# Submit one cw array, doing in ONE step what the global Cannon protocol requires: stage games.dat,
# capture the array id with sbatch --parsable into JOBID_<tag>, write the array's own index set to
# IDX_<tag>, and append one line to the cross-project ledger ~/jobs.tsv
#     date  PartnersRivalsM1f  cw  jobid  ntasks  note          (tab-separated)
#
#   bash submit.sh probe     5 tasks, the rate probe   (cwp.x, itend 1e6, partition test)  -> probe/
#   bash submit.sh prod    512 tasks, the production   (cw.x, --array=3-1536:3, sapphire)  -> run1/
#   bash submit.sh kappa     1 task,  the kernel masses (kappa.x, 512 cells in one loop)    -> kappa/
#   run F2 (N=100 beta=3 u=1e-4 itend=1e8; cw2.f = mkcw2.py(cw.f); kappa is shared with F1):
#   bash submit.sh probe2    5 tasks, the F2 rate probe (cwp2.x, itend 1e7, partition test) -> probe2/
#   bash submit.sh prod2   512 tasks, F2 production     (cw2.x, --array=3-1536:3, sapphire) -> run2/
#   F2 writes JOBID_probe2 / JOBID_prod2 only: the plain JOBID file is F1's campaign id (48366968) and is
#   never touched by an F2 submission.
#
# Order: probe (and kappa, independent) -> bash verify.sh probe -> set prod.slurm's -t from the probe's
# rate -> prod.  Build first on a login node: bash compile.sh; bash compile.sh probe.
#
# games.dat is staged by write-to-temp-then-rename (cp -f truncates its target first and a task starting
# in that instant would read a short list; mv on one filesystem is atomic) -- mw's submit.sh, kept.
set -e
cd "$(dirname "$0")"
case "$1" in
  probe) SCRIPT=probe.slurm; RD=probe; EXE=cwp.x;   NT=5;   IDX="3 1527 1530 1533 1536"
         NOTE="[probe] continuous-M1 WF rate probe: sf1 F1 process (N=1000 beta=100 u=1e-2 eps=1e-4) at itend=1e6 (1/10), 5 cells = mw probe1's" ;;
  prod)  SCRIPT=prod.slurm;  RD=run1;  EXE=cw.x;    NT=512; IDX="$(seq 3 3 1536)"
         NOTE="[prod] continuous-M1 WF (exact-alpha kernel, exact eps->0 atoms): sf1 F1 process N=1000 beta=100 u=1e-2 itend=1e7 eps=1e-4, 512 cells x 10 reps -- Fig4c/5c" ;;
  kappa) SCRIPT=kappa.slurm; RD=kappa; EXE=kappa.x; NT=1;   IDX="1"
         NOTE="[kappa] kernel atom/stratum masses, M=1e7 per cell, 512 cells in one task -- Fig4c enrichment reference" ;;
  probe2) SCRIPT=probe2.slurm; RD=probe2; EXE=cwp2.x; NT=5; IDX="3 1527 1530 1533 1536"
         NOTE="[probe2] continuous-M1 WF run F2 rate probe: sf2 F2 process (N=100 beta=3 u=1e-4 itend=1e8 eps=1e-4) at itend=1e7 (1/10), 5 cells = mw probe2's, seeds 42000000+" ;;
  prod2) SCRIPT=prod2.slurm; RD=run2;   EXE=cw2.x;  NT=512; IDX="$(seq 3 3 1536)"
         NOTE="[prod2] continuous-M1 WF run F2 (exact-alpha kernel, exact eps->0 atoms): sf2 F2 process N=100 beta=3 u=1e-4 itend=1e8 eps=1e-4, 512 cells x 10 reps, seeds 42000000+ -- SI Fig4c/5c (N=100)" ;;
  *) echo "usage: bash submit.sh {probe|prod|kappa|probe2|prod2}"; exit 1 ;;
esac
TAG=$1

[ -f "$SCRIPT" ] || { echo "$SCRIPT not found"; exit 1; }
[ -x "$EXE" ] || { echo "$EXE not built -- run: bash compile.sh (and bash compile.sh probe)"; exit 1; }
# never a stub build: the real classifier exports m1cls_detail, the quarantined stub does not
nm "$EXE" 2>/dev/null | grep -q m1cls_detail || { echo "REFUSING: $EXE is not linked with the real m1cls.c"; exit 1; }
[ "$TAG" = kappa ] || [ -x task.sh ] || { echo "task.sh missing"; exit 1; }
# F2: the binary must be an F2 build of cw2.f (its w-file header names maxc 202 = 2*100+2; cw.x says 2002),
# and cw2.f must still be exactly mkcw2.py's output of cw.f
case "$TAG" in probe2|prod2)
  # ifx stores the w-file header's literals interleaved ('...high water (m' ... 'axc 202); ...'), so match the
  # fragment 'axc 202)', which an F2 build has under gfortran and ifx and an F1 build ('axc 2002)') never has
  { strings "$EXE" | grep -q 'axc 202)' && ! strings "$EXE" | grep -q 'axc 2002)'; } \
    || { echo "REFUSING: $EXE is not a cw2.f (N = 100) build"; exit 1; }
  python3 -B mkcw2.py --check >/dev/null || { echo "REFUSING: cw2.f is not mkcw2.py's output of cw.f"; exit 1; } ;;
esac

# no second copy of a live array of the same tag: the two would write the same files concurrently
if [ -f "JOBID_$TAG" ]; then
  live=$(squeue -j "$(cat JOBID_$TAG)" -h -o '%i' 2>/dev/null | wc -l | tr -d ' ')
  if [ "${live:-0}" -gt 0 ]; then
    echo "REFUSING: $live task(s) of $TAG array $(cat JOBID_$TAG) are still queued or running."; exit 1
  fi
  # an empty squeue means NOT RUNNING and nothing more (global protocol): a finished array vanishes from
  # it whether every task succeeded or every task died.  Resubmitting on that basis is forbidden.
  if [ "${RESUBMIT:-0}" != 1 ]; then
    echo "REFUSING: JOBID_$TAG exists ($(cat JOBID_$TAG)) and is not in the queue.  Decide what happened from"
    echo "sacct per task first:  bash verify.sh $TAG   -- then, only if a resubmission is really wanted,"
    echo "RESUBMIT=1 bash submit.sh $TAG  (task.sh skips cells that are already complete)."; exit 1
  fi
  echo "note: RESUBMIT=1: JOBID_$TAG ($(cat JOBID_$TAG)) is kept as JOBID_$TAG.prev"
  cp -f "JOBID_$TAG" "JOBID_$TAG.prev"
fi

# the ledger's task count must be the script's own, not a number typed twice
spec=$(awk -F= '/^#SBATCH --array=/{print $2}' "$SCRIPT")
n=$(echo "$IDX" | wc -w | tr -d ' ')
[ "$n" -eq "$NT" ] || { echo "internal: $TAG index set has $n entries, NT says $NT"; exit 1; }
if [ -n "$spec" ]; then
  nspec=$(python3 - "$spec" <<'PY'
import sys
t = 0
for part in sys.argv[1].split(','):
    if '-' in part:
        rng, _, st = part.partition(':')
        a, b = rng.split('-')
        t += len(range(int(a), int(b) + 1, int(st or 1)))
    else:
        t += 1
print(t)
PY
)
  [ "$nspec" -eq "$NT" ] || { echo "internal: --array=$spec is $nspec tasks, NT says $NT"; exit 1; }
fi

# courtesy between projects (global protocol): before a large array, report what is already committed
if [ "$NT" -gt 100 ]; then
  inq=$(squeue -u "$USER" -h -t PENDING,RUNNING -o '%A' | wc -l | tr -d ' ')
  echo "the user already has $inq job(s)/task(s) pending or running (all projects)"
  if [ "${QUEUE_OK:-0}" != 1 ]; then
    echo "REFUSING without QUEUE_OK=1: read the ledger breakdown of that number, report it, and ask Martin"
    echo "before adding a $NT-task array to it (global protocol, courtesy between projects)."; exit 1
  fi
fi
echo "$SCRIPT: --array=${spec:-none}, $NT task(s), $RD/, $EXE"

mkdir -p "$RD" log
cp -f games.dat "$RD/.games.dat.new" && mv -f "$RD/.games.dat.new" "$RD/games.dat"
echo "$IDX" | tr ' ' '\n' | grep -v '^$' > "IDX_$TAG"

jid=$(sbatch --parsable "$SCRIPT")
jid=${jid%%;*}
echo "$jid" > "JOBID_$TAG"
[ "$TAG" = prod ] && echo "$jid" > JOBID     # the campaign's id under the protocol's own name
printf '%s\tPartnersRivalsM1f\tcw\t%s\t%s\t%s\n' "$(date +%F)" "$jid" "$NT" "$NOTE" >> ~/jobs.tsv
echo "submitted $SCRIPT as $jid ($NT task(s)); JOBID_$TAG and IDX_$TAG written, ~/jobs.tsv appended"
