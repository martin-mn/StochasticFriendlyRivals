#!/usr/bin/env python3
"""runstep.py -- run one step of a reproduce_*.sh script: time it, measure its memory, judge it, and compare what it
wrote with the deposited files.

    python3 runstep.py --name NAME --label LABEL --cwd DIR --logdir LOGS --summary SUMMARY.jsonl [--run-id ID]
                       [--scan WORKDIR=DEPOSITDIR ...] [--strip PATH ...] [--need PROGRAM] [--expect REGEX ...]
                       [--stdout-vs DEPOSITED_LOG] [--must-write FILE ...] [--ignore REGEX ...] [--normalize REGEX ...]
                       [--nocounts] [--no-fail-words] -- COMMAND ...

The command runs under "nice -n 10" in DIR, its output going to LOGS/NAME.out and LOGS/NAME.err.  The step FAILS if
the command exits non-zero, if its output contains FAIL, FAILED, DIFFERS, NOT CERTIFIED, Traceback or "rows_failing
N" / "N failed" with N > 0 (--no-fail-words, --nocounts switch these readings off where a step reports other things'
states), if a line matching an --expect pattern is missing, or if a --must-write file was not rewritten.  After it,
every file under each
WORKDIR written by the step is compared with the file at the same relative path under DEPOSITDIR (compare.py;
--ignore leaves a path out, --normalize masks a pattern in both texts, said in the verdict);
a file that differs is reported but does not fail the step by itself: the script's final verdict counts them.
Wall-clock time and the peak resident set size of the largest process the step started (getrusage of the
children) are recorded.  One JSON line per step is appended to SUMMARY.jsonl.  Exit status: 0 passed, 1 failed,
3 skipped (a required program or deposited input is missing).
"""
import argparse
import json
import os
import re
import resource
import shutil
import subprocess
import sys
import time

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare                                                    # noqa: E402

FAIL_RX = re.compile(r"(^|[^A-Za-z_])(FAIL|FAILED|DIFFERS|NOT CERTIFIED|NOT_CERTIFIED|Traceback)([^A-Za-z_]|$)")
COUNT_RX = re.compile(r"(rows_failing|violations?|mismatch(es)?|failed|failing|NOT_EXCLUDED) ?[:=]? ?([1-9]\d*)\b"
                      r"|\b([1-9]\d*) (failed|mismatch(es)?|violations?)\b", re.I)


def rss_bytes(ru):
    return ru * (1 if sys.platform == "darwin" else 1024)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--label", default="")
    ap.add_argument("--cwd", required=True)
    ap.add_argument("--logdir", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--run-id", default="", help="the invocation of the script this step belongs to, kept in its record "
                    "(a summary that keeps the records of earlier --only runs tells them apart by it)")
    ap.add_argument("--scan", action="append", default=[])
    ap.add_argument("--strip", action="append", default=[])
    ap.add_argument("--need", action="append", default=[])
    ap.add_argument("--need-file", action="append", default=[], help="a deposited input the step needs (relative to "
                    "--cwd or absolute); if one is missing the step is skipped and says so (see --need-what)")
    ap.add_argument("--need-what", default="", help="where the missing inputs are found, for the skip message")
    ap.add_argument("--expect", action="append", default=[])
    ap.add_argument("--nocounts", action="store_true", help="do not read 'N failed' style counts as failures")
    ap.add_argument("--no-fail-words", action="store_true", help="do not read FAIL/FAILED/DIFFERS in the output as failures "
                    "(for a step that reports other runs' states, e.g. the Slurm accounting); --expect decides")
    ap.add_argument("--stdout-vs", action="append", default=[], help="a deposited log to compare the step's output with")
    ap.add_argument("--must-write", action="append", default=[], help="a file the step must (re)write, or it fails")
    ap.add_argument("--normalize", action="append", default=[], help="a regex replaced by a placeholder in both files "
                    "before text files are compared (e.g. an arg-min among exact ties, which depends on the process order)")
    ap.add_argument("--ignore", action="append", default=[], help="a regex of written paths (relative to WORKDIR) not to "
                    "compare, e.g. a log that every run appends to")
    ap.add_argument("cmd", nargs=argparse.REMAINDER)
    a = ap.parse_args()
    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == "--" else a.cmd
    compare.EXTRA_NORM[:] = a.normalize
    os.makedirs(a.logdir, exist_ok=True)
    rec = dict(name=a.name, label=a.label, cwd=a.cwd, command=" ".join(cmd), run=a.run_id)
    print("\n[%s] %s\n     (in %s: %s)" % (a.name, a.label, a.cwd, " ".join(cmd)), flush=True)
    for prog in a.need:
        if shutil.which(prog) is None:
            rec.update(status="skipped", why="%s is not installed" % prog)
            print("     SKIPPED: %s is not installed" % prog)
            with open(a.summary, "a") as fh:
                fh.write(json.dumps(rec) + "\n")
            return 3
    missing = [f for f in a.need_file if not os.path.exists(f if os.path.isabs(f) else os.path.join(a.cwd, f))]
    if missing:
        why = "missing %d input(s), e.g. %s%s" % (len(missing), missing[0], ("; " + a.need_what) if a.need_what else "")
        rec.update(status="skipped", why=why)
        print("     SKIPPED: " + why)
        with open(a.summary, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        return 3
    t0 = time.time()
    time.sleep(1.05)                                                # every file the step writes is then newer than t0
    t0 = time.time()
    fo = open(os.path.join(a.logdir, a.name + ".out"), "w")
    fe = open(os.path.join(a.logdir, a.name + ".err"), "w")
    st = subprocess.call(["nice", "-n", "10"] + cmd, cwd=a.cwd, stdout=fo, stderr=fe)
    fo.close()
    fe.close()
    dt = time.time() - t0
    peak = rss_bytes(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    out = open(os.path.join(a.logdir, a.name + ".out"), errors="replace").read()
    err = open(os.path.join(a.logdir, a.name + ".err"), errors="replace").read()
    why = []
    if st != 0:
        why.append("exit status %d" % st)
    for txt in (out, err):
        for line in txt.splitlines():
            if not a.no_fail_words and FAIL_RX.search(line) and not re.search(r"\bFAIL=0\b|\b0 FAIL", line):
                why.append("printed: " + line.strip()[:160])
                break
            if not a.nocounts and not line.lstrip().startswith("PASS"):
                m = COUNT_RX.search(line)
                if m:
                    why.append("reports: " + line.strip()[:160])
                    break
    for pat in a.expect:
        if not re.search(pat, out, re.M):
            why.append("missing expected output /%s/" % pat)
    for f in a.must_write:
        p = f if os.path.isabs(f) else os.path.join(a.cwd, f)
        if not os.path.exists(p) or os.path.getmtime(p) < t0:
            why.append("did not write " + f)
    # what did the step write, and is it the deposited file?
    comp = []
    for sc in a.scan:
        work, dep = sc.split("=", 1)
        for d, dirs, files in os.walk(work):
            dirs[:] = [x for x in dirs if x not in ("__pycache__",)]
            for f in files:
                p = os.path.join(d, f)
                try:
                    if os.path.getmtime(p) < t0 or os.path.islink(p):
                        continue
                except OSError:
                    continue
                rel = os.path.relpath(p, work)
                if any(re.search(x, rel) for x in a.ignore):
                    ok, v = None, "not compared (--ignore)"
                else:
                    ok, v = compare.same(p, os.path.join(dep, rel), a.strip + [work, dep])
                comp.append(dict(file=os.path.relpath(p, a.cwd) if p.startswith(a.cwd) else p, rel=rel, ok=ok, verdict=v,
                                 bytes=os.path.getsize(p)))
    for dep in a.stdout_vs:                                          # the printed output against a deposited log
        ok, v = compare.same(os.path.join(a.logdir, a.name + ".out"), dep, a.strip + [a.cwd])
        comp.append(dict(file=dep, rel="(printed output) vs " + os.path.basename(dep), ok=ok, verdict=v, bytes=len(out)))
    comp.sort(key=lambda c: c["rel"])
    ndiff = sum(1 for c in comp if c["ok"] is False)
    nsame = sum(1 for c in comp if c["ok"] is True)
    nnew = sum(1 for c in comp if c["ok"] is None)
    status = "failed" if why else "passed"
    rec.update(status=status, why=why, seconds=round(dt, 1), peak_rss_mib=round(peak / 2**20, 1), exit=st,
               written=len(comp), same=nsame, differ=ndiff, new=nnew,
               differing=[c["rel"] + ": " + c["verdict"] for c in comp if c["ok"] is False][:50])
    last = [l for l in out.splitlines() if l.strip()][-1:] or [""]
    print("     %s in %.0f s, peak memory %.0f MiB. Last line: %s" % (status.upper() if why else "ok", dt, peak / 2**20, last[0][:150]))
    for w in why[:4]:
        print("     why: " + w)
    if comp:
        print("     wrote %d file(s): %d the same as deposited, %d different, %d not compared (no deposited copy, or a rebuilt program)" % (len(comp), nsame, ndiff, nnew))
        for c in comp:
            if c["ok"] is False:
                print("       %s: %s" % (c["rel"], c["verdict"]))
    if why:
        for l in (out.splitlines()[-12:] + err.splitlines()[-12:]):
            print("       | " + l[:200])
    with open(a.summary, "a") as fh:
        fh.write(json.dumps(rec) + "\n")
    return 1 if why else 0


if __name__ == "__main__":
    sys.exit(main())
