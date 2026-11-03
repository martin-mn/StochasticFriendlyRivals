#!/usr/bin/env python3
"""trace_reads.py -- run one Python script and record every file it opens under a root directory.

    python3 trace_reads.py --root ROOT --log LOG.json [--cwd DIR] SCRIPT [ARGS ...]

The script runs in this interpreter (runpy, as __main__) with an audit hook (PEP 578) that records the
events 'open' (every open(), io.open_code() for imports, numpy.load, ...), 'os.listdir', 'os.scandir',
'glob.glob', 'ctypes.dlopen' and 'subprocess.Popen'.  Only paths under ROOT are kept.  The log is JSON:
{"script", "argv", "seconds", "exit", "reads": [...], "writes": [...], "listed": [...], "dlopen": [...],
"subprocess": [...]}.  A read is a path opened for reading that exists as a regular file, given as the script sees it
(symbolic links are not resolved; "realpath" maps each read through a link to its target); a write is a
path opened with w, a, x or +.  --alias FROM=TO reports a shadow copy's paths as the originals'.  assemble_release.py uses these logs to find every data file a figure
script reads; nothing here decides what the release contains by itself.
"""
import argparse
import json
import os
import runpy
import sys
import time

ap = argparse.ArgumentParser()
ap.add_argument("--root", required=True)
ap.add_argument("--log", required=True)
ap.add_argument("--cwd", default=None)
ap.add_argument("--alias", action="append", default=[], help="FROM=TO: report paths under FROM as under TO (a shadow copy)")
ap.add_argument("script")
ap.add_argument("args", nargs=argparse.REMAINDER)
a = ap.parse_args()
ROOT = os.path.abspath(a.root)
ALIAS = [tuple(os.path.abspath(x) for x in s.split("=", 1)) for s in a.alias]
REC = {"reads": set(), "writes": set(), "listed": set(), "dlopen": set(), "subprocess": []}


def _keep(p):
    try:
        p = os.fsdecode(p)
    except Exception:
        return None
    if not isinstance(p, str) or not p:
        return None
    p = os.path.abspath(p)                  # the logical path, symbolic links NOT resolved (see realpath below)
    for f, t in ALIAS:
        if p == f or p.startswith(f + os.sep):
            p = t + p[len(f):]
            break
    return p if (p == ROOT or p.startswith(ROOT + os.sep)) else None


def hook(ev, args):
    try:
        if ev == "open":
            path, mode = args[0], args[1] if len(args) > 1 else "r"
            if isinstance(path, int):
                return
            p = _keep(path)
            if p is None:
                return
            mode = mode or "r"
            flags = args[2] if len(args) > 2 and isinstance(args[2], int) else 0
            wr = any(c in str(mode) for c in "wax+") or (flags & (os.O_WRONLY | os.O_RDWR | os.O_CREAT))
            (REC["writes"] if wr else REC["reads"]).add(p)
        elif ev in ("os.listdir", "os.scandir"):
            p = _keep(args[0] if args and args[0] is not None else ".")
            if p:
                REC["listed"].add(p)
        elif ev == "glob.glob":
            p = _keep(os.path.dirname(os.fsdecode(args[0])) or ".")
            if p:
                REC["listed"].add(p)
        elif ev == "ctypes.dlopen":
            p = _keep(args[0]) if args and args[0] else None
            if p:
                REC["dlopen"].add(p)
        elif ev == "subprocess.Popen":
            REC["subprocess"].append([os.fsdecode(x) for x in (args[1] if isinstance(args[1], (list, tuple)) else [args[1]])])
    except Exception:
        pass


if a.cwd:
    os.chdir(a.cwd)
script = os.path.abspath(a.script)
sys.argv = [script] + a.args
sys.path.insert(0, os.path.dirname(script))
sys.addaudithook(hook)
t0 = time.time()
code = 0
try:
    runpy.run_path(script, run_name="__main__")
except SystemExit as e:
    code = e.code if isinstance(e.code, int) else (0 if e.code is None else 1)
except BaseException as e:                                      # noqa: BLE001
    import traceback
    traceback.print_exc()
    code = 1
dt = time.time() - t0
sys.stdout.flush()
out = {"script": script, "argv": a.args, "seconds": round(dt, 1), "exit": code,
       "reads": sorted(p for p in REC["reads"] if os.path.isfile(p)),
       "realpath": {p: os.path.realpath(p) for p in sorted(REC["reads"]) if os.path.isfile(p) and os.path.realpath(p) != p},
       "writes": sorted(REC["writes"]), "listed": sorted(REC["listed"]), "dlopen": sorted(REC["dlopen"]),
       "subprocess": REC["subprocess"]}
with open(a.log, "w") as fh:
    json.dump(out, fh, indent=1)
sys.__stderr__.write("trace_reads: %s exit %s in %.1f s, %d reads, %d writes under %s\n"
                     % (os.path.basename(script), code, dt, len(out["reads"]), len(out["writes"]), ROOT))
sys.exit(code)
