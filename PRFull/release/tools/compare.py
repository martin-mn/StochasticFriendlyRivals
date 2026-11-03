#!/usr/bin/env python3
"""compare.py -- is a regenerated file the same as the deposited one?

    python3 compare.py NEW OLD [--strip PATH ...]        one pair; prints the verdict, exit 0 if the same, 1 if not

"The same" is decided by the kind of file, and only what cannot be reproduced is set aside:
  PDF        byte-identical once /CreationDate, /ModDate and /ID are removed (matplotlib writes the date)
  PNG        byte-identical, or the same pixels
  JSON       the same data once every key that records a time, a date or a command line is dropped
             (seconds, timing, generated, date, regenerate, elapsed, ...) and paths are replaced (--strip)
  .gz        compared after decompression, as the inner file
  .npz/.npy  the same arrays, exactly
  .pkl       the same data once loaded (plain data only), whatever the pickle protocol or memo layout
  text       the same lines once run times ("12.3 s", "in 45 s", "0.8 ms", "row 0.15"), dates, clock times and the
             paths given by --strip are replaced by placeholders; lines written by several threads may come in
             another order (said so in the verdict)
  executable a compiled program or library (Mach-O, ELF): not compared -- a rebuild is not byte-identical in
             general (linker UUIDs, paths, compiler versions); what it computes is checked by the steps that use it
  otherwise  byte-identical
Used by runstep.py after every step of the reproduce_*.sh scripts.
"""
import gzip
import io
import json
import os
import pickle
import re
import sys

KEY_RX = re.compile(r"(^|_)(seconds?|secs?|time|times|timing|timings|generated|date|regenerate|elapsed|wall|host|"
                    r"hostname|cwd|command|cmd|started|finished|user|nproc|threads?)(_|$)|seconds|timing|^_seconds$|"
                    r"(^|_)(us|ms|ns)_per_|_per_(pair|point|call|second)|(^|_)rate_s$", re.I)
TIME_RX = re.compile(r"(?<![\w.])\d+(\.\d+)?(e[+-]?\d+)?\s?(s|sec|secs|seconds|ms|us|µs|min|mins|minutes|h)\b(?!\w)"
                     r"(?![ \t]*[=:]?[ \t]*\d)")          # "q=13362 seconds 85.7": "seconds" is a key there, not a unit
CLOCK_RX = re.compile(r"\b\d{4}-\d\d-\d\d([ T]\d\d:\d\d(:\d\d(\.\d+)?)?Z?)?\b|\b\d\d:\d\d:\d\d(\.\d+)?\b")
KEYTIME_RX = re.compile(r"\b(row|dp|mark|elapsed|wall|time|timing|seconds|secs|us_per_pair|ms_per_pair|full-search)([ =:]+)\d+(\.\d+)?\b")
TIMING_LINE_RX = re.compile(r"^(\s*)(TIMING|timing|Timing)\b.*$", re.M)
PDF_RX = [re.compile(rb"/(CreationDate|ModDate) ?\(D:[^)]*\)"), re.compile(rb"/ID ?\[ ?<[0-9A-Fa-f]*> ?<[0-9A-Fa-f]*> ?\]")]
TEXT_EXT = (".txt", ".out", ".tsv", ".csv", ".log", ".summary", ".md", ".tex", ".json", ".res", ".err", ".dat")


EXTRA_NORM = []                                              # --normalize: step-specific patterns (runstep.py)


def _norm_text(s, strip):
    for rx in EXTRA_NORM:
        s = re.sub(rx, "<normalized>", s)
    for p in sorted(strip, key=len, reverse=True):
        if p:
            s = s.replace(p, "<ROOT>")
    s = TIMING_LINE_RX.sub(lambda m: m.group(1) + "<timing line>", s)
    s = CLOCK_RX.sub("<clock>", s)
    s = TIME_RX.sub("<t>", s)
    s = KEYTIME_RX.sub(lambda m: m.group(1) + m.group(2) + "<t>", s)
    return s


def _drop(obj, strip):
    if isinstance(obj, dict):
        return {k: _drop(v, strip) for k, v in obj.items() if not KEY_RX.search(str(k))}
    if isinstance(obj, list):
        return [_drop(v, strip) for v in obj]
    if isinstance(obj, str):
        return _norm_text(obj, strip)
    return obj


def _read(path):
    b = open(path, "rb").read()
    if path.endswith(".gz"):
        b = gzip.decompress(b)
        path = path[:-3]
    return path, b


MAGIC = (b"\xcf\xfa\xed\xfe", b"\xce\xfa\xed\xfe", b"\xfe\xed\xfa\xcf", b"\xfe\xed\xfa\xce", b"\xca\xfe\xba\xbe",
         b"\x7fELF")


def is_binary_program(path):
    try:
        with open(path, "rb") as fh:
            return fh.read(4) in MAGIC
    except OSError:
        return False


class _Plain(pickle.Unpickler):
    """reads a pickle of plain data only (lists, tuples, dicts, strings, numbers): no class is looked up, so no code runs"""
    def find_class(self, module, name):
        raise pickle.UnpicklingError("holds %s.%s, not plain data" % (module, name))


def same(new, old, strip=()):
    """(bool or None, verdict); None: nothing to compare with, or a rebuilt program (not compared)"""
    if os.path.exists(new) and is_binary_program(new):
        return None, "rebuilt program (compiled code is not compared)"
    if not os.path.exists(old):
        return None, "new file (no deposited copy)"
    if not os.path.exists(new):
        return False, "NOT WRITTEN"
    rn, a = _read(new)
    ro, b = _read(old)
    if a == b:
        return True, "identical"
    ext = os.path.splitext(rn)[1].lower()
    if ext == ".pdf":
        x, y = a, b
        for rx in PDF_RX:
            x, y = rx.sub(b"", x), rx.sub(b"", y)
            if x == y:
                return True, "identical apart from the PDF creation date" + ("" if rx is PDF_RX[0] else " and /ID")
        return False, "DIFFERS (PDF content)"
    if ext == ".png":
        try:
            import matplotlib.image as mpimg
            import numpy as np
            if np.array_equal(mpimg.imread(io.BytesIO(a), format="png"), mpimg.imread(io.BytesIO(b), format="png")):
                return True, "identical pixels"
        except Exception as e:                                   # noqa: BLE001
            return False, "DIFFERS (PNG, %s)" % e
        return False, "DIFFERS (PNG pixels)"
    if ext in (".npz", ".npy"):
        import numpy as np
        A, B = np.load(io.BytesIO(a), allow_pickle=False), np.load(io.BytesIO(b), allow_pickle=False)
        if ext == ".npy":
            return (True, "identical arrays") if np.array_equal(A, B) else (False, "DIFFERS (array)")
        if sorted(A.files) != sorted(B.files):
            return False, "DIFFERS (npz members)"
        bad = [k for k in A.files if not np.array_equal(A[k], B[k])]
        return (True, "identical arrays") if not bad else (False, "DIFFERS (npz arrays %s)" % ",".join(bad[:4]))
    if ext == ".pkl":
        try:
            A, B = _Plain(io.BytesIO(a)).load(), _Plain(io.BytesIO(b)).load()
        except Exception as e:                                   # noqa: BLE001
            return False, "DIFFERS (pickle, %s)" % e
        return (True, "identical data (pickle)") if A == B else (False, "DIFFERS (pickle data)")
    if ext == ".json":
        try:
            ja, jb = json.loads(a), json.loads(b)
        except ValueError:
            ja = jb = None
        if ja is not None:
            if _drop(ja, strip) == _drop(jb, strip):
                return True, "identical apart from times, dates and paths"
            return False, "DIFFERS (JSON data)"
    if ext in TEXT_EXT or ext == "":
        try:
            ta, tb = a.decode("utf-8"), b.decode("utf-8")
        except UnicodeDecodeError:
            return False, "DIFFERS (bytes)"
        la, lb = _norm_text(ta, strip).splitlines(), _norm_text(tb, strip).splitlines()
        if la == lb:
            return True, "identical apart from times, dates and paths"
        if sorted(la) == sorted(lb):
            return True, "the same lines in another order (threads), apart from times, dates and paths"
        nd = sum(1 for x, y in zip(la, lb) if x != y) + abs(len(la) - len(lb))
        first = next((i for i, (x, y) in enumerate(zip(la, lb)) if x != y), min(len(la), len(lb)))
        return False, "DIFFERS (%d of %d lines; first at line %d)" % (nd, max(len(la), len(lb)), first + 1)
    return False, "DIFFERS (bytes)"


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("new")
    ap.add_argument("old")
    ap.add_argument("--strip", action="append", default=[])
    a = ap.parse_args()
    ok, v = same(a.new, a.old, a.strip)
    print(v)
    sys.exit(0 if ok else 1)
