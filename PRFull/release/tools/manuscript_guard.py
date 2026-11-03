#!/usr/bin/env python3
"""manuscript_guard.py -- find manuscript text in the files of a release.

    python3 manuscript_guard.py --root DIR --list FILES.txt [--sources manuscript_sources.json] [--min-run 10]

Compares every text file named in FILES.txt (paths relative to DIR) with manuscripts and drafts, which are not part of
a release, and reports, per file and per group of sources, the words it shares with
them in runs of at least --min-run consecutive words and its longest such run (after stripping LaTeX markup and
lower-casing; \\input and \\include are resolved).  Every .tex file and every other text file of up to 1 MB is compared
whole, so that a legend amid a table of numbers or formulas is found; a larger file is compared only if it is mostly
prose (is_prose), which spares the large tables of numbers.  Also reports files whose names alone mark them as
manuscripts and PDFs whose first page carries a listed title.

Which manuscripts to compare with is read from a local file (default: manuscript_sources.json beside this script),
which names them by glob patterns relative to DIR in groups, and lists titles; it is kept on the author's machine and
is not deposited.  Without it the guard checks file names only.  assemble_release.py calls scan() on its selection.
"""
import argparse
import glob
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
NAME_RE = re.compile(r"(^|/)(ms|cs|Main|SI)\.(tex|pdf)$|(^|/)refs\.tex$")
TEXT_EXT = (".md", ".txt", ".tex", ".py", ".sh", ".c", ".h", ".f", ".f90", ".json", ".jsonl", ".tsv", ".csv", ".sbatch",
            ".slurm", ".out", ".log", ".rst", ".html", ".summary", ".patch", ".diff")
PROSE_GATE_BYTES = 1_000_000      # a text file larger than this (and not .tex) is compared only if it is mostly prose


def words(text, latex=True):
    if latex:
        text = re.sub(r"(?<!\\)%.*", " ", text)
        text = re.sub(r"\\(cite[pt]?|ref|eqref|label|SIsec)\*?(\[[^\]]*\])?\{[^}]*\}", " ", text)
        text = re.sub(r"\\[A-Za-z]+", " ", text)
    text = re.sub(r"\[[^\]\n]{1,30}\]", " ", text)          # editorial insertions in quotations, "[error]"
    text = text.lower().replace("\u2019", "'")
    return re.sub(r"[^a-z0-9']+", " ", text).split()


def read_tex(path, depth=0):
    """the text of a .tex file with \\input and \\include resolved (relative to its directory)"""
    t = open(path, encoding="utf-8", errors="replace").read()
    if depth > 5:
        return t
    def sub(m):
        name = m.group(2)
        p = os.path.join(os.path.dirname(path), name if name.endswith(".tex") else name + ".tex")
        return read_tex(p, depth + 1) if os.path.isfile(p) else m.group(0)
    return re.sub(r"\\(input|include)\{([^}]+)\}", sub, t)


def is_prose(ws):
    """at least a third of the words alphabetic and longer than two letters: skips number tables and code dumps"""
    if len(ws) < 30:
        return False
    return sum(1 for w in ws if w.isalpha() and len(w) > 2) > len(ws) / 3


def load_config(path):
    if path and os.path.isfile(path):
        return json.load(open(path))
    return {"groups": {}, "exclude": [], "titles": []}


def load_refs(root, n, cfg):
    """{n-gram: (group, source)} over every source file of every group, and the list of sources read"""
    grams, found = {}, []
    for group, pats in cfg.get("groups", {}).items():
        for pat in pats:
            for p in sorted(glob.glob(os.path.join(root, pat), recursive=True)):
                rel = os.path.relpath(p, root)
                if not os.path.isfile(p) or any(glob.fnmatch.fnmatch(rel, x) for x in cfg.get("exclude", [])):
                    continue
                found.append(rel)
                ws = words(read_tex(p) if p.endswith(".tex") else open(p, encoding="utf-8", errors="replace").read(),
                           latex=p.endswith(".tex"))
                for i in range(len(ws) - n + 1):
                    g = tuple(ws[i:i + n])
                    if sum(1 for w in g if w.isalpha() and len(w) > 2) >= (n + 1) // 2:   # prose, not numbers or strings
                        grams.setdefault(g, (group, rel))
    return grams, found


def looks_text(p):
    try:
        with open(p, "rb") as fh:
            b = fh.read(4096)
        b.decode("utf-8")
        return b"\0" not in b
    except Exception:
        return False


def scan(root, rels, min_run=10, max_bytes=8_000_000, sources=None, refs_root=None):
    """hits: list of (path, kind, detail, info); info for kind 'text' is {group: {"words", "longest", "sources"}}"""
    cfg = load_config(sources or os.path.join(HERE, "manuscript_sources.json"))
    grams, found = load_refs(refs_root or root, min_run, cfg)   # refs_root: where the manuscripts are, if not under root
    titles = [t.lower() for t in cfg.get("titles", [])]
    hits = []
    for r in rels:
        p = os.path.join(root, r)
        if NAME_RE.search(r):
            hits.append((r, "name", "a manuscript file name", {}))
        if r.lower().endswith(".pdf") and os.path.isfile(p) and titles:
            try:
                t = subprocess.run(["pdftotext", "-l", "1", p, "-"], capture_output=True, text=True, timeout=60).stdout
            except Exception:
                t = ""
            t = " ".join(t.lower().split())
            for ti in titles:
                if ti in t:
                    hits.append((r, "pdf-title", ti, {}))
        if not grams or not os.path.isfile(p) or os.path.getsize(p) > max_bytes:
            continue
        ext = os.path.splitext(r)[1].lower()
        if not (r.lower().endswith(TEXT_EXT) or (ext == "" and looks_text(p))):
            continue
        txt = open(p, encoding="utf-8", errors="replace").read()
        ws = words(txt, latex=r.endswith(".tex"))
        if not r.endswith(".tex") and os.path.getsize(p) > PROSE_GATE_BYTES and not is_prose(ws):
            continue
        info, i = {}, 0
        while i <= len(ws) - min_run:
            g = tuple(ws[i:i + min_run])
            if g in grams:
                grp, src = grams[g]
                j = i + min_run
                while j < len(ws) and tuple(ws[j - min_run + 1:j + 1]) in grams:
                    j += 1
                d = info.setdefault(grp, {"words": 0, "longest": 0, "sources": set()})
                d["words"] += j - i
                d["longest"] = max(d["longest"], j - i)
                d["sources"].add(src)
                i = j
            else:
                i += 1
        if info:
            total = sum(d["words"] for d in info.values())
            hits.append((r, "text", "%d words in runs of >= %d shared with %s" % (
                total, min_run, "; ".join("%s (%d, longest %d)" % (g, d["words"], d["longest"]) for g, d in sorted(info.items()))),
                info))
    return hits, found


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--root", required=True)
    ap.add_argument("--list", required=True)
    ap.add_argument("--sources", default=None)
    ap.add_argument("--refs-root", default=None, help="where the manuscripts are (default --root), e.g. the author's tree "
                    "when --root is an assembled release")
    ap.add_argument("--min-run", type=int, default=10)
    a = ap.parse_args()
    rels = [l.strip() for l in open(a.list) if l.strip()]
    hits, found = scan(a.root, rels, a.min_run, sources=a.sources, refs_root=a.refs_root)
    print("sources compared: %d files%s" % (len(found), "" if found else " (no manuscript_sources.json: names only)"))
    print("files checked: %d; hits: %d" % (len(rels), len(hits)))
    for r, kind, d, _ in sorted(hits, key=lambda h: (h[1], h[0])):
        print("  [%s] %s: %s" % (kind, r, d))
    return 1 if hits else 0


if __name__ == "__main__":
    sys.exit(main())
