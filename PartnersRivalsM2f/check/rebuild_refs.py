#!/usr/bin/env python3
"""rebuild_refs.py -- Text/refs_cited.tex from Text/refs.tex, in the order of first citation in the typeset document:
ms.tex with every \\input expanded where it stands (methods, tables, figure captions), so that the numbered list follows
the text.  Refuses to write if a cited key has no bibitem."""
import os, re, sys
T = "/Users/martin/Documents/CL1/PartnersRivalsM2f/Text"
def expand(path, seen):
    src = open(path, encoding="utf-8").read()
    src = re.sub(r"(?<!\\)%.*", "", src)                          # comments do not cite
    out = []; pos = 0
    for m in re.finditer(r"\\input\{([^}]+)\}", src):
        out.append(src[pos:m.start()]); pos = m.end()
        name = m.group(1)
        if name == "refs_cited": continue
        f = os.path.normpath(os.path.join(T, name if name.endswith(".tex") else name + ".tex"))
        out.append(expand(f, seen))
    out.append(src[pos:]); return "".join(out)
doc = expand(os.path.join(T, "ms.tex"), set())
order = []
for m in re.finditer(r"\\cite[pt]?\*?(?:\[[^\]]*\])*\{([^}]+)\}", doc):
    for k in m.group(1).split(","):
        k = k.strip()
        if k and k not in order: order.append(k)
refs = open(os.path.join(T, "refs.tex"), encoding="utf-8").read()
items = {m.group(1): m.group(2).strip() for m in re.finditer(r"\\bibitem\{([^}]+)\}\s*(.*?)(?=\\bibitem\{|\\end\{thebibliography\}|\Z)", refs, re.S)}
missing = [k for k in order if k not in items]
if missing: sys.exit("no bibitem for: " + ", ".join(missing))
with open(os.path.join(T, "refs_cited.tex"), "w", encoding="utf-8") as f:
    f.write("% refs_cited.tex -- written by check/rebuild_refs.py: the cited entries of refs.tex, in order of first citation\n\n")
    for k in order: f.write("\\bibitem{%s}\n%s\n\n" % (k, items[k]))
print("%d references, in order of first citation:" % len(order)); print("  " + ", ".join("%d %s" % (i + 1, k) for i, k in enumerate(order)))
