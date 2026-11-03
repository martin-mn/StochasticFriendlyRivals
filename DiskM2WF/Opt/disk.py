"""
NE.pdf's disk, so that DiskM2IN's panels can be laid beside NE's.

    (u,v)  ->  (u,v)/sqrt(k^2+u^2+v^2),   k = 4

carries R^2 onto the open unit disk, the directions at infinity onto the rim,
fixes the origin and keeps every line through the origin straight; every other
line becomes a conic arc.  Games with |(u,v)| <= k take exactly half the disk's
area.
"""
import math
import numpy as np

K = 4.0


def to_disk(u, v, k=K):
    u = np.asarray(u, dtype=float); v = np.asarray(v, dtype=float)
    f = np.sqrt(k * k + u * u + v * v)
    return u / f, v / f


def from_disk(x, y, k=K):
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    r2 = x * x + y * y
    f = k / np.sqrt(np.maximum(1.0 - r2, 1e-300))
    return x * f, y * f


def line_on_disk(a, b, c, k=K, n=4001):
    """the image of { a u + b v = c }, a conic arc from rim to rim.

    The line is parameterised t = k*tan(s), s in (-pi/2, pi/2), so the
    samples are dense near the origin -- where the disk map has all its
    curvature -- and still reach the rim as s -> +-pi/2.  Scaling t by k
    rather than by a huge constant is the whole point: with span = 1e7 the
    sample nearest the origin sat at |t| = 7854, whose image is already on
    the rim, so the "arc" was drawn as two straight chords and missed the
    true point (0.2425, 0) on u+v=1 by 0.17 disk radii.
    """
    s = np.linspace(-math.pi / 2 + 1e-12, math.pi / 2 - 1e-12, n)
    t = k * np.tan(s)
    nrm = math.hypot(a, b)
    ux, uy = -b / nrm, a / nrm             # direction along the line
    px, py = a * c / nrm**2, b * c / nrm**2   # the point closest to the origin
    return to_disk(px + t * ux, py + t * uy, k)


WEDGES = {
    "W": dict(cond=lambda u, v: (u + v < 1) & (u - v < 1),
              arc=(135.0, 225.0), nfr=8,    eff="reach R",   riv="Riv+"),
    "S": dict(cond=lambda u, v: (u + v < 1) & (u - v > 1),
              arc=(225.0, 315.0), nfr=1519, eff="reach R",   riv="Riv-"),
    "E": dict(cond=lambda u, v: (u + v > 1) & (u - v > 1),
              arc=(-45.0, 45.0),  nfr=80,   eff="alternate", riv="Riv-"),
    "N": dict(cond=lambda u, v: (u + v > 1) & (u - v < 1),
              arc=(45.0, 135.0),  nfr=80,   eff="alternate", riv="Riv+"),
}


def wedge_of(u, v):
    for k, d in WEDGES.items():
        if d["cond"](np.asarray(u), np.asarray(v)):
            return k
    return "boundary"


def draw_frame(ax, k=K, lines=True, labels=True, quad=True, lw=1.0):
    """the rim, the two diameters, the switch line and T = S"""
    th = np.linspace(0, 2 * math.pi, 721)
    ax.plot(np.cos(th), np.sin(th), color="0.15", lw=1.2, zorder=5)
    ax.plot([-1, 1], [0, 0], color="0.55", lw=0.6, zorder=4)
    ax.plot([0, 0], [-1, 1], color="0.55", lw=0.6, zorder=4)
    if lines:
        x, y = line_on_disk(1, 1, 1, k)         # u + v = 1, the switch
        ax.plot(x, y, color="#c1121f", lw=lw + 0.4, zorder=6)
        x, y = line_on_disk(1, -1, 1, k)        # u - v = 1, T = S
        ax.plot(x, y, color="#0353a4", lw=lw + 0.4, ls=(0, (5, 2)), zorder=6)
    if quad:
        for lab, (px, py) in (("PD", (-0.72, 0.72)), ("SD", (0.72, 0.72)),
                              ("SH", (-0.72, -0.72)), ("HA", (0.72, -0.72))):
            ax.text(px, py, lab, ha="center", va="center", fontsize=7,
                    color="0.45", zorder=7)
    if labels:
        for w, d in WEDGES.items():
            a = math.radians(0.5 * (d["arc"][0] + d["arc"][1]))
            ax.text(1.13 * math.cos(a), 1.13 * math.sin(a),
                    "%s\n%d" % (w, d["nfr"]), ha="center", va="center",
                    fontsize=7.5, color="#333333", zorder=7)
    ax.set_xlim(-1.28, 1.28); ax.set_ylim(-1.28, 1.28)
    ax.set_aspect("equal"); ax.axis("off")
