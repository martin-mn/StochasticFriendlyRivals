#!/usr/bin/env python3
"""Shared furniture for the projective-disk figures -- p3disk.py, p3ray.py.

THE MAP.      (u,v)  ->  (u,v) / sqrt(k^2 + u^2 + v^2)

carries the whole plane onto the open unit disk, the rim being the circle of
directions at infinity.  Lines through the origin stay straight, so u = 0 and
v = 0 remain diameters and the four game families remain clean quadrants --
which is what lets each family be named on the rim of its own quadrant.
"""
import numpy as np
from matplotlib import font_manager, ft2font
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D

# the four families, by the sign quadrant of (u,v) -- the payoffs are
# (R,S,T,P) = (1, u, 1+v, 0), so u > 0 is S > P and v > 0 is T > R.  Same
# convention as NashGrid/census.py.  Each name is centred on its diagonal.
NAMES = ((45.0, "Snowdrift"),            # u > 0, v > 0
         (135.0, "Prisoner's Dilemma"),  # u < 0, v > 0
         (225.0, "Stag Hunt"),           # u < 0, v < 0
         (315.0, "Harmony"))             # u > 0, v < 0
NM_FS, NM_TRK, NM_GAP = 42.0, 0.09, 0.028      # points, em of tracking, radii


def phi(p, k):
    p = np.atleast_2d(np.asarray(p, dtype=float))
    return p / np.sqrt(k * k + (p ** 2).sum(axis=1))[:, None]


def arc(A, B, f, tol, depth=0):
    """A straight edge in the plane is a conic in the disk: subdivide the
    chord A--B until the drawn polyline is within `tol` of its true image.

    `f` is the map to use -- passed in rather than assumed, because p3disk.py
    wraps disklib.phi to offer the donation-cost radius as well, and that
    profile does NOT send lines to conics (see project_state.md section
    23.13).  The subdivision is agnostic: it measures the sag of the actual
    image, whatever the map.
    """
    a, b = np.asarray(A, float), np.asarray(B, float)
    fa, fb = f(a)[0], f(b)[0]
    m = 0.5 * (a + b)
    if depth >= 16 or np.hypot(*(f(m)[0] - 0.5 * (fa + fb))) < tol:
        return [fa]
    return arc(a, m, f, tol, depth + 1) + arc(m, b, f, tol, depth + 1)


def conic(a, b, c, f, rb, arms=(-1, 1), npt=20000):
    """The image of the line a*u + b*v = c, as one polyline per arm.

    The plane is cut at |t| = rb, whose image sits inside the rim by ~k^2/rb^2.
    arms=(1,) draws a half-line -- which is what the donation family is, since
    a cost is positive.
    """
    d = np.array([-b, a], dtype=float)
    p0 = np.array([a, b], dtype=float) * c / (a * a + b * b)
    t = np.linspace(0, rb, npt)[:, None]
    return [f(p0[None, :] + s * t * d[None, :]) for s in arms]


def discrete_bands(vals, weights, cmap="magma", vmin=None, vmax=None,
                   minshare=0.010, group=str, by="value", share=True):
    """One colour band per distinct level of a piecewise-constant field.

    The bar then ENUMERATES instead of implying a continuum, which is the only
    reason a discrete version exists (project_state.md section 23.18a).

    `vals` are the field's values over the cells that carry a value and
    `weights` how much of the drawn disk each cell occupies -- face area for a
    vector page, pixel count for a raster one.  That is the only difference
    between the two families, so the logic lives here and not in either.

    Band boundaries sit at the GEOMETRIC midpoints of consecutive levels,
    which is the right midpoint for a quantity read on a log scale, with the
    outer two placed symmetrically so the extreme bands are not half-width.

    WHICH LEVELS GET A TICK.  Ranked by the share of the disk they cover and
    taken greedily from the top, skipping any within MINSEP bands of one
    already taken; the two extremes are always labelled.  Each label carries
    the level AND its share, because "important" here means "occupies area"
    and that is the one number a reader cannot get from the picture.  MINSEP
    is derived from the band count rather than picked -- at a hand-chosen 3
    the third-largest level by area was silently crowded out by an extreme --
    and the receipt names the largest levels that missed a tick, so a loss of
    that kind is visible.

    `by` chooses what the band COLOURS encode.  "value" is the default and the
    original behaviour -- colour is the log of the level, so a level twice
    another is a fixed distance along the ramp.  "rank" spreads the colormap
    uniformly over the bands instead, so adjacent levels are always as distinct
    as the map allows.  Use it when the levels bunch: on Fig 3's limit panel
    two thirds of the disk's levels land in the orange third of magma under
    "value", because the field's levels are not uniform in log either.  Under
    "rank" colour no longer encodes magnitude proportionally -- the bar's own
    labels do that, and the bar has uniform band heights anyway, so on a
    discrete bar "rank" simply makes colour agree with position.

    `share=False` drops the "(x %)" from every label, leaving the bare level.

    Returns (ListedColormap, BoundaryNorm, ticks, ticklabels, receipt).
    """
    assert by in ("value", "rank"), "by must be value or rank"
    from matplotlib.colors import BoundaryNorm, ListedColormap, LogNorm
    from matplotlib import pyplot as plt

    vals = np.asarray(vals)
    lv, inv = np.unique(vals, return_inverse=True)
    ash = np.zeros(len(lv))
    np.add.at(ash, inv, np.asarray(weights, dtype=float))
    ash /= ash.sum()
    g = np.sqrt(lv[:-1] * lv[1:])
    bnd = np.concatenate(([lv[0] ** 2 / g[0]], g, [lv[-1] ** 2 / g[-1]]))
    lo = float(vals.min()) if vmin is None else vmin
    hi = float(vals.max()) if vmax is None else vmax
    if by == "rank":
        t = np.linspace(0.0, 1.0, len(lv)) if len(lv) > 1 else np.zeros(1)
        cols = plt.get_cmap(cmap)(t)
    else:
        cols = plt.get_cmap(cmap)(LogNorm(lo, hi)(lv))

    minsep = max(2, int(round(len(lv) / 38.0)))
    keep, lost = [0, len(lv) - 1], []
    for i in np.argsort(ash)[::-1]:
        if ash[i] < minshare:
            break
        if all(abs(i - j) >= minsep for j in keep):
            keep.append(int(i))
        elif int(i) not in keep:
            lost.append(int(i))
    keep = sorted(keep)
    receipt = ("discrete bar: %d bands, %d labelled (>= %.1f%% of the disk, at "
               "least %d bands apart); labelled levels cover %.1f%% of it%s"
               % (len(lv), len(keep), 100 * minshare, minsep,
                  100 * sum(ash[i] for i in keep),
                  "" if not lost else
                  ".  Largest unlabelled: %s"
                  % ", ".join("%g (%.1f%%)" % (lv[i], 100 * ash[i])
                              for i in lost[:3])))
    return (ListedColormap(cols), BoundaryNorm(bnd, len(lv)),
            [lv[i] for i in keep],
            [("%s  (%.1f%%)" % (group(lv[i]), 100 * ash[i])) if share
             else group(lv[i]) for i in keep],
            receipt)


def advances(s, fp, fontsize, tracking):
    """Pen advance of every character of s, in points, kerned, plus tracking."""
    face = ft2font.FT2Font(font_manager.findfont(fp))
    face.set_size(fontsize, 72.0)            # 72 dpi, so one unit is one point
    adv, prev = [], None
    for ch in s:
        idx = face.get_char_index(ord(ch))
        if prev is not None:
            k = face.get_kerning(prev, idx, ft2font.KERNING_DEFAULT)
            adv[-1] += k / 64.0 if k else 0.0
        adv.append(face.load_char(ord(ch)).linearHoriAdvance / 65536.0
                   + tracking * fontsize)
        prev = idx
    return adv


def arc_text(ax, s, theta_c, r_in, fontsize, dpu, gap, tracking,
             color="0.10", zorder=12):
    """Set s on an arc outside r_in, centred on theta_c (degrees).

    dpu is data units per point.  On the upper half the letters stand radially
    outward and the reading direction is decreasing theta; on the lower half
    they stand radially inward and it is increasing theta, so both halves read
    the right way up and the nearest ink is `gap` from the rim either way.
    Returns (theta_lo, theta_hi, r_outer) for the receipts.
    """
    fp = FontProperties(size=fontsize)
    adv = advances(s, fp, fontsize, tracking)
    total = sum(adv) - tracking * fontsize            # drop trailing tracking
    ink = TextPath((0, 0), s, size=fontsize, prop=fp).get_extents()
    asc, desc = ink.y1, max(-ink.y0, 0.0)             # points above / below
    up = (theta_c % 360.0) < 180.0
    r0 = r_in + gap + (0.0 if up else asc * dpu)      # the baseline radius
    half = 0.5 * total * dpu / r0                     # radians
    th_c, x = np.deg2rad(theta_c), 0.0
    for ch, a in zip(s, adv):
        if ch != " ":
            d = (x + 0.5 * a) * dpu / r0              # centre of this slot
            th = th_c + (half - d) if up else th_c - half + d
            tr = (Affine2D().translate(-0.5 * a, 0.0).scale(dpu)
                  .rotate_deg(np.rad2deg(th) + (-90.0 if up else 90.0))
                  .translate(r0 * np.cos(th), r0 * np.sin(th)))
            ax.add_patch(PathPatch(
                tr.transform_path(TextPath((0, 0), ch, size=fontsize, prop=fp)),
                fc=color, ec="none", lw=0, zorder=zorder))
        x += a
    return (np.rad2deg(th_c - half), np.rad2deg(th_c + half),
            r0 + (asc if up else desc) * dpu)
