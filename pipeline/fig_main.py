#!/usr/bin/env python3
"""Main figures: expression of every panel gene across the stages, one figure per gene family and unit.

Usage: fig_main.py <config.yaml> <directory with the panel tables> <figure directory>

Each figure is a grid with one panel per gene of the family. A panel shows the mean normalized
TPM at each stage with its standard error: as bars when the unit's stages are categories (x_axis:
categorical), as a line over a numeric axis when the stages are time points (x_axis: numeric,
the number read from the stage name with numeric_from). Above a regulated gene the compact
letters from the letters table are placed over the error bars; stages that share a letter do not
differ significantly. Two styles are drawn: "manuscript" (serif font, gray bars, few ticks, and
optionally fixed y-axis maxima so that related panels share a scale) and "enhanced" (sans-serif
font, gray bars, one dot per animal). Family order, figure numbers, filenames, style
directories, display names and the species silhouette all come from the figures block of the
config. The silhouette sits in an empty grid cell when the family leaves one, otherwise in the
top-right margin. Output: one PNG (300 dpi) and one PDF (vector, fonts embedded) per figure.

BASE_LOCAL names the project directory under which the panel file and the silhouette are found;
by default it is the directory two levels above the directory that holds the config (a config at
<root>/pipeline/configs/x.yaml resolves to <root>). Set BASE_LOCAL=. when running from
the project root, as in the README.
"""
import os, sys, math, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib import font_manager as fm
from matplotlib.ticker import FuncFormatter, MaxNLocator, MultipleLocator
from PIL import Image

import yaml
if len(sys.argv) != 4: sys.exit("usage: fig_main.py <config.yaml> <output_dir with tables> <figure_dir>")
CFG = yaml.safe_load(open(sys.argv[1])); DATA, OUT = sys.argv[2], sys.argv[3]
ROOT = os.environ.get("BASE_LOCAL", os.path.abspath(os.path.join(os.path.dirname(sys.argv[1]), "..", "..")))
FIG = CFG["figures"]; UCOL = CFG["unit_col"]; UNITS = list(CFG["units"].keys())

# Register the Liberation fonts wherever the system keeps them, so that the figures look the same
# on every machine that has them. Without them matplotlib falls back to the next font in the lists
# set by style_rc (the data drawn are unaffected).
for f in fm.findSystemFonts():
    if os.path.basename(f).startswith(("LiberationSerif", "LiberationSans")):
        try: fm.fontManager.addfont(f)
        except Exception: pass
plt.rcParams.update({"pdf.fonttype":42, "ps.fonttype":42, "svg.fonttype":"none"})  # embed fonts as editable text in the PDFs

# ---- data ----
# The panel gives the family of each gene; the tables written by dge_panel.R give the means, the
# standard errors, the letters and the per-animal values.
panel = pd.read_csv(os.path.join(ROOT, CFG["panel"]), sep="\t").drop_duplicates("gene")
fam_order = FIG["family_order"]
NCOLS = {f: (FIG.get("ncols") or {}).get(f, 2) for f in fam_order}
DISPLAY_REMAP = FIG.get("display_names") or {}
def disp(g): return DISPLAY_REMAP.get(g, g)
genes_by_fam = {f: [g for g in panel.gene if panel.set_index("gene").family[g] == f] for f in fam_order}
if CFG["tables"]["order"] == "family_name":
    genes_by_fam = {f: sorted(gs, key=disp) for f, gs in genes_by_fam.items()}
TAG = CFG.get("panel_tag", "panel")   # names the per-unit panel tables
tpm = {u: pd.read_csv(f"{DATA}/{u}_{TAG}_tpm_by_stage.tsv", sep="\t") for u in UNITS}
let = {u: pd.read_csv(f"{DATA}/{u}_{TAG}_letters.tsv", sep="\t") for u in UNITS}
pts = pd.read_csv(f"{DATA}/{CFG['name']}_{TAG}_persample_tpm.tsv", sep="\t")
SILP = os.path.join(ROOT, FIG["silhouette"]); SIL = Image.open(SILP); SIL_LABEL = FIG["silhouette_label"]
SIL_ASPECT = SIL.size[0] / SIL.size[1]           # width / height of the silhouette image
SIL_W_IN   = 0.66                                 # silhouette width in inches, the same on every figure
FIG_NUM = {tuple(k.split("/")): v for k, v in FIG["figure_numbers"].items()}
FIXED_YLIM = {tuple(k.split("/")): v for k, v in (FIG.get("fixed_ylim") or {}).items()}
STYLE_DIRS = FIG.get("style_dirs") or {"manuscript": "manuscript", "enhanced": "enhanced"}
FNAME = FIG.get("filename", "Fig_{n:02d}_{family}_{unit}")
# Two styles written to the same directory under the same filename would overwrite each other.
_dirs = [STYLE_DIRS[s] for s in FIG["styles"]]
if len(set(_dirs)) != len(_dirs): sys.exit(f"figures.style_dirs: styles {FIG['styles']} share an output directory")

SZ = {"manuscript": dict(tick=12, axlab=13, title=14.5, letter=11, panel=20),
      "enhanced":   dict(tick=8.5, axlab=9.5, title=10.5, letter=8, panel=12)}

def nice_ceil(x):
    """Round a y-axis maximum up to a round number (1, 1.5, 2, 2.5, 3, 4, 5, 6, 8 or 10 times a power of ten)."""
    if not np.isfinite(x) or x <= 0: return 1.0
    mag = 10 ** math.floor(math.log10(x))
    for m in [1,1.5,2,2.5,3,4,5,6,8,10]:
        if x <= m*mag: return m*mag
    return 10*mag
def fmt_letters(s): return ",".join(list(str(s)))          # "ab" is shown as "a,b"
def is_true(v):     return str(v).strip().upper() == "TRUE"  # the tables spell booleans TRUE/FALSE
def dec_for(ymax):  return 1 if ymax >= 1 else (2 if ymax >= 0.1 else 3)   # decimals on the y-axis ticks

# ---- compact-letter (significance) placement ----
# Each letter starts a few points above its own error-bar cap (mean + SEM). It is then moved up by
# the smallest amount that clears any real overlap with a dot, a neighboring bar or letter, or the
# top of the frame. If the letters still do not fit, the y-axis is extended a little and the
# placement is solved again.
LETTER_BASE_PT = 3.0     # letter baseline seats this many points above the error-bar cap
LETTER_STEP_PT = 1.0     # minimal upward nudge increment (points)
LETTER_PAD_PT  = 1.5     # clearance padding around obstacles (points)

def _pad_box(e, p): return (e.x0 - p, e.y0 - p, e.x1 + p, e.y1 + p)
def _overlap(a, b): return not (a[2] <= b[0] or a[0] >= b[2] or a[3] <= b[1] or a[1] >= b[3])

def place_compact_letters(fig, ax, xs, caps, letters, dot_xy, bar_boxes_data, s):
    """Seat each compact-letter LETTER_BASE_PT points above its error-bar cap (mean+SEM),
    then nudge it up by the minimal amount needed to clear per-animal dots, neighboring
    bars/letters, and the top frame. Collisions are resolved in display space so the offset
    stays a fixed few points."""
    texts = [ax.text(xi, cap, fmt_letters(lt), ha="center", va="bottom",
                     fontsize=s["letter"], zorder=6, clip_on=False)
             for xi, cap, lt in zip(xs, caps, letters)]
    def set_off(t, pt):
        t.set_transform(mtransforms.offset_copy(ax.transData, fig=fig, x=0, y=pt, units="points"))
    dpi = fig.dpi
    pad_px = LETTER_PAD_PT * dpi / 72.0
    mrad   = math.sqrt(9.0/math.pi) * dpi / 72.0      # scatter s=9 marker radius in px
    for _ in range(5):                                # re-solve if the top frame needs more headroom
        for t in texts: set_off(t, LETTER_BASE_PT)
        fig.canvas.draw(); r = fig.canvas.get_renderer()
        dot_boxes = []
        for dx, dy in dot_xy:
            px, py = ax.transData.transform((dx, dy))
            dot_boxes.append((px - mrad, py - mrad, px + mrad, py + mrad))
        bar_boxes = []
        for (x0d, x1d, ytop) in bar_boxes_data:
            X0, Y0 = ax.transData.transform((x0d, 0.0)); X1, Y1 = ax.transData.transform((x1d, ytop))
            bar_boxes.append((min(X0, X1), min(Y0, Y1), max(X0, X1), max(Y0, Y1)))
        axtop = ax.get_window_extent(r).y1
        placed = []
        for i in sorted(range(len(texts)), key=lambda k: caps[k]):   # settle lowest cap first
            off = LETTER_BASE_PT
            for _step in range(400):
                b = _pad_box(texts[i].get_window_extent(r), pad_px)
                if not (any(_overlap(b, db) for db in dot_boxes) or
                        any(_overlap(b, bb) for bb in bar_boxes) or
                        any(_overlap(b, pb) for pb in placed)):
                    break
                off += LETTER_STEP_PT; set_off(texts[i], off)
            placed.append(_pad_box(texts[i].get_window_extent(r), pad_px))
        top_px = max(texts[i].get_window_extent(r).y1 for i in range(len(texts)))
        if top_px <= axtop - pad_px: break                           # fits under the frame -> done
        need = ax.transData.inverted().transform((0, top_px + pad_px))[1]   # grow ylim to fit the letters
        ax.set_ylim(0, need * 1.01)
    return texts

def style_rc(style):
    """Matplotlib settings for one style: serif and heavier lines for the manuscript, sans-serif and lighter lines otherwise."""
    if style == "manuscript":
        plt.rcParams.update({"font.family":"serif","font.serif":["Liberation Serif","Times New Roman","DejaVu Serif"],
            "axes.linewidth":1.2,"xtick.direction":"out","ytick.direction":"out","mathtext.default":"regular",
            "xtick.major.size":6,"ytick.major.size":6,"xtick.major.width":1.2,"ytick.major.width":1.2,
            "xtick.labelsize":SZ["manuscript"]["tick"],"ytick.labelsize":SZ["manuscript"]["tick"]})
    else:
        plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Liberation Sans","Arial","DejaVu Sans"],
            "axes.linewidth":0.8,"xtick.direction":"out","ytick.direction":"out","axes.grid":False,
            "xtick.major.size":3.5,"ytick.major.size":3.5,"xtick.major.width":0.8,"ytick.major.width":0.8,
            "xtick.labelsize":SZ["enhanced"]["tick"],"ytick.labelsize":SZ["enhanced"]["tick"]})

def silhouette_in_slot(fig, ax, s, figw, figh):
    """Seat the uniform-size silhouette centered in an empty panel slot (a family whose gene count
    leaves the last grid cell free). Same rendered size (SIL_W_IN) as the top-right placement."""
    ax.axis("off")
    box = ax.get_position()
    W_in = SIL_W_IN; H_in = W_in / SIL_ASPECT
    bw, bh = W_in/figw, H_in/figh
    cx, cy = box.x0 + box.width/2, box.y0 + box.height/2
    pax = fig.add_axes([cx - bw/2, cy - bh/2 + 0.06*bh, bw, bh])   # small up-nudge centers image+label as a pair
    pax.imshow(SIL); pax.axis("off")
    pax.text(0.5, -0.10, SIL_LABEL, transform=pax.transAxes, ha="center", va="top",
             fontstyle="italic", fontsize=s["tick"])

def silhouette_topright(fig, axB, s, figw, figh):
    """Seat the uniform-size silhouette high in the top-right margin above panel B (top edge near
    0.983), used when the grid has no empty slot."""
    box = axB.get_position()                            # panel B
    x0 = box.x1 + 0.004                                 # just right of panel B -> no interior overlap
    W_in = SIL_W_IN                                     # fixed size -> identical on every figure
    H_in = W_in / SIL_ASPECT
    bw, bh = W_in/figw, H_in/figh
    y0 = (0.915 + 0.068) - bh                           # top edge at the high corner
    pax = fig.add_axes([x0, y0, bw, bh]); pax.imshow(SIL); pax.axis("off")
    pax.text(0.5, -0.12, SIL_LABEL, transform=pax.transAxes, ha="center", va="top",
             fontstyle="italic", fontsize=s["tick"])

def draw(family, unit, style):
    """Draw and save one figure: the genes of one family in one unit, in one style."""
    genes = genes_by_fam[family]; n = len(genes)
    nc = NCOLS[family]; nr = math.ceil(n / nc)
    figw, figh = (8.3 if nc >= 4 else 7.0), 2.7*nr   # widen the compact 4-column grid (still one page wide)
    ws = 0.70 if nc >= 4 else 0.36                    # open the column gap so each y-axis label clears its neighbor
    s = SZ[style]; style_rc(style)
    unit_cfg = CFG["units"][unit]; stage_order = unit_cfg["levels"]
    fig, axes = plt.subplots(nr, nc, figsize=(figw, figh))
    fig.subplots_adjust(top=0.90, hspace=0.70, wspace=ws)
    axes = np.array(axes).reshape(-1)
    tp, le = tpm[unit], let[unit]
    x = np.arange(len(stage_order))
    for i in range(nr*nc):
        ax = axes[i]
        if i >= n: ax.axis("off"); continue
        gene = genes[i]
        sub  = tp[tp.gene == gene].set_index("stage").reindex(stage_order)
        lsub = le[le.gene == gene].set_index("stage").reindex(stage_order)
        means = sub.mean_ntpm.values.astype(float); sems = sub.sem_ntpm.values.astype(float)
        molt_reg = is_true(lsub.molt_regulated.iloc[0]); letters = lsub.letter.values
        pg = pts[(pts[UCOL] == unit) & (pts.gene == gene)]
        # The y-axis maximum: fixed by the config for the manuscript style when given (so related
        # panels share a scale), otherwise a round number above the tallest bar (plus dots in the
        # enhanced style).
        mm = np.nanmax(means + np.nan_to_num(sems))
        fixed = FIXED_YLIM.get((family, unit))
        if style == "manuscript" and fixed is not None:
            ymax = float(fixed[i])
        elif style == "manuscript":
            ymax = nice_ceil(mm * 1.05)
        else:
            dmax = float(np.nanmax(pg.ntpm.values)) if len(pg) else 0.0
            ymax = nice_ceil(max(mm, dmax) * 1.05)

        # Two kinds of x-axis. Categorical stages: bars at 0, 1, 2, ... Numeric stages: a line over the
        # numbers read from the stage names. The jitter of the per-animal dots uses a fixed seed so
        # that the figure is identical on every run.
        dot_xy = []                                             # (x,y) of every plotted per-animal dot
        bar_w = 0.72 if style == "manuscript" else 0.70
        if unit_cfg.get("x_axis", "categorical") == "categorical":
            if style == "manuscript":                               # ---- gray bars, black error bars ----
                ax.bar(x, means, color="#585858", edgecolor="black", lw=0.9, width=bar_w, zorder=2)
                ax.errorbar(x, means, yerr=sems, fmt="none", ecolor="black", capsize=3.5, lw=1.1, zorder=3)
            else:                                                   # ---- gray bars, per-animal dots ----
                ax.bar(x, means, color="#9E9E9E", edgecolor="none", width=bar_w, zorder=1)
                for xi, st in enumerate(stage_order):
                    yv = pg[pg.stage == st].ntpm.values
                    jit = (np.random.RandomState(xi+7).rand(len(yv)) - 0.5) * 0.28
                    xj = xi + jit
                    ax.scatter(xj, yv, s=9, color="black", alpha=0.75, edgecolor="none", zorder=3)
                    dot_xy.extend(zip(xj.tolist(), yv.tolist()))
                ax.errorbar(x, means, yerr=sems, fmt="none", ecolor="#2b2b2b", capsize=2, lw=0.9, zorder=4)
            ax.set_xticks(x); ax.set_xticklabels(stage_order, fontsize=s["tick"])
            ax.set_xlabel(unit_cfg.get("x_label", "Stage"), fontsize=s["axlab"])
        else:                                                   # ---- numeric axis: line chart ----
            import re as _re
            xd = np.array([int(_re.search(unit_cfg["numeric_from"], st).group(1)) for st in stage_order])
            if style == "manuscript":
                ax.plot(xd, means, color="#666666", lw=1.7, zorder=2)
                ax.errorbar(xd, means, yerr=sems, fmt="o", ms=6, color="black", ecolor="black", capsize=3.5, lw=1.1, zorder=3)
            else:
                for j, st in enumerate(stage_order):
                    yv = pg[pg.stage == st].ntpm.values
                    jit = (np.random.RandomState(j+7).rand(len(yv)) - 0.5) * 0.5
                    ax.scatter(xd[j]+jit, yv, s=8, color="#808080", alpha=0.45, edgecolor="none", zorder=2)
                    dot_xy.extend(zip((xd[j]+jit).tolist(), yv.tolist()))
                ax.plot(xd, means, color="black", lw=1.5, zorder=3)
                ax.errorbar(xd, means, yerr=sems, fmt="o", ms=4.5, color="black", ecolor="black", capsize=2, lw=0.9, zorder=4)
            ax.set_xlim(-1, xd.max() + 1); ax.set_xticks([t for t in range(0, int(xd.max()) + 1, 5)])
            ax.xaxis.set_minor_locator(MultipleLocator(1))
            ax.tick_params(axis="x", which="minor", length=(3 if style == "manuscript" else 2.2), width=(1.0 if style == "manuscript" else 0.7))
            ax.set_xlabel(unit_cfg.get("x_label", "Time"), fontsize=s["axlab"]); x = xd; bar_w = 0.0
        ax.set_ylim(0, ymax)
        if style == "manuscript":
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1,2,2.5,3,5,10]))   # fewer, cleaner ticks
            d = dec_for(ymax)
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _p, d=d: f"{v:.{d}f}"))
        ax.set_ylabel("Normalized TPM", fontsize=s["axlab"])
        ax.set_title(disp(gene), fontstyle="italic", fontsize=s["title"], pad=6)

        if molt_reg:                                            # compact-letter display over regulated panels
            caps = [m + (0.0 if np.isnan(se) else se) for m, se in zip(means, sems)]   # error-bar caps
            if unit_cfg.get("x_axis", "categorical") == "categorical":
                bar_boxes = [(xi - bar_w/2, xi + bar_w/2, m) for xi, m in zip(x, means)]
            else:                                               # line branch: only dots + neighbors obstruct
                bar_boxes = []
            place_compact_letters(fig, ax, x, caps, letters, dot_xy, bar_boxes, s)
        ax.spines[["top","right"]].set_visible(False)
        ax.text(-0.26, 1.22, chr(65+i), transform=ax.transAxes, fontsize=s["panel"],
                fontweight=("normal" if style == "manuscript" else "bold"), va="top", ha="left")

    # Silhouette placement: inside an empty grid slot when the family leaves one, otherwise snug to
    # panel B's top-right margin. The image and its label come from the config's figures block.
    if nr*nc - n > 0:
        silhouette_in_slot(fig, axes[nr*nc - 1], s, figw, figh)
    else:
        silhouette_topright(fig, axes[1], s, figw, figh)

    os.makedirs(os.path.join(OUT, STYLE_DIRS[style]), exist_ok=True)
    base = os.path.join(OUT, STYLE_DIRS[style], FNAME.format(n=FIG_NUM[(family, unit)], family=family, unit=unit))
    fig.savefig(base+".png", dpi=300, bbox_inches="tight")
    fig.savefig(base+".pdf", dpi=600, bbox_inches="tight")       # vector; embedded silhouette raster at 600 dpi
    plt.close(fig)

if __name__ == "__main__":
    for fam in fam_order:
        for unit in UNITS:
            for style in FIG["styles"]:
                draw(fam, unit, style)
    print("done:", OUT)
