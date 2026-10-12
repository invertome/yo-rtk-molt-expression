#!/usr/bin/env python3
"""Supplementary figures: volcano and MA plots of the pairwise comparisons between stages.

Usage: python fig_supp.py <config.yaml> <directory with the *_dge_*.tsv tables> <output directory>

For each entry of figures.supplement in the config, one figure is drawn for one unit and saved as
<filename>.png and <filename>.pdf (the entry's filename key, by default Fig<tag>_<unit>_dge). Each
row of the figure is one comparison (each stage against the previous one, the closing step of a
cycle when configured, or the extra comparisons of the config); the left column is a volcano plot
(log2 fold change against -log10 of the adjusted p-value) and the right column an MA plot (log2
fold change against mean expression). Every feature of the transcriptome is a gray point; the
panel genes are red points; a panel gene whose transcriptome-wide adjusted p-value is below alpha
in that comparison (the rtk_sig column of the table) is drawn as a star with its name.

The gray background clouds are rasterized inside the PDF to keep the files small; everything
else stays vector. The script also writes <name>_stars.json, which lists for each figure (by filename)
the (comparison, gene) pairs that received a star, so the starred genes can be checked against
the rtk_sig column of the tables. A missing table is an error; nothing is drawn from made-up data.
"""
import os, sys, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D

import json, yaml
if len(sys.argv) != 4: sys.exit("usage: fig_supp.py <config.yaml> <output_dir with *_dge_*.tsv> <figure_dir>")
CFG = yaml.safe_load(open(sys.argv[1])); SP, OUT = sys.argv[2], sys.argv[3]; os.makedirs(OUT, exist_ok=True)
DISPLAY_REMAP = CFG["figures"].get("display_names") or {}
PANEL_LABEL = CFG["figures"].get("panel_label") or "panel gene"   # legend wording for the panel genes
def contrasts_for(unit_cfg, which):
    """The comparisons of a unit, in the same order as supp_contrasts.R writes them."""
    lv = [l for l in unit_cfg["levels"] if l not in (unit_cfg.get("outside_stage_test") or [])]
    trans = [(lv[i+1], lv[i]) for i in range(len(lv)-1)]
    if unit_cfg.get("closing_contrast"): trans.append((lv[0], lv[-1]))
    extra = [tuple(c) for c in (unit_cfg.get("extra_contrasts") or [])]
    return trans if which == "transitions" else extra
GROUPS = [(g["unit"], contrasts_for(CFG["units"][g["unit"]], g["contrasts"]), g["contrasts"] == "transitions", g["tag"], g["title"],
           g.get("filename") or f"Fig{g['tag']}_{g['unit']}_dge")
          for g in CFG["figures"]["supplement"]]
for g in CFG["figures"]["supplement"]:                       # a misspelled key would otherwise be ignored
    bad = set(g) - {"unit", "contrasts", "tag", "title", "filename"}
    if bad: sys.exit(f"figures.supplement entry {g.get('tag')}: unknown keys {sorted(bad)}")
_names = [grp[-1] for grp in GROUPS]                            # two entries with one filename would overwrite each other
if len(set(_names)) != len(_names): sys.exit(f"figures.supplement: duplicate filenames {_names}")
STARS = {}
CUR = None; RECORD_STARS = False
for f in fm.findSystemFonts():                      # register Liberation Sans wherever the system keeps it
    if os.path.basename(f).startswith("LiberationSans"):
        try: fm.fontManager.addfont(f)
        except Exception: pass
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Liberation Sans","Arial","DejaVu Sans"],
    "pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none","axes.linewidth":0.8,
    "xtick.direction":"out","ytick.direction":"out","xtick.labelsize":7,"ytick.labelsize":7})

GRAY="#C7CDD6"; RED="#C0392B"; SIG_EDGE="#111111"; ARROW="→"

def load(tis, a, b):
    """Read the table of one comparison and normalize the columns the plots use."""
    f = os.path.join(SP, f"{tis}_dge_{a}_vs_{b}.tsv")
    if not os.path.exists(f):
        raise FileNotFoundError(f)
    d = pd.read_csv(f, sep="\t", low_memory=False)
    for c in ("padj","log2FoldChange","baseMean","lfcSE"):
        if c in d: d[c] = pd.to_numeric(d[c], errors="coerce")
    d["is_rtk"]  = d.get("is_rtk","").astype(str).str.upper().eq("TRUE")
    d["rtk_sig"] = d.get("rtk_sig","").astype(str).str.upper().eq("TRUE")
    if "gene" not in d: d["gene"] = ""
    d["gene"] = d["gene"].map(lambda g: DISPLAY_REMAP.get(g, g))   # apply display-name remaps for the label
    return d

def _split(d):
    """Separate the background features from the panel genes and precompute the plot coordinates."""
    rk = d[d.is_rtk].copy()
    rk["l2"]  = rk.log2FoldChange.fillna(0.0)
    rk["nlp"] = -np.log10(np.clip(rk.padj.fillna(1.0), 1e-300, None))
    rk["lm"]  = np.log10(np.clip(rk.baseMean, 1, None))
    rk["sig"] = rk["rtk_sig"] if "rtk_sig" in rk else (rk.padj.notna() & (rk.padj < 0.05))  # per-contrast sig flag
    return d[~d.is_rtk], rk

def _label(ax, sig, xcol, ycol):
    """label significant panel genes, nudging vertically to reduce overlap."""
    if sig.empty: return
    yr = ax.get_ylim(); dx = (ax.get_xlim()[1]-ax.get_xlim()[0])*0.012; dy = (yr[1]-yr[0])*0.045
    used = []
    for _, r in sig.sort_values(ycol, ascending=False).iterrows():
        x, y = r[xcol], r[ycol]; oy = dy
        while any(abs((y+oy)-u) < dy*0.9 for u in used): oy += dy*0.9
        used.append(y+oy)
        ax.annotate(r.gene, (x,y), (x+dx, y+oy), fontsize=5.6, color=RED, fontstyle="italic",
                    ha="left", va="bottom", zorder=6,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color=RED, shrinkA=0, shrinkB=1.5))
        if RECORD_STARS:   # record each drawn star once, from the volcano panel only (MA is guarded off)
            STARS.setdefault(CUR[0], []).append([CUR[1], r.gene])

def volcano(ax, d):
    """One volcano panel: background in gray, panel genes in red, significant panel genes starred."""
    bg, rk = _split(d)
    bg_nlp = -np.log10(np.clip(bg.padj.fillna(1),1e-300,None))
    ax.scatter(bg.log2FoldChange, bg_nlp, s=3, c=GRAY, alpha=0.5, edgecolor="none", rasterized=True)
    ns, sg = rk[~rk.sig], rk[rk.sig]
    ax.scatter(ns.l2, ns.nlp, s=13, c=RED, alpha=0.85, edgecolor="none", zorder=3)
    ax.scatter(sg.l2, sg.nlp, s=48, c=RED, marker="*", edgecolor=SIG_EDGE, lw=0.4, zorder=4)
    ax.axhline(-np.log10(0.05), ls="--", lw=0.6, c="#999"); ax.axvline(0, lw=0.5, c="#ccc")
    # fit both axes to each panel's own data so the cloud fills the box (not a flat band in an empty
    # frame), while clipping only the extreme tail and always keeping the panel genes in view
    rk_max = float(rk.nlp.max()) if len(rk) else 0.0
    ytop = min(float(np.nanpercentile(np.concatenate([bg_nlp.values, rk.nlp.values]), 99)), 25.0)
    ytop = max(ytop, rk_max*1.2, 4.0)
    allx = np.abs(np.concatenate([bg.log2FoldChange.values, rk.l2.values]))
    xr = max(float(np.nanpercentile(allx, 99)), (float(abs(rk.l2).max())*1.15 if len(rk) else 0), 3.0)
    ax.set_ylim(-0.02*ytop, ytop*1.06); ax.set_xlim(-xr, xr)
    _label(ax, sg, "l2", "nlp")
    ax.set_xlabel("log$_2$ fold change", fontsize=8)
    ax.set_ylabel("−log$_{10}$ adj. $p$", fontsize=8)
    ax.spines[["top","right"]].set_visible(False)

def ma(ax, d):
    """One MA panel: the same points as the volcano, against mean expression instead of significance."""
    bg, rk = _split(d)
    ax.scatter(np.log10(np.clip(bg.baseMean,1,None)), bg.log2FoldChange,
               s=3, c=GRAY, alpha=0.5, edgecolor="none", rasterized=True)
    ns, sg = rk[~rk.sig], rk[rk.sig]
    ax.scatter(ns.lm, ns.l2, s=13, c=RED, alpha=0.85, edgecolor="none", zorder=3)
    ax.scatter(sg.lm, sg.l2, s=48, c=RED, marker="*", edgecolor=SIG_EDGE, lw=0.4, zorder=4)
    ax.axhline(0, lw=0.6, c="#999")
    _label(ax, sg, "lm", "l2")
    ax.set_xlabel("log$_{10}$ mean expression", fontsize=8)
    ax.set_ylabel("log$_2$ fold change", fontsize=8)
    ax.spines[["top","right"]].set_visible(False)

def dge_figure(tis, pairs, is_transition, title, fname):
    """Assemble one supplementary figure: one row per comparison, volcano left, MA right."""
    global CUR, RECORD_STARS
    nr = len(pairs)
    fig, axes = plt.subplots(nr, 2, figsize=(6.8, 2.35*nr+0.7), squeeze=False)
    STARS.setdefault(fname, [])   # ensure a figure with no starred gene still gets an (empty) entry
    for i,(a,b) in enumerate(pairs):
        d = load(tis, a, b)
        CUR = (fname, f"{a}_vs_{b}")
        RECORD_STARS = True;  volcano(axes[i,0], d)   # record stars from the volcano panel
        RECORD_STARS = False; ma(axes[i,1], d)        # MA panel draws the same stars but must not double-count
        lab = f"{b} {ARROW} {a}" if is_transition else f"{a} vs {b}"
        axes[i,0].text(-0.36, 0.5, lab, transform=axes[i,0].transAxes,
                       rotation=90, ha="center", va="center", fontsize=9, fontweight="bold")
    axes[0,0].set_title("Volcano", fontsize=10, fontweight="bold")
    axes[0,1].set_title("MA", fontsize=10, fontweight="bold")
    for j,ax in enumerate(axes.reshape(-1)):
        ax.text(-0.22, 1.07, chr(65+j), transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")
    leg = [Line2D([0],[0], marker='o', ls='none', markerfacecolor=GRAY, markeredgecolor='none', markersize=5, label='all transcripts'),
           Line2D([0],[0], marker='o', ls='none', markerfacecolor=RED, markeredgecolor='none', markersize=6, label=PANEL_LABEL + ' of interest'),
           Line2D([0],[0], marker='*', ls='none', markerfacecolor=RED, markeredgecolor=SIG_EDGE, markersize=12, label=PANEL_LABEL + ' significant ($p_{adj}$<0.05)')]
    fig.legend(handles=leg, loc="upper center", ncol=3, frameon=False, fontsize=7.5, bbox_to_anchor=(0.5, 1.002))
    fig.suptitle(title, y=1.014, fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0.03,0,1,0.974])
    fig.savefig(os.path.join(OUT, fname + ".png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(OUT, fname + ".pdf"), dpi=600, bbox_inches="tight")  # vector annotations; rasterized bg at 600 dpi
    plt.close(fig)
    return nr

if __name__ == "__main__":
    for (tis, pairs, is_tr, tag, title, fname) in GROUPS:
        n = dge_figure(tis, pairs, is_tr, title, fname)
        print(f"  {fname}: {n} rows  ({title})")
    json.dump({k: sorted(set(map(tuple, v))) for k, v in STARS.items()}, open(os.path.join(OUT, f"{CFG['name']}_stars.json"), "w"))
    print(f"done -> {OUT}")
