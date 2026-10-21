#!/usr/bin/env python3
"""Summary table of the stage test for every panel gene (the supplementary significance table).

Usage: table_panel.py <config.yaml> <directory with the *_<panel_tag>_letters.tsv files> <output directory>

Reads the letters table of every unit written by dge_panel.R and writes one row per panel gene
and unit with: the gene, its family, the unit, the p-value of the likelihood-ratio test, the
targeted false discovery rate (adjusted across the panel), the transcriptome-wide false discovery rate,
and whether the gene is called regulated (targeted rate below alpha). Two files come out: a
tab-separated table (.tsv) for machines, and a Markdown table (.md) with the numbers in
scientific notation for the manuscript. Rows follow either the order of the panel file
(tables.order: panel) or family order and then display name (tables.order: family_name).
BASE_LOCAL names the project directory under which the panel file is found; by default it is the
directory two levels above the directory that holds the config (a config at
<root>/pipeline/configs/x.yaml resolves to <root>). Set BASE_LOCAL=. when running from the project root, as
in the README."""
import math, os, re, sys
import pandas as pd, yaml

SUP = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")   # superscript digits for the Markdown exponents
def sci(x):
    """Format a p-value for the Markdown table, e.g. 3.1 × 10⁻¹¹. NA stays the text NA."""
    if pd.isna(x): return "NA"
    if x == 0:     return "0"
    e = math.floor(math.log10(abs(x))); m = round(x / 10**e, 1)
    if m >= 10: m /= 10; e += 1
    return f"{m:.1f} × 10{str(e).translate(SUP)}"

def main(cfg_path, data, out):
    cfg = yaml.safe_load(open(cfg_path)); os.makedirs(out, exist_ok=True)
    tag = cfg.get("panel_tag", "panel")
    base = os.environ.get("BASE_LOCAL", os.path.dirname(os.path.dirname(os.path.abspath(cfg_path))) + "/..")
    panel = pd.read_csv(os.path.join(base, cfg["panel"]), sep="\t") if os.path.exists(os.path.join(base, cfg["panel"])) else None
    units = list(cfg["units"].keys()); ucol = cfg["unit_col"]; ulabel = ucol.capitalize()
    disp = cfg.get("figures", {}).get("display_names", {}) or {}
    fam_order = cfg["figures"]["family_order"]
    # Collect the test statistics of every gene in every unit from the letters tables (one row per
    # gene is enough, because the statistics are repeated on each stage row).
    rows_by = {}; families = {}
    for u in units:
        d = pd.read_csv(f"{data}/{u}_{tag}_letters.tsv", sep="\t").drop_duplicates("gene")
        for _, r in d.iterrows():
            rows_by[(r.gene, u)] = (float(r.lrt_p), float(r.q_targeted), float(r.q_transcriptomewide))
    genes_in_order = list(dict.fromkeys(g for g, _ in rows_by))
    if panel is not None:
        families = dict(zip(panel.gene, panel.family)); panel_order = list(dict.fromkeys(panel.gene))
    else:
        families = {g: re.sub(r"[0-9].*$", "", re.sub(r"^[A-Za-z]+-", "", disp.get(g, g))) for g in genes_in_order}; panel_order = genes_in_order
    # Row order, then the rows themselves. The regulated call is recomputed here from q_targeted so
    # that the table cannot disagree with the letters tables.
    if cfg["tables"]["order"] == "panel": genes = [g for g in panel_order if (g, units[0]) in rows_by]
    else: genes = sorted(genes_in_order, key=lambda g: (fam_order.index(families[g]), disp.get(g, g)))
    rows = []
    for g in genes:
        for u in units:
            p, qt, qg = rows_by[(g, u)]
            rows.append(dict(gene=disp.get(g, g), family=families[g], **{ucol: u}, LRT_p=p, q_targeted=qt, q_transcriptomewide=qg,
                             molt_regulated=bool(qt < cfg["alpha"]) if not pd.isna(qt) else False))
    tbl = pd.DataFrame(rows); name = cfg["tables"]["name"]
    tsv = tbl.copy(); tsv["molt_regulated"] = tsv["molt_regulated"].map(lambda b: "TRUE" if b else "FALSE")
    tsv.to_csv(f"{out}/{name}.tsv", sep="\t", index=False, na_rep="NA")
    CHECK, DASH = "✓", "—"
    lines = [f"| Gene | Family | {ulabel} | LRT *p* | *q* (targeted) | *q* (transcriptome-wide) | Molt-regulated |",
             "| :--- | :--- | :--- | ---: | ---: | ---: | :---: |"]
    for _, r in tbl.iterrows():
        lines.append(f"| {r.gene} | {r.family} | {r[ucol]} | {sci(r.LRT_p)} | {sci(r.q_targeted)} | {sci(r.q_transcriptomewide)} | {CHECK if r.molt_regulated else DASH} |")
    open(f"{out}/{name}.md", "w").write("\n".join(lines) + "\n")
    print(f"{name}: {len(tbl)} rows, {int(tbl.molt_regulated.sum())} molt-regulated, {int(tbl.LRT_p.isna().sum())} NA LRT p -> {out}")

if __name__ == "__main__":
    if len(sys.argv) != 4: sys.exit(__doc__)
    main(*sys.argv[1:])
