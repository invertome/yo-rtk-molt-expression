#!/usr/bin/env python3
"""Compare two output directories of the pipeline, table by table and cell by cell.

Usage: compare_outputs.py <ref_dir> <new_dir> --units U1,U2 --name <run> [--tag panel] [--tol 1e-6] [--allow-missing]

Purpose: answer "did this rerun give the same results as an existing set of tables?" with a yes or
no that does not depend on harmless differences. Text columns (letters, calls, sample names) must
match exactly. Numeric columns may differ by floating-point noise, which different machines and
library versions produce in the last decimal places, but not by more than the tolerance. P-values
get one extra rule, because very small p-values carry their information in the exponent: two
p-values must also agree within 0.01 on the log10 scale. The result is one line per expected
table (PASS, FAIL with the offending cells, SKIPPED, EXTRA or N/A) and a one-line summary; the
exit status is 0 only when every expected table passed. The details follow.

--tag is the config's panel_tag and selects the per-unit panel tables.
Canonicalizes column names (sample->sampleid, rep_ids->id, experiment|tissue->unit,
first column of *_deseq2_LRT_all.tsv and *_dge_*.tsv -> feature), then compares per file, row by key:
- strings exactly (only the boolean spellings TRUE/True/true and FALSE/False/false are unified; NA must
  be NA on both sides);
- numerics with identical NA patterns; finite pairs within |ref - new| <= tol * max(1, |ref|); non-finite
  values (Inf, -Inf) must match exactly, including sign;
- p-value-like columns (pvalue, padj, lrt_p, q_targeted, q_transcriptomewide) within tol relative to
  max(1, |ref|) and within 0.01 on the log10 scale, everywhere (values are floored at 1e-300, so two exact zeros are
  equal). Reason: cross-machine floating-point noise moves extreme tail p-values by more than 1e-6 relative,
  while a changed tail computation (e.g. 1 - pchisq) moves them by orders of magnitude; only the second
  must fail;
- duplicated keys on either side fail the file.
Every expected file is printed with one of: PASS, FAIL (itemized), SKIPPED (present in ref, missing in new),
EXTRA (present only in new, not compared), N/A (absent from both sides, not compared); every other *.tsv in
either directory is listed as UNCHECKED (never gating). Exit 2 on a usage error (missing directory, empty
unit or run name, tol outside [0, 0.1]); 1 if any file fails, if nothing was compared, or (without
--allow-missing) if a reference file is skipped; 0 otherwise. The tool never writes into either directory.
"""
import argparse, glob, math, os, re, sys
import numpy as np, pandas as pd

RENAME = {"sample": "sampleid", "rep_ids": "id", "experiment": "unit", "tissue": "unit"}
SPEC = {  # filename pattern -> (key columns, compared columns)
    "{u}_{t}_letters.tsv":      (["gene", "stage"], ["letter", "lrt_p", "q_targeted", "q_transcriptomewide", "molt_regulated"]),
    "{u}_{t}_tpm_by_stage.tsv": (["gene", "stage"], ["id", "n", "mean_tpm", "sem_tpm"]),
    "{u}_deseq2_LRT_all.tsv":   (["feature"], ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]),
    "{u}_dge_*.tsv":            (["feature"], ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "is_rtk", "rtk_sig"]),
    "{u}_pca.tsv":              (["sampleid"], ["PC1", "PC2", "stage"]),
    "{n}_{t}_persample_tpm.tsv": (["gene", "sampleid"], ["stage", "tpm"]),
}
PVAL_COLS = {"pvalue", "padj", "lrt_p", "q_targeted", "q_transcriptomewide"}
PVAL_LOG10_TOL, PVAL_FLOOR = 0.01, 1e-300
BOOL_TOKENS = {"TRUE": "TRUE", "True": "TRUE", "true": "TRUE", "FALSE": "FALSE", "False": "FALSE", "false": "FALSE"}

def is_na(v):
    """True for any representation of a missing value."""
    return v is None or v is pd.NA or (isinstance(v, (float, np.floating)) and np.isnan(v))

def plain(v):
    """Render a cell for messages: TRUE/FALSE for booleans, NA for missing, Python scalars otherwise."""
    if isinstance(v, (bool, np.bool_)): return "TRUE" if v else "FALSE"
    if is_na(v): return "NA"
    if isinstance(v, np.integer): return int(v)
    if isinstance(v, np.floating): return float(v)
    return v

def key(i):
    return tuple(plain(v) for v in i) if isinstance(i, tuple) else plain(i)

def norm_token(v):
    """Value for exact string comparison: booleans unified to TRUE/FALSE, NA kept as None, all else verbatim."""
    if isinstance(v, (bool, np.bool_)): return "TRUE" if v else "FALSE"
    if is_na(v): return None
    s = str(v)
    return BOOL_TOKENS.get(s, s)

def side_missing(col, ref, new):
    a, b = col not in ref.columns, col not in new.columns
    return "both" if a and b else "ref" if a else "new" if b else None

def load(path, first_to_feature):
    """Read one table and bring its column names to the canonical spelling used by SPEC."""
    df = pd.read_csv(path, sep="\t", keep_default_na=True, na_values=["NA", ""], low_memory=False)
    if first_to_feature: df = df.rename(columns={df.columns[0]: "feature"})
    return df.rename(columns={k: v for k, v in RENAME.items() if k in df.columns and v not in df.columns})

def compare_numeric(a, b, c, tol, label):
    """Compare one numeric column: same missing values, same infinities, finite values within tol."""
    problems = []
    na_a, na_b = a.isna(), b.isna(); bad_na = na_a != na_b
    if bad_na.any():
        problems.append(f"{label}: {c}: NA pattern differs for {int(bad_na.sum())} rows (NA only in ref: {int((na_a & ~na_b).sum())}, NA only in new: {int((~na_a & na_b).sum())}), e.g. {[key(i) for i in a.index[bad_na][:3]]}")
    m = ~na_a & ~na_b
    x, y = a[m].astype(float), b[m].astype(float)
    fin = np.isfinite(x.to_numpy()) & np.isfinite(y.to_numpy())
    xn, yn = x[~fin], y[~fin]
    if not np.array_equal(xn.to_numpy(), yn.to_numpy()):
        badn = xn.to_numpy() != yn.to_numpy(); w = int(np.argmax(badn))
        problems.append(f"{label}: {c}: {int(badn.sum())} non-finite values differ (worst at {key(xn.index[w])}: ref={plain(xn.iloc[w])} new={plain(yn.iloc[w])})")
    x, y = x[fin], y[fin]
    diff = (x - y).abs().to_numpy(); lim = tol * np.maximum(1.0, x.abs().to_numpy())
    bad = diff > lim
    if bad.any():
        ratio = np.divide(diff, lim, out=np.full_like(diff, np.inf), where=lim > 0)
        w = int(np.argmax(np.where(bad, ratio, -1.0)))
        problems.append(f"{label}: {c}: {int(bad.sum())} values beyond tol (worst at {key(x.index[w])}: ref={plain(x.iloc[w])} new={plain(y.iloc[w])})")
    if c in PVAL_COLS:
        ld = np.abs(np.log10(np.maximum(x.to_numpy(), PVAL_FLOOR)) - np.log10(np.maximum(y.to_numpy(), PVAL_FLOOR)))
        badl = ld > PVAL_LOG10_TOL
        if badl.any():
            w = int(np.argmax(np.where(badl, ld, -1.0)))
            problems.append(f"{label}: {c}: {int(badl.sum())} p-value-like values differ on the log10 scale beyond {PVAL_LOG10_TOL} (worst at {key(x.index[w])}: ref={plain(x.iloc[w])} new={plain(y.iloc[w])})")
    return problems

def compare_strings(a, b, c, label):
    """Compare one text column exactly (after unifying the spellings of TRUE and FALSE)."""
    sa, sb = a.map(norm_token), b.map(norm_token)
    bad = ~((sa == sb) | (sa.isna() & sb.isna()))
    if not bad.any(): return []
    pos = np.flatnonzero(bad.to_numpy())[:3]
    return [f"{label}: {c}: {int(bad.sum())} values differ, e.g. {[(key(a.index[p]), plain(a.iloc[p]), plain(b.iloc[p])) for p in pos]}"]

def compare_frame(ref, new, keys, cols, tol, label):
    """Compare two tables row by key: duplicated keys, differing row sets, then each listed column."""
    for k in keys:
        side = side_missing(k, ref, new)
        if side: return [f"{label}: key column {k!r} missing in {side}"]
    ref = ref.set_index(keys).sort_index(); new = new.set_index(keys).sort_index()
    problems = []
    for side, df in (("ref", ref), ("new", new)):
        dup = df.index[df.index.duplicated()].unique()
        if len(dup): problems.append(f"{label}: {len(dup)} duplicated keys in {side}, e.g. {[key(i) for i in dup[:3]]}")
    if problems: return problems
    if not ref.index.equals(new.index):
        only_ref = ref.index.difference(new.index); only_new = new.index.difference(ref.index)
        problems.append(f"{label}: row sets differ (only in ref: {len(only_ref)}, e.g. {[key(i) for i in only_ref[:3]]}; only in new: {len(only_new)}, e.g. {[key(i) for i in only_new[:3]]})")
        common = ref.index.intersection(new.index); ref = ref.loc[common]; new = new.loc[common]
    for c in cols:
        side = side_missing(c, ref, new)
        if side: problems.append(f"{label}: column {c!r} missing in {side}"); continue
        a, b = ref[c], new[c]
        numeric = all(pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s) for s in (a, b))
        problems += compare_numeric(a, b, c, tol, label) if numeric else compare_strings(a, b, c, label)
    return problems

def die(msg):
    print(msg, file=sys.stderr); sys.exit(2)

def main():
    # Command-line arguments, then one pass over SPEC: every expected filename is built from the unit
    # names, the run name and the panel tag, and classified by where it exists before being compared.
    ap = argparse.ArgumentParser(); ap.add_argument("ref"); ap.add_argument("new")
    ap.add_argument("--units", required=True); ap.add_argument("--name", required=True); ap.add_argument("--tag", default="panel")
    ap.add_argument("--tol", type=float, default=1e-6); ap.add_argument("--allow-missing", action="store_true")
    a = ap.parse_args()
    for side, d in (("ref", a.ref), ("new", a.new)):
        if not os.path.isdir(d): die(f"{side} directory not found: {d}")
    units = [u.strip() for u in a.units.split(",")]
    if "" in units: die(f"--units must be a comma-separated list of non-empty unit names, got {a.units!r}")
    name = a.name.strip()
    if name == "": die(f"--name must be a non-empty run name, got {a.name!r}")
    if not re.fullmatch(r"[A-Za-z0-9_]+", a.tag): die(f"--tag must be one word of letters, digits or underscores, got {a.tag!r}")
    if not (math.isfinite(a.tol) and 0 <= a.tol <= 0.1): die(f"--tol must be finite and within [0, 0.1], got {a.tol!r}")
    n_pass = n_fail = n_skip = n_extra = n_na = 0; report = []; expected = set()
    for pat, (keys, cols) in SPEC.items():
        pt = pat.replace("{t}", a.tag)
        names = []
        if "{u}" in pt:
            for u in units:
                p = pt.replace("{u}", u)
                if "*" in p:
                    found = sorted({os.path.basename(x) for d in (a.ref, a.new) for x in glob.glob(os.path.join(d, p))})
                    if not found: n_na += 1; report.append(f"N/A {p}: no matches on either side")
                    names += found
                else:
                    names.append(p)
        else:
            names = [pt.replace("{n}", name)]
        expected.update(names)
        for fn in names:
            rp, np_ = os.path.join(a.ref, fn), os.path.join(a.new, fn)
            in_ref, in_new = os.path.exists(rp), os.path.exists(np_)
            if not in_ref and not in_new:
                n_na += 1; report.append(f"N/A {fn}: absent from both sides (not compared)"); continue
            if in_ref and not in_new:
                n_skip += 1; report.append(f"SKIPPED {fn}: missing in new"); continue
            if in_new and not in_ref:
                n_extra += 1; report.append(f"EXTRA {fn}: present only in new (not compared)"); continue
            f2f = fn.endswith("_deseq2_LRT_all.tsv") or "_dge_" in fn
            probs = compare_frame(load(rp, f2f), load(np_, f2f), keys, cols, a.tol, fn)
            if probs: n_fail += 1; report += ["FAIL " + p for p in probs]
            else: n_pass += 1; report.append(f"PASS {fn}")
    present = {}
    for side, d in (("ref", a.ref), ("new", a.new)):
        for f in os.listdir(d):
            if f.endswith(".tsv") and os.path.isfile(os.path.join(d, f)): present.setdefault(f, set()).add(side)
    unchecked = sorted(f for f in present if f not in expected)
    for f in unchecked: report.append(f"UNCHECKED {f} ({'both' if len(present[f]) == 2 else next(iter(present[f]))})")
    print("\n".join(report))
    nothing = n_pass + n_fail == 0
    if nothing: print("FAIL: no file was compared (check the directories, --units, --name and --tag)")
    print(f"files compared: pass={n_pass} fail={n_fail} skipped={n_skip} extra={n_extra} n/a={n_na} unchecked={len(unchecked)} (tol={a.tol})")
    sys.exit(1 if nothing or n_fail else (0 if n_skip == 0 or a.allow_missing else 1))

if __name__ == "__main__":
    main()
