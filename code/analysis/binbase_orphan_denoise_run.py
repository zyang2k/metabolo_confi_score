"""binbase_orphan_denoise_run.py — reverse/containment denoise on real LCB studies + validation.

PRIMARY PIPELINE (--unconfirmed): of the generated-but-unaccepted UNCONFIRMED candidate
bins, which are an ISF / adduct / 13C-isotope of an existing CONFIRMED compound (->
relational artifact, should NOT be promoted) vs genuinely novel? These candidates passed
all of CARROT's QC gates, so the artifact-vs-novel call is the real promote/don't-promote
decision. Method-level: each candidate matched vs confirmed bins by retention_index.

    source .venv_bench/bin/activate
    python code/analysis/binbase_orphan_denoise_run.py --unconfirmed
    # inputs: data/lcb_hilicneg_unconfirmed.csv (SQL stmt 7) + data/lcb_hilicneg_bins.csv (stmt 3)

LEGACY / OPTIONAL MODES (kept for the record, not the pipeline):
  --validate-bins   confirmed-vs-confirmed library audit (separate QC deliverable).
  --from-csv        per-sample INVALID_TARGET orphan denoise. NOTE: INVALID_TARGET is mostly
                    run-level ISTD-coverage rejection (real compounds from QC-failed runs),
                    NOT spectral noise -> wrong population to denoise. Optionally + --same-sample.
  (no flag)         connect to carrot-prod directly (needs lab-network/VPN DNS).
"""
import os, re, sys
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "code", "analysis"))
import isf_orphan_denoise as isf

DATA = os.path.join(ROOT, "data")
BINS_CSV = os.path.join(DATA, "lcb_hilicneg_bins.csv")
ORPH_CSV = os.path.join(DATA, "lcb_hilicneg_orphans.csv")
DET_CSV  = os.path.join(DATA, "lcb_sample_detections.csv")
UNCONF_CSV = os.path.join(DATA, "lcb_hilicneg_unconfirmed.csv")
RI_WIN = 4.0
ORPHAN_TYPES = ("INVALID_TARGET", "UNCONFIRMED")

_argv = [a for a in sys.argv[1:] if not a.startswith("--")]
METHOD = _argv[0] if _argv else "5m hilic premier | orbitrap | beh amide | negative"


def parse_msms(s):
    out = []
    for tok in str(s).split():
        if ":" in tok:
            a, b = tok.split(":", 1)
            try:
                out.append([float(a), float(b)])
            except ValueError:
                pass
    return out or None


def peaks_dict(*dfs):
    pk = {}
    for df in dfs:
        for _, r in df.iterrows():
            p = parse_msms(r["msms"])
            if p:
                pk[r["wiki_id"]] = p
    return pk


def _rate(orphans, bins, peaks):
    res = isf.denoise(orphans[["wiki_id", "rt", "precursor_mz"]],
                      bins[["wiki_id", "rt", "precursor_mz"]], peaks, rt_win=RI_WIN)
    return res, float(res["explained"].fillna(False).mean())


def denoise_same_sample(orphans, bins, peaks, detected):
    """Restrict each orphan's candidate parents to bins detected in its own injection."""
    parts = []
    binids = set(bins["wiki_id"])
    for s, og in orphans.groupby("sample"):
        ids = detected.get(s, set()) & binids
        refs = bins[bins["wiki_id"].isin(ids)]
        if not len(refs):
            continue
        parts.append(isf.denoise(og[["wiki_id", "rt", "precursor_mz"]],
                                 refs[["wiki_id", "rt", "precursor_mz"]], peaks, rt_win=RI_WIN))
    return pd.concat(parts) if parts else pd.DataFrame(columns=["wiki_id", "explained", "parent", "relation", "containment"])


# ─────────────────── loaders ───────────────────
def from_csv():
    bins = pd.read_csv(BINS_CSV).rename(columns={"ri": "rt"})
    orphans = pd.read_csv(ORPH_CSV).rename(columns={"ri": "rt"})
    if "target_type" in orphans:
        orphans = orphans[orphans["target_type"].isin(ORPHAN_TYPES)].copy()
    print(f"loaded {len(bins)} bins, {len(orphans)} orphans "
          f"({orphans['sample'].nunique() if 'sample' in orphans else 1} injections)\n")
    return bins, orphans


def from_db():
    import psycopg2
    def creds():
        src = open(os.path.join(ROOT, "code", "r_eda", "get_binbase_data.R")).read()
        g = lambda k: re.findall(rf'{k}\s*=\s*"([^"]+)"', src)[-1]
        return dict(host=g("host"), dbname=g("dbname"), user=g("user"), password=g("password"), port=5432)
    con = psycopg2.connect(connect_timeout=15, **creds())
    bins = pd.read_sql("""SELECT id AS wiki_id, accurate_mass AS precursor_mz, retention_index AS rt,
             name, fragment_of, msms FROM compound
             WHERE method=%s AND target_type='CONFIRMED' AND msms IS NOT NULL AND msms<>''""", con, params=(METHOD,))
    orphans = pd.read_sql("""WITH samp AS (SELECT sample FROM compound
               WHERE method=%s AND target_type IN ('INVALID_TARGET','UNCONFIRMED') AND msms IS NOT NULL AND msms<>''
               GROUP BY sample ORDER BY count(*) DESC LIMIT 50)
             SELECT c.id AS wiki_id, c.sample, c.target_type, c.reason,
                    c.accurate_mass AS precursor_mz, c.retention_index AS rt, c.fragment_of, c.msms
             FROM compound c WHERE c.method=%s AND c.target_type IN ('INVALID_TARGET','UNCONFIRMED')
               AND c.msms IS NOT NULL AND c.msms<>'' AND c.sample IN (SELECT sample FROM samp)""",
             con, params=(METHOD, METHOD))
    con.close()
    return bins, orphans


# ─────────────────── validation #1: labeled concordance on CONFIRMED bins ───────────────────
def validate_bins():
    bins = pd.read_csv(BINS_CSV).rename(columns={"ri": "rt"})
    if "fragment_of" not in bins:
        print("bins.csv has no fragment_of column — re-export statement (3) with the new columns.")
        return
    n_lab = int(bins["fragment_of"].notna().sum())
    print(f"CONFIRMED bins: {len(bins)}   CARROT-labeled ISF (fragment_of set): {n_lab}")
    if n_lab == 0:
        print("\nNo CARROT ISF labels on confirmed bins for this method -> labeled concordance not possible here.")
        print("(Same as the orphans: CARROT's in-source-fragment annotation is absent in this method —")
        print(" i.e. the very gap this tool fills. Validation must lean on null controls + same-injection")
        print(" + curator spot-check, or a method where fragment_of IS populated.)")
        return
    peaks = peaks_dict(bins)
    res = isf.denoise(bins[["wiki_id", "rt", "precursor_mz"]],
                      bins[["wiki_id", "rt", "precursor_mz"]], peaks, rt_win=RI_WIN)
    bb = bins.merge(res, on="wiki_id", how="left"); bb["explained"] = bb["explained"].fillna(False)
    lab = bb["fragment_of"].notna()
    rec = bb.loc[lab, "explained"].mean()
    sub = bb[lab & bb["explained"]]
    agree = (sub["parent"].astype(float) == sub["fragment_of"].astype(float)).mean() if len(sub) else float("nan")
    print("\n=== validation #1: concordance vs CARROT fragment_of ===")
    print(f"  recall  (we flag a CARROT-labeled ISF):           {100*rec:.1f}%  ({int(bb.loc[lab,'explained'].sum())}/{n_lab})")
    print(f"  parent-id agreement (of those we flag):           {100*agree:.1f}%")
    print(f"  we also flag (CARROT did NOT label):              {int(bb.loc[~lab,'explained'].sum())} "
          f"({100*bb.loc[~lab,'explained'].mean():.1f}% of unlabeled) -> candidate ISF CARROT missed")
    print("  NB precision-vs-CARROT is a LOWER bound — CARROT under-labels ISF.")


# ─────────────────── UNCONFIRMED candidate-bin audit ───────────────────
def audit_unconfirmed():
    """Are the generated-but-unaccepted candidate bins actually ISF/adduct/isotope of an
    existing CONFIRMED compound (-> shouldn't be promoted), or genuinely novel?
    Method-level: match each UNCONFIRMED candidate vs CONFIRMED bins by retention_index.
    --tag <name> selects data/<name>_bins.csv + data/<name>_unconfirmed.csv (default lcb_hilicneg)."""
    tag = sys.argv[sys.argv.index("--tag") + 1] if "--tag" in sys.argv else "lcb_hilicneg"
    bins_csv = os.path.join(DATA, f"{tag}_bins.csv")
    if not os.path.exists(bins_csv) and os.path.exists(os.path.join(DATA, f"{tag}_confirmed.csv")):
        bins_csv = os.path.join(DATA, f"{tag}_confirmed.csv")   # accept _confirmed naming too
    unc_csv = os.path.join(DATA, f"{tag}_unconfirmed.csv")
    for p in (bins_csv, unc_csv):
        if not os.path.exists(p):
            print(f"Missing {p} — export it first (SQL stmt 3 for bins, stmt 7 for unconfirmed)."); return
    if "lipid" in tag:
        isf.NEUTRAL_LOSSES = {**isf.NEUTRAL_LOSSES, **isf.LIPID_LOSSES}
        print(f"[tag={tag}]  [lipid mode: +{len(isf.LIPID_LOSSES)} headgroup/fatty-acyl losses]")
    else:
        print(f"[tag={tag}]")
    bins = pd.read_csv(bins_csv).rename(columns={"ri": "rt"})
    unc = pd.read_csv(unc_csv).rename(columns={"ri": "rt"})
    print(f"UNCONFIRMED candidate bins: {len(unc)}   CONFIRMED reference: {len(bins)}\n")
    peaks = peaks_dict(bins, unc)
    res, rate = _rate(unc, bins, peaks)
    o = unc.merge(res, on="wiki_id", how="left"); o["explained"] = o["explained"].fillna(False)
    n = int(o.explained.sum())
    print("=" * 70)
    print(f"{n}/{len(o)} ({100*rate:.1f}%) UNCONFIRMED candidates are an ISF/adduct/isotope of a")
    print("  co-eluting CONFIRMED bin -> relational artifact, likely should NOT be promoted.")
    print(f"  residual {len(o)-n} ({100*(1-rate):.1f}%) -> genuine novel candidates.")
    print("=" * 70)
    o["reltype"] = o["relation"].fillna("").str.split(":").str[0]
    print("relation classes among flagged:")
    print("  " + o[o.explained]["reltype"].value_counts().to_string().replace("\n", "\n  "))

    rng = np.random.default_rng(0)
    u_mz = unc.copy(); u_mz["precursor_mz"] = rng.permutation(u_mz["precursor_mz"].values)
    u_rt = unc.copy(); u_rt["rt"] = rng.permutation(u_rt["rt"].values)
    _, r_mz = _rate(u_mz, bins, peaks); _, r_rt = _rate(u_rt, bins, peaks)
    save = (isf.NEUTRAL_LOSSES, isf.ADDUCT_DELTAS, isf.ISOTOPE)
    isf.NEUTRAL_LOSSES, isf.ADDUCT_DELTAS, isf.ISOTOPE = isf.DECOY_LOSSES, {}, {}
    _, r_dec = _rate(unc, bins, peaks)
    isf.NEUTRAL_LOSSES, isf.ADDUCT_DELTAS, isf.ISOTOPE = save
    print(f"\nNULL CONTROLS:  real {100*rate:.1f}%  |  random-Δm/z {100*r_mz:.1f}%  |  "
          f"RI-shuffle {100*r_rt:.1f}%  |  decoy-loss {100*r_dec:.1f}%")

    nm = bins[["wiki_id", "name"]].rename(columns={"wiki_id": "parent", "name": "pn"})
    e = o[o.explained].merge(nm, on="parent", how="left")
    e = e[~e["pn"].fillna("").str.startswith("unknown_")]
    print(f"\nflagged candidates linking to a NAMED confirmed bin: {len(e)}")
    print("examples (unknown candidate -> named confirmed parent it's a related ion of):")
    for _, r in e.sort_values("containment", ascending=False).head(12).iterrows():
        print(f"  cand m/z={r.precursor_mz:9.4f} RI={r.rt:6.1f}  {str(r.relation):14s} "
              f"cont={r.containment:.2f}  <- {str(r['pn'])[:36]}")
    out = os.path.join(DATA, f"{tag}_unconfirmed_denoise.csv")
    o.drop(columns=[c for c in ("msms",) if c in o]).to_csv(out, index=False)
    print(f"\nper-candidate table -> {out}")


# ─────────────────── main orphan analysis ───────────────────
def analyze(bins, orphans, detected=None):
    if not len(orphans) or not len(bins):
        print("Empty bins or orphans."); return
    peaks = peaks_dict(bins, orphans)
    res, rate = _rate(orphans, bins, peaks)
    o = orphans.merge(res, on="wiki_id", how="left"); o["explained"] = o["explained"].fillna(False)

    print("=" * 70)
    print(f"NOISE-REDUCTION (method-wide parents): {int(o.explained.sum())}/{len(o)} ({100*rate:.1f}%)")
    print("=" * 70)
    o["reltype"] = o["relation"].fillna("").str.split(":").str[0]
    print("relation classes among explained:")
    print("  " + o[o.explained]["reltype"].value_counts().to_string().replace("\n", "\n  "))
    if "sample" in o:
        per = o.groupby("sample")["explained"].mean() * 100
        print(f"\nper-injection explained-rate ({per.size}): median {per.median():.1f}%  "
              f"IQR [{per.quantile(.25):.1f},{per.quantile(.75):.1f}]%  range [{per.min():.1f},{per.max():.1f}]%")

    # null controls (method-wide)
    rng = np.random.default_rng(0)
    o_mz = orphans.copy(); o_mz["precursor_mz"] = rng.permutation(o_mz["precursor_mz"].values)
    o_rt = orphans.copy(); o_rt["rt"] = rng.permutation(o_rt["rt"].values)
    _, r_mz = _rate(o_mz, bins, peaks); _, r_rt = _rate(o_rt, bins, peaks)
    save = (isf.NEUTRAL_LOSSES, isf.ADDUCT_DELTAS, isf.ISOTOPE)
    isf.NEUTRAL_LOSSES, isf.ADDUCT_DELTAS, isf.ISOTOPE = isf.DECOY_LOSSES, {}, {}
    _, r_dec = _rate(orphans, bins, peaks)
    isf.NEUTRAL_LOSSES, isf.ADDUCT_DELTAS, isf.ISOTOPE = save
    print(f"\nNULL CONTROLS:  real {100*rate:.1f}%  |  random-Δm/z {100*r_mz:.1f}%  |  "
          f"RI-shuffle {100*r_rt:.1f}%  |  decoy-loss {100*r_dec:.1f}%")

    # validation #2: same-injection parents (compared apples-to-apples on the SAME
    # injections that have a detections list)
    if detected is not None:
        cov = [s for s in orphans["sample"].unique() if s in detected]
        osub = orphans[orphans["sample"].isin(cov)].copy()
        print("\n=== validation #2: same-injection parent (parent must be detected in the orphan's run) ===")
        print(f"  injections with detections: {len(cov)}/{orphans['sample'].nunique()}  "
              f"(n={len(osub)} orphans in those)")
        if len(osub):
            _, mw_sub = _rate(osub, bins, peaks)                    # method-wide on same subset
            ss = denoise_same_sample(osub, bins, peaks, detected)
            ss_rate = float(ss["explained"].fillna(False).mean()) if len(ss) else float("nan")
            osub_rt = osub.copy(); osub_rt["rt"] = rng.permutation(osub_rt["rt"].values)
            ss_shuf = denoise_same_sample(osub_rt, bins, peaks, detected)
            ss_shuf_rate = float(ss_shuf["explained"].fillna(False).mean()) if len(ss_shuf) else float("nan")
            ndet = float(np.mean([len(detected.get(s, set())) for s in cov]))
            print(f"  avg confirmed bins detected per injection: {ndet:.0f}  (of {len(bins)} method-wide)")
            print(f"  explained, method-wide parents:    {100*mw_sub:.1f}%")
            print(f"  explained, SAME-INJECTION parents: {100*ss_rate:.1f}%   <- defensible noise-reduction")
            print(f"  RI-shuffle background: method-wide {100*r_rt:.1f}%  ->  same-injection {100*ss_shuf_rate:.1f}%")
            print("  (lower same-injection background = the co-elution+presence link is real, not coincidental.)")

    if "reason" in o:
        o["reason_short"] = o["reason"].fillna("").str.extract(r"\(([A-Za-z]+)")[0].fillna("(none)")
        print("\nexplained-rate by BinBase invalidation reason:")
        xt = o.groupby("reason_short")["explained"].agg(["size", "mean"]).sort_values("size", ascending=False)
        for rs, row in xt.head(8).iterrows():
            print(f"  {rs:28s} n={int(row['size']):6d}  explained {100*row['mean']:.1f}%")

    out = os.path.join(DATA, "lcb_orphan_denoise.csv")
    o.drop(columns=[c for c in ("msms",) if c in o]).to_csv(out, index=False)
    print(f"\nper-orphan table -> {out}")


def load_detected():
    d = pd.read_csv(DET_CSV)
    m = {}
    for s, g in d.groupby("sample"):
        m[s] = set(g["compound_id"].astype("int64"))
    print(f"loaded same-injection detections: {len(m)} injections, "
          f"{int(d.shape[0])} (sample,compound) rows\n")
    return m


def main():
    if "--unconfirmed" in sys.argv:
        audit_unconfirmed(); return
    if "--validate-bins" in sys.argv:
        validate_bins(); return
    bins, orphans = from_csv() if "--from-csv" in sys.argv else from_db()
    detected = load_detected() if ("--same-sample" in sys.argv and os.path.exists(DET_CSV)) else None
    analyze(bins, orphans, detected)


if __name__ == "__main__":
    main()
