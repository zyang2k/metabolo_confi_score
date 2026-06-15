"""
Validate the oliver-curate skill against the validated subset.

Validated subset = bins under names where Oliver actually had to choose
(name has ≥1 TP bin AND ≥1 FP bin — see project_oliver_curation_method).

For each bin in the subset, pick the top-entropy candidate hit (the
pipeline's natural pick), apply the skill rules, and compare verdict
against spectrum_label.
"""

import re
import pandas as pd

CONTAM_PATTERNS = [
    r"\b(hexa|hepta|octa|nona|deca)ethylene glycol\b",
    r"\bpolyethylene glycol\b",
    r"\b(tri|tetra)methylamine\b",
    r"\btriethylamine\b",
    r"\boctadecyl[- ]?trimethylammonium\b",
    r"\bmethylguanidine\b",
    r"\bpipenzolate\b",
    r"\birganox\b",
    r"\bbht\b",
    r"\b2,6-di-tert-butyl",
    r"\blauryl diethanolamide\b",
    r"\b2-\(octadecylamino\)ethanol\b",
    r"\bsulfooxy\b",
    r"\bdiethanolamine\b",
]
CONTAM_RE = re.compile("|".join(CONTAM_PATTERNS), re.IGNORECASE)

# Adducts treated as monoisotopic (Δm extreme veto applies)
MONOISOTOPIC_ADDUCTS = {"[M+H]+", "[M-H]-", "[M+Na]+", "[M+K]+", "[M+NH4]+", "M+H", "M-H", "M+Na", "M+K", "M+NH4"}

N_ACETYL_RE = re.compile(r"\bn[- ]?acetyl", re.IGNORECASE)


def is_monoisotopic(adduct):
    if not isinstance(adduct, str):
        return False
    return adduct.strip() in MONOISOTOPIC_ADDUCTS


def evaluate_row(row):
    """Return (verdict, signals_list)."""
    signals = []  # list of (name, sign, weight_tier, weight_val)

    def add(name, sign, tier):
        weight_map = {"veto": 100, "high": 3, "med": 2, "low": 1}
        signals.append((name, sign, tier, weight_map[tier]))

    # --- Anchoring ---
    n_adducts = row.get("n_candidate_adducts", 1) or 1
    has_ok_adduct = bool(row.get("compound_has_ok_adduct", 0))
    if n_adducts >= 2 and has_ok_adduct:
        add("multi_adduct_anchor", "+", "veto")

    if bool(row.get("hit_isf_no_ok", 0)):
        esim_chk = row.get("entropy_similarity") or 0
        dmda_chk = row.get("delta_mda")
        drt_chk = row.get("signed_delta_rt")
        # Demote ISF veto when all other evidence is clean (validated 2026-05-21).
        clean_other = (
            esim_chk >= 0.85
            and pd.notna(dmda_chk) and abs(dmda_chk) < 3
            and pd.notna(drt_chk) and abs(drt_chk) < 30
        )
        add("isf_no_anchor", "-", "high" if clean_other else "veto")

    # --- Mass error ---
    dmda = row.get("delta_mda")
    adduct = row.get("adduct", "")
    if pd.notna(dmda):
        adm = abs(dmda)
        # Validated 2026-05-21: TP rows have median |Δm| ~0.58 mDa for ALL adduct classes.
        # The mono vs non-mono distinction in Δm tolerance was artificial. Apply Oliver's
        # ~10 mDa threshold universally; only modulate the positive (tight) by mono.
        if adm > 10:
            add("delta_mda_extreme", "-", "veto")
        elif 5 <= adm <= 10:
            add("delta_mda_gray", "-", "low")
        elif adm < 2:
            tier = "med" if is_monoisotopic(adduct) else "low"
            add("delta_mda_tight", "+", tier)

    # --- RT error ---
    drt = row.get("signed_delta_rt")
    name = row.get("name", "") or ""
    kong_off = bool(N_ACETYL_RE.search(name))
    if pd.notna(drt):
        adrt = abs(drt)
        if adrt < 10:
            add("delta_rt_tight", "+", "med")
        elif 20 <= adrt <= 50:
            add("delta_rt_gray", "-", "low")
        elif adrt > 50:
            if kong_off:
                add("kong_rt_known_off", "neutral", "med")
                add("delta_rt_gray", "-", "low")  # downgrade veto to gray
            else:
                add("delta_rt_extreme", "-", "veto")

    # --- MS2 ---
    esim = row.get("entropy_similarity", 0)
    rcos = row.get("reverse_cosine", 0)
    if pd.notna(esim) and esim >= 0.75:
        add("forward_sim_strong", "+", "med")
    if pd.notna(esim) and esim < 0.50:
        add("forward_sim_weak", "-", "low")
    if pd.notna(rcos) and rcos >= 0.80:
        add("reverse_cosine_strong", "+", "high")
    if pd.notna(esim) and pd.notna(rcos) and esim < 0.50 and rcos < 0.60:
        add("forward_reverse_both_weak", "-", "veto")
    sentropy = row.get("spectral_entropy")
    if pd.notna(sentropy) and sentropy < 1.0:
        add("spectral_entropy_sparse", "-", "low")
    sgap = row.get("sim_gap")
    if pd.notna(sgap) and sgap < 0.10:
        add("sim_gap_narrow", "-", "med")

    # --- Adduct taxonomy ---
    cat = row.get("hit_adduct_cat")
    if cat == "ok":
        add("adduct_ok_cat", "+", "low")
    elif cat == "dubious":
        add("adduct_dubious_cat", "-", "med")

    # --- Library trust ---
    db = row.get("db", "")
    hit_ik = row.get("hit_ik14")
    if isinstance(db, str) and db.lower().startswith("nist"):
        add("db_nist_trusted", "+", "low")
    if (not isinstance(hit_ik, str) or not hit_ik.strip()) and isinstance(db, str) and db.lower() in {"gnps", "massbank", "mona"}:
        add("db_empty_ik14_untrusted", "-", "low")

    # --- Contamination ---
    if isinstance(name, str) and CONTAM_RE.search(name):
        # Spectrum is annotated here (validated set), so weight lower than for blanks
        add("contaminant_pattern", "-", "med")

    # --- Noise gate (Orbitrap Snorm) ---
    snorm = row.get("normalized_entropy")
    polarity = row.get("polarity", "")
    if pd.notna(snorm) and snorm >= 0.987:
        add("snorm_failed_noise_gate", "-", "veto")

    # --- Verdict logic ---
    vetos_pos = [s for s in signals if s[2] == "veto" and s[1] == "+"]
    vetos_neg = [s for s in signals if s[2] == "veto" and s[1] == "-"]

    if vetos_neg:
        verdict = "yy"
    elif vetos_pos:
        neg_high = [s for s in signals if s[2] == "high" and s[1] == "-"]
        verdict = "dubious" if neg_high else "keep"
    else:
        net = sum((+1 if s[1] == "+" else (-1 if s[1] == "-" else 0)) * s[3]
                  for s in signals)
        if net >= 3:
            verdict = "keep"
        elif net >= 1:
            verdict = "probable_ok"
        elif net >= -1:
            verdict = "dubious"
        else:
            # net <= -2
            any_pos = any(s[1] == "+" for s in signals)
            verdict = "dubious" if any_pos else "yy"

    # --- Top 3 signals by |weight| ---
    signals_sorted = sorted(signals, key=lambda s: (-s[3], 0 if s[1] == "+" else 1))
    top3 = signals_sorted[:3]
    signals_str = " | ".join(f"{s[0]}:{s[1]}" for s in top3)
    return verdict, signals_str


def main():
    df = pd.read_csv("data/feature_table_v2.csv")
    bins = df.drop_duplicates("wiki_id")[["wiki_id", "anno_name_lower", "spectrum_label"]]
    name_labels = bins.groupby("anno_name_lower")["spectrum_label"].agg(set)
    validated_names = name_labels[name_labels.apply(lambda s: "TP" in s and "FP" in s)].index
    validated_bins = bins[bins["anno_name_lower"].isin(validated_names)]
    v_rows = df[df["wiki_id"].isin(validated_bins["wiki_id"])]

    # Top-entropy candidate per bin (the pipeline's natural pick)
    top = (v_rows.sort_values("entropy_similarity", ascending=False)
                 .drop_duplicates("wiki_id", keep="first")
                 .copy())

    # Score
    verdicts, signals = [], []
    for _, r in top.iterrows():
        v, s = evaluate_row(r)
        verdicts.append(v)
        signals.append(s)
    top["oliver_verdict"] = verdicts
    top["oliver_signals"] = signals

    # Confusion table
    print(f"=== Oliver skill validation on validated subset ===")
    print(f"Validated names: {len(validated_names)}")
    print(f"Validated bins:  {len(top)}")
    print()
    print("Confusion (rows=verdict, cols=spectrum_label):")
    confusion = pd.crosstab(top["oliver_verdict"], top["spectrum_label"], margins=True)
    confusion = confusion.reindex(["keep", "probable_ok", "dubious", "yy", "All"])
    print(confusion)
    print()

    # Reduce to binary call: keep / probable_ok → predict TP; yy → predict FP; dubious abstains
    def to_pred(v):
        if v in ("keep", "probable_ok"):
            return "TP"
        if v == "yy":
            return "FP"
        return "abstain"

    top["pred"] = top["oliver_verdict"].apply(to_pred)
    decided = top[top["pred"] != "abstain"]
    if len(decided):
        n_correct = (decided["pred"] == decided["spectrum_label"]).sum()
        print(f"Decided rows: {len(decided)}/{len(top)} ({100*len(decided)/len(top):.0f}%)")
        print(f"  Accuracy among decided: {n_correct}/{len(decided)} = {100*n_correct/len(decided):.1f}%")
        # Per-class
        for lbl in ["TP", "FP"]:
            sub = decided[decided["spectrum_label"] == lbl]
            if len(sub):
                acc = (sub["pred"] == lbl).sum() / len(sub)
                print(f"  {lbl}: {(sub['pred']==lbl).sum()}/{len(sub)} = {100*acc:.1f}%")
    print()
    abstain = top[top["pred"] == "abstain"]
    if len(abstain):
        print(f"Abstain (dubious) rows: {len(abstain)}")
        print(f"  TP within abstain: {(abstain['spectrum_label']=='TP').sum()}")
        print(f"  FP within abstain: {(abstain['spectrum_label']=='FP').sum()}")

    # Misclassified examples
    print("\n=== Worst misses ===")
    misses_keep_on_fp = top[(top["oliver_verdict"].isin(["keep", "probable_ok"])) &
                            (top["spectrum_label"] == "FP")]
    misses_yy_on_tp = top[(top["oliver_verdict"] == "yy") &
                          (top["spectrum_label"] == "TP")]
    print(f"\n'keep'/'probable_ok' on FP bins (skill kept what Oliver rejected): {len(misses_keep_on_fp)}")
    if len(misses_keep_on_fp):
        cols = ["wiki_id", "anno_name_lower", "name", "adduct", "entropy_similarity",
                "delta_mda", "signed_delta_rt", "oliver_verdict", "oliver_signals"]
        print(misses_keep_on_fp[cols].head(10).to_string(index=False))
    print(f"\n'yy' on TP bins (skill rejected what Oliver kept): {len(misses_yy_on_tp)}")
    if len(misses_yy_on_tp):
        cols = ["wiki_id", "anno_name_lower", "name", "adduct", "entropy_similarity",
                "delta_mda", "signed_delta_rt", "oliver_verdict", "oliver_signals"]
        print(misses_yy_on_tp[cols].head(10).to_string(index=False))

    top.to_csv("data/oliver_skill_validation.csv", index=False)
    print("\nFull table → data/oliver_skill_validation.csv")


if __name__ == "__main__":
    main()
