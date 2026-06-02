#!/usr/bin/env python3
"""Build the boundary-weighted case set, COMPUTE gold, and assert it matches intent.

Each case: a measured neutral monoisotopic mass + a candidate formula (and, for
isotope cases, measured M+1/M+2 relative %). Verdict = pass (retain) / reject.

Design: weighted toward (a) valid-but-extreme formulas a naive model wrongly
REJECTS, and (b) ratio/probability/TMS/isotope decoys a naive model wrongly PASSES.
A few easy controls anchor the set.
"""
import pandas as pd
from golden_rules import parse_formula, mono_mass, check_rules, theoretical_isotopes

# (case_id, formula, derivatization, meas_M+1%, meas_M+2%, intended, note)
# meas values None -> isotope channel not provided for that case.
CASES = [
    # --- easy controls ---
    ("T01", "C6H12O6",   "none", None, None, "pass",   "glucose, textbook"),
    ("T02", "C9H11NO2",  "none", None, None, "pass",   "phenylalanine"),
    ("T03", "C2H8N",     "none", None, None, "reject", "RDBE<0 (R2)"),
    ("T04", "CH6",       "none", None, None, "reject", "H/C=6 (R4)"),

    # --- valid-but-extreme PASS (naive model over-rejects) ---
    ("T05", "C8H10N4O2", "none", None, None, "pass",   "caffeine; 4 N looks suspicious but N/C=0.5"),
    ("T06", "C5H5N5",    "none", None, None, "pass",   "adenine; N/C=1.0 <=1.3"),
    ("T07", "C6H8O7",    "none", None, None, "pass",   "citric acid; O/C=1.17, near cap but valid"),
    ("T08", "C2H6O",     "none", None, None, "pass",   "ethanol; H/C=3.0 boundary"),
    ("T09", "C8H14O2S2", "none", None, None, "pass",   "alpha-lipoic acid; 2 S, S/C=0.25"),
    ("T10", "C10H2",     "none", None, None, "pass",   "decadiyne-ish; H/C=0.2 lower boundary"),

    # --- ratio-violating REJECT decoys (naive model over-accepts) ---
    ("T11", "C2H8N4",    "none", None, None, "reject", "N/C=2.0 (R5)"),
    ("T12", "C2H4O4",    "none", None, None, "reject", "O/C=2.0 (R5)"),
    ("T13", "C2H6S3",    "none", None, None, "reject", "S/C=1.5 (R5)"),
    ("T14", "C3H9P3",    "none", None, None, "reject", "P/C=1.0 (R5)"),
    ("T15", "C20H2",     "none", None, None, "reject", "H/C=0.1 (R4)"),

    # --- valence / SENIOR REJECT (subtle) ---
    ("T16", "C2H5O",     "none", None, None, "reject", "non-integer RDBE -> radical, not neutral (R2)"),
    ("T17", "C3H10",     "none", None, None, "reject", "H/C=3.33 >3.1 AND RDBE<0 (R2/R4)"),

    # --- Rule 6 combined-probability decoys ---
    ("T18", "C18H30N4O4P4S4", "none", None, None, "reject", "P=4,S=4 NOPS combo (R6) though per-ratios pass"),
    ("T19", "C8H12N12",  "none", None, None, "reject", "N/C=1.5 (R5)"),

    # --- Rule 1 element-max ---
    ("T20", "C3H9N25",   "none", None, None, "reject", "25 N > 20 cap (R1) + N/C huge"),

    # --- isotope (Rule 3) consistency ---
    ("T21", "C7H7NO2",   "none", 7.8, 33.0, "reject", "M+2=33% implies 1 Cl, formula has none (R3)"),
    ("T22", "C6H4Cl2",   "none", 6.6, 64.0, "pass",   "M+2=64% consistent with 2 Cl"),
    ("T23", "C6H12O6",   "none", 22.0, 1.4, "reject", "M+1=22% implies ~20 C, formula has 6 (R3)"),
    ("T24", "C9H9BrO",   "none", 9.7, 97.0, "pass",   "M+2=97% consistent with 1 Br"),
    ("T25", "C10H14N2",  "none", 11.5, 0.8, "pass",   "nicotine; M+1 consistent (~11% for 10 C + 2 N)"),

    # --- Rule 7: TMS derivatization (GC-MS) ---
    ("T26", "C9H22O3Si2",  "tms",  None, None, "pass",   "lactic acid 2TMS; Si valid under R7"),
    ("T27", "C12H30O4Si3", "tms",  None, None, "pass",   "glyceric acid 3TMS; high Si ok under R7"),
    ("T28", "C9H22O3Si2",  "none", None, None, "reject", "same formula but NOT flagged TMS -> Si unexpected (R7 context)"),

    # --- more controls / fillers ---
    ("T29", "C16H34O",   "none", None, None, "pass",   "hexadecanol; H/C=2.1"),
    ("T30", "C44H10",    "none", None, None, "reject", "C44 > 39 cap (R1) + H/C=0.23 low"),
]


def isotope_consistent(counts, m1, m2):
    """True if measured M+1/M+2 (%) are within tolerance of theoretical."""
    t1, t2 = theoretical_isotopes(counts)
    def ok(meas, theo):
        return abs(meas - theo) <= max(2.0, 0.30 * theo)
    return ok(m1, t1) and ok(m2, t2)


def compute_verdict(formula, deriv, m1, m2):
    counts = parse_formula(formula)
    passes, failed = check_rules(counts, deriv)
    reasons = list(failed)
    # Rule 3 isotope channel
    if m1 is not None and m2 is not None:
        if not isotope_consistent(counts, m1, m2):
            passes = False
            reasons.append("R3:isotope_mismatch")
    # Rule 7 context: Si present but not declared TMS -> unexpected silicon
    if counts.get("Si", 0) > 0 and deriv.lower() != "tms":
        passes = False
        reasons.append("R7:undeclared_Si")
    return ("pass" if passes else "reject"), reasons


rows = []
mismatches = []
for cid, f, deriv, m1, m2, intended, note in CASES:
    counts = parse_formula(f)
    verdict, reasons = compute_verdict(f, deriv, m1, m2)
    if verdict != intended:
        mismatches.append((cid, f, intended, verdict, reasons, note))
    rows.append({
        "case_id": cid, "query_mass": round(mono_mass(counts), 4),
        "candidate_formula": f, "derivatization": deriv,
        "meas_mplus1_pct": m1, "meas_mplus2_pct": m2,
        "gold": verdict, "computed_reasons": ";".join(reasons), "note": note,
    })

df = pd.DataFrame(rows)
print(df[["case_id", "candidate_formula", "derivatization", "gold", "computed_reasons"]].to_string(index=False))
print(f"\nPASS: {(df.gold=='pass').sum()}  REJECT: {(df.gold=='reject').sum()}  total: {len(df)}")
if mismatches:
    print("\n!!! INTENT/COMPUTED MISMATCHES (fix before shipping):")
    for m in mismatches:
        print("  ", m)
else:
    print("\nAll computed verdicts match intent. ✓")

df.to_csv("cases_full.csv", index=False)
