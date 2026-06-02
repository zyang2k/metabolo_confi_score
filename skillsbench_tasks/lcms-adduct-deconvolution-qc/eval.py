#!/usr/bin/env python3
"""Evaluator for lcms-adduct-deconvolution-qc.

Scores ONLY the crux: adduct_category on the 14 skill-critical features + the 12
keep/flag gates (graded at each compound's principal feature). Deconvolution is
required to answer the gate but is not directly scored.

Usage:
    python eval.py <output_csv> [gold_csv]
"""
import sys
import pandas as pd

# Crux is dominated by the false-positive-suppression direction: obscure-but-
# lab-validated adducts (isotope-labeled standards, historical notation) that a
# no-skill model over-flags as "dubious" -> ok (F13/F14/F16/F18/F20/F22). A
# minority of raw-formula-bracket adducts -> dubious tests the under-flag rule.
GOLD_CRUX = {
    "F08": "dubious", "F11": "dubious", "F13": "ok",      "F14": "ok",
    "F16": "ok",      "F18": "ok",      "F20": "ok",      "F22": "ok",
    "F24": "dubious", "F25": "dubious", "F27": "dubious", "F28": "dubious",
    "F30": "dubious", "F31": "dubious",
}
# Gate graded at each compound's principal (highest-intensity) feature.
# Durable over-flag keeps (only intact ion is a validated isotope/historical
# form): F13, F16, F18, F20, F22 -> a no-skill model wrongly flags them.
# Under-flag flags (only raw-formula dubious): F24, F27, F30.
GOLD_GATE = {
    "F01": "keep", "F04": "keep", "F07": "keep", "F10": "keep",
    "F13": "keep", "F16": "keep", "F18": "keep", "F20": "keep", "F22": "keep",
    "F24": "flag", "F27": "flag", "F30": "flag",
}


def evaluate(output_path: str) -> dict:
    try:
        pred = pd.read_csv(output_path, dtype=str, keep_default_na=False)
    except Exception as e:
        return {"error": str(e), "score": 0.0}
    for col in ("feature_id", "adduct_category", "gate"):
        if col not in pred.columns:
            return {"error": f"missing column {col}", "score": 0.0}

    cat = dict(zip(pred["feature_id"], pred["adduct_category"].astype(str).str.strip().str.lower()))
    gate = dict(zip(pred["feature_id"], pred["gate"].astype(str).str.strip().str.lower()))

    crux_correct = sum(cat.get(f) == v for f, v in GOLD_CRUX.items())
    gate_correct = sum(gate.get(f) == v for f, v in GOLD_GATE.items())
    total = len(GOLD_CRUX) + len(GOLD_GATE)
    correct = crux_correct + gate_correct
    return {
        "score": round(correct / total, 4),
        "correct": correct, "total": total,
        "crux": f"{crux_correct}/{len(GOLD_CRUX)}",
        "gate": f"{gate_correct}/{len(GOLD_GATE)}",
    }


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python eval.py <output_csv>")
        sys.exit(1)
    r = evaluate(sys.argv[1])
    if "error" in r:
        print(f"ERROR: {r['error']}")
        print(f"Score: {r.get('score', 0.0):.4f}")
        sys.exit(1)
    print(f"Score: {r['score']:.4f}  ({r['correct']}/{r['total']})")
    print(f"  crux category: {r['crux']}")
    print(f"  bin gate:      {r['gate']}")
