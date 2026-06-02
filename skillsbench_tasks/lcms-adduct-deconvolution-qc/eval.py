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

GOLD_CRUX = {
    "F08": "dubious", "F11": "dubious", "F14": "dubious", "F17": "dubious",
    "F20": "ok",      "F22": "dubious", "F24": "dubious", "F25": "dubious",
    "F27": "dubious", "F28": "dubious", "F30": "dubious", "F31": "dubious",
    "F33": "dubious", "F34": "dubious",
}
GOLD_GATE = {
    "F01": "keep", "F04": "keep", "F07": "keep", "F10": "keep", "F13": "keep",
    "F16": "keep", "F19": "keep", "F21": "keep",
    "F24": "flag", "F27": "flag", "F30": "flag", "F33": "flag",
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
