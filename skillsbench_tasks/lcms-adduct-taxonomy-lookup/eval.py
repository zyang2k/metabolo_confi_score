#!/usr/bin/env python3
"""Evaluator for lcms-adduct-taxonomy-lookup task.

Usage:
    python eval.py <output_csv> [gold_csv]

Reads agent output and gold answers, reports accuracy.
"""
import sys
import pandas as pd


def evaluate(output_path: str, gold_path: str = "data/gold.csv") -> dict:
    gold = pd.read_csv(gold_path)
    try:
        pred = pd.read_csv(output_path)
    except Exception as e:
        return {"error": str(e), "score": 0.0, "correct": 0, "total": len(gold)}

    required_cols = {"case_id", "category"}
    if not required_cols.issubset(pred.columns):
        return {
            "error": f"output.csv must have columns {required_cols}, got {set(pred.columns)}",
            "score": 0.0,
            "correct": 0,
            "total": len(gold),
        }

    pred["category"] = pred["category"].str.strip().str.lower()
    valid = {"ok", "isf", "dubious"}
    invalid = set(pred["category"]) - valid
    if invalid:
        return {
            "error": f"Invalid category values: {invalid}. Must be 'ok', 'isf', or 'dubious'.",
            "score": 0.0,
            "correct": 0,
            "total": len(gold),
        }

    merged = gold.merge(pred[["case_id", "category"]], on="case_id", suffixes=("_gold", "_pred"))
    if len(merged) < len(gold):
        missing = set(gold["case_id"]) - set(pred["case_id"])
        return {
            "error": f"Missing predictions for case_ids: {missing}",
            "score": 0.0,
            "correct": 0,
            "total": len(gold),
        }

    correct = (merged["category_gold"] == merged["category_pred"]).sum()
    total = len(merged)
    score = correct / total

    by_class = {}
    for cat in ["ok", "isf", "dubious"]:
        subset = merged[merged["category_gold"] == cat]
        c = (subset["category_gold"] == subset["category_pred"]).sum()
        by_class[cat] = {"correct": int(c), "total": len(subset), "accuracy": c / len(subset) if len(subset) > 0 else 0.0}

    return {
        "score": round(score, 4),
        "correct": int(correct),
        "total": total,
        "by_class": by_class,
    }


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python eval.py <output_csv> [gold_csv]")
        sys.exit(1)

    output_path = sys.argv[1]
    gold_path = sys.argv[2] if len(sys.argv) > 2 else "data/gold.csv"

    result = evaluate(output_path, gold_path)

    if "error" in result:
        print(f"ERROR: {result['error']}")
        print(f"Score: {result.get('score', 0.0):.4f}")
        sys.exit(1)

    print(f"Score: {result['score']:.4f}  ({result['correct']}/{result['total']} correct)")
    print()
    for cat, stats in result.get("by_class", {}).items():
        print(f"  {cat}: {stats['correct']}/{stats['total']} = {stats['accuracy']:.1%}")
