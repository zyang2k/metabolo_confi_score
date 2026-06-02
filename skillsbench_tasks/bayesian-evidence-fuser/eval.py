#!/usr/bin/env python3
"""Standalone evaluator for bayesian-evidence-fuser.

Scores the agent's /root/ranked_posteriors.json against the gold posteriors
computed under the lab QC SOP (m/z-zoned Gaussian sigma + MS/MS 0.50 floor).

Score = fraction of checks passed, where the checks mirror the pytest suite:
  5 per-candidate posterior matches (tol 1e-4) + normalization + ranking = 7
plus the 2 format gates handled implicitly (structural failure -> score 0).

Usage:
    python eval.py <ranked_posteriors.json>
"""
import json
import sys

TOL = 1e-4

GOLD = {
    "L-Phenylalanine": 0.3548799017,
    "Citric_acid":     0.3274002778,
    "Adenosine":       0.1988559215,
    "LysoPC_16:0":     0.1143642591,
    "D-Glucose":       0.0044996399,
}
GOLD_ORDER = [
    "L-Phenylalanine", "Citric_acid", "Adenosine", "LysoPC_16:0", "D-Glucose",
]


def evaluate(path: str) -> dict:
    try:
        with open(path) as fh:
            data = json.load(fh)
    except Exception as e:
        return {"error": str(e), "score": 0.0}
    if not isinstance(data, list) or len(data) != 5:
        return {"error": "output must be a JSON array of 5 objects", "score": 0.0}
    try:
        parsed = [(int(o["rank"]), str(o["Candidate_Name"]), float(o["Posterior_Probability"]))
                  for o in data]
    except (KeyError, TypeError, ValueError) as e:
        return {"error": f"malformed object: {e}", "score": 0.0}

    post = {n: p for _, n, p in parsed}
    if set(post) != set(GOLD):
        return {"error": "candidate names do not match input set", "score": 0.0}

    total = len(GOLD) + 2  # 5 values + sum + ranking
    correct = sum(abs(post[n] - GOLD[n]) <= TOL for n in GOLD)
    value_correct = correct

    sum_ok = abs(sum(post.values()) - 1.0) <= TOL
    correct += int(sum_ok)

    by_rank = sorted(parsed, key=lambda t: t[0])
    order = [n for _, n, _ in by_rank]
    ranks = [r for r, _, _ in by_rank]
    rank_ok = (ranks == [1, 2, 3, 4, 5] and order == GOLD_ORDER)
    correct += int(rank_ok)

    return {
        "score": round(correct / total, 4),
        "correct": correct, "total": total,
        "values": f"{value_correct}/{len(GOLD)}",
        "normalized": sum_ok, "ranking": rank_ok,
    }


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python eval.py <ranked_posteriors.json>")
        sys.exit(1)
    r = evaluate(sys.argv[1])
    if "error" in r:
        print(f"ERROR: {r['error']}")
        print(f"Score: {r.get('score', 0.0):.4f}")
        sys.exit(1)
    print(f"Score: {r['score']:.4f}  ({r['correct']}/{r['total']})")
    print(f"  posterior values: {r['values']}")
    print(f"  normalized:       {r['normalized']}")
    print(f"  ranking:          {r['ranking']}")
