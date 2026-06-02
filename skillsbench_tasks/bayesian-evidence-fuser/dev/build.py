#!/usr/bin/env python3
"""Build / regenerate the bayesian-evidence-fuser task data and gold values.

Run from anywhere:  python dev/build.py

Writes the input CSV to environment/data/input.csv and prints the
gold posteriors (paste into tests/test_outputs.py GOLD and eval.py GOLD if the
input ever changes). Also prints the no-skill "textbook" answer so the
discrimination margin is visible at a glance.

The two non-derivable cruxes encoded in the gold:
  Channel B  m/z-zoned Gaussian sigma:  mz<250 -> 4.0; 250..450 -> 2.0; mz>450 -> 3.0
  Channel C  MS/MS reliability floor:   cosine < 0.50 -> 0.01, else raw
"""
import csv
import math
import os

# name, prior, observed_mz, mass_error_ppm, msms_similarity
ROWS = [
    ("D-Glucose",        0.30, 179.0561, 3.0, 0.48),
    ("Citric_acid",      0.20, 191.0197, 1.0, 0.85),
    ("L-Phenylalanine",  0.20, 164.0717, 0.5, 0.90),
    ("LysoPC_16:0",      0.15, 480.3096, 4.0, 0.70),
    ("Adenosine",        0.15, 266.0894, 2.0, 0.55),
]

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(HERE, "..", "environment", "data", "input.csv")


def npdf(x, s, mu=0.0):
    return (1.0 / (s * math.sqrt(2 * math.pi))) * math.exp(-((x - mu) ** 2) / (2 * s * s))


def sigma_zone(mz):
    if mz < 250:
        return 4.0
    if mz <= 450:
        return 2.0
    return 3.0


def ms2(s):
    return 0.01 if s < 0.50 else s


def fuse(sigma_fn, ms2_fn):
    u = [(n, p * npdf(ppm, sigma_fn(mz)) * ms2_fn(s)) for n, p, mz, ppm, s in ROWS]
    z = sum(v for _, v in u)
    return sorted(((n, v / z) for n, v in u), key=lambda t: -t[1])


def write_csv():
    with open(CSV_PATH, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Candidate_Name", "Prior_Probability", "Observed_mz",
                    "Mass_Error_ppm", "MSMS_Similarity_Score"])
        for n, p, mz, ppm, s in ROWS:
            w.writerow([n, p, mz, ppm, s])
    print(f"wrote {CSV_PATH}")


if __name__ == "__main__":
    write_csv()
    print("\nGOLD (lab SOP — paste into tests/eval):")
    for n, p in fuse(sigma_zone, ms2):
        print(f'    "{n}": {p:.10f},')
    print("\nNo-skill textbook (constant sigma=2.0, raw MS/MS) — should DIFFER:")
    for n, p in fuse(lambda mz: 2.0, lambda s: s):
        print(f"    {n:18s} {p:.10f}")
