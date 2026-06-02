#!/usr/bin/env python3
"""Fused task probe: deconvolution skeleton + EMPIRICAL taxonomy crux.

Each feature carries an upstream adduct_annotation string (messy, multi-source).
Agent must: group co-eluting features -> infer neutral mass (PLUMBING, derivable),
then triage each annotation ok/isf/dubious per the lab taxonomy and decide a
per-compound keep/flag gate (CRUX, non-derivable).

Crux is concentrated in counter-intuitive cases:
  - look-alike DUBIOUS: raw-formula brackets ([M+C2H3O2]- = acetate, [M+COOH]- =
    formate) are chemically identical to OK shorthand but lab-flagged dubious.
  - weird OK: M-H1, M+TFA-H, [2M-H+2i]- look wrong but are lab-validated.
Gate flips on exactly these: a flag-bin's only ions are look-alike-ok (model keeps);
a keep-bin's only ions look weird (model flags).
"""
import pandas as pd

PROTON = 1.007276; E = 0.000549
H2O = 18.010565; CO2 = 43.989830
CL = 34.968853 + E; FORMATE = 44.998203; ACETATE = 59.013851
TFA = 112.985586; NEUTRON = 1.008665

def mz_of(M, ion):
    return {
        "MH":      M - PROTON,
        "acetate": M + ACETATE,
        "formate": M + FORMATE,
        "Cl":      M + CL,
        "TFA":     M + TFA,
        "dimer2i": 2*M - PROTON + 2*NEUTRON,
        "MH-H2O":  M - PROTON - H2O,
        "MH-CO2":  M - PROTON - CO2,
        "catdub":  M - PROTON - 2.0,        # [Cat-2H]- placeholder, dubious data
        "catisf":  M - PROTON - 176.0110,   # loss of C9H4O4
    }[ion]

# (compound, M, rt, gate, [(annotation, ion, category, is_isf, is_crux)])
FAMILIES = [
    ("citric acid", 192.027003, 9.21, "keep", [
        ("[M-H]-",        "MH",      "ok",      False, False),
        ("[M+CH3COO]-",   "acetate", "ok",      False, True),   # shorthand acetate = ok
        ("[M-H-CO2]-",    "MH-CO2",  "isf",     True,  False),
    ]),
    ("malic acid", 134.021524, 8.47, "keep", [
        ("M-H",           "MH",      "ok",      False, False),
        ("[M+C2H3O2]-",   "acetate", "dubious", False, True),   # raw-formula acetate = dubious
        ("[M-H-H2O]-",    "MH-H2O",  "isf",     True,  False),
    ]),
    ("glucose", 180.063388, 10.12, "keep", [
        ("[M-H]-",        "MH",      "ok",      False, False),
        ("[M+COOH]-",     "formate", "dubious", False, True),   # raw-formula formate = dubious
        ("M+Cl",          "Cl",      "ok",      False, False),
    ]),
    ("phenylalanine", 165.078979, 5.22, "keep", [
        ("M-H1",          "MH",      "ok",      False, True),   # weird historical notation = ok
        ("M+TFA-H",       "TFA",     "ok",      False, True),   # HILIC buffer adduct = ok
    ]),
    ("lactic acid", 90.031694, 6.51, "keep", [
        ("[M-H]-",        "MH",      "ok",      False, False),
        ("[2M-H+2i]-",    "dimer2i", "ok",      False, True),   # isotope-labeled dimer = ok
    ]),
    ("caffeine", 194.080376, 12.03, "keep", [
        ("[M-H]-",        "MH",      "ok",      False, False),
        ("[Cat-2H]-",     "catdub",  "dubious", False, True),   # no organic formula = dubious
        ("[Cat-2H-C9H4O4]-","catisf","isf",     True,  False),
    ]),
    # ---- FLAG bins: only look-alike-ok ions -> model wrongly keeps ----
    ("succinic acid", 118.026609, 7.88, "flag", [
        ("[M+C2H3O2]-",   "acetate", "dubious", False, True),
        ("[M+COOH]-",     "formate", "dubious", False, True),
        ("[M-H-H2O]-",    "MH-H2O",  "isf",     True,  False),
    ]),
    ("fumaric acid", 116.010959, 7.05, "flag", [
        ("[M+COOH]-",     "formate", "dubious", False, True),
        ("[M-H-H2O]-",    "MH-H2O",  "isf",     True,  False),
    ]),
]

INTERFERENTS = [("noise1", 250.1812, 4.10), ("noise2", 311.0827, 13.40),
                ("noise3", 88.9934, 3.05)]

rows, gold = [], []
fid = 1
for comp, M, rt, gate, feats in FAMILIES:
    for j, (ann, ion, cat, isf, crux) in enumerate(feats):
        rows.append({"feature_id": f"F{fid:02d}", "mz": round(mz_of(M, ion), 4),
                     "rt": round(rt + (j-1)*0.012, 3), "intensity": 100000//(j+1),
                     "adduct_annotation": ann})
        gold.append({"feature_id": f"F{fid:02d}", "compound": comp,
                     "neutral_mass": round(M, 4), "adduct_category": cat,
                     "is_isf": isf, "gate": gate, "is_crux": crux})
        fid += 1
for name, mz, rt in INTERFERENTS:
    rows.append({"feature_id": f"F{fid:02d}", "mz": mz, "rt": rt,
                 "intensity": 50000, "adduct_annotation": "unknown"})
    gold.append({"feature_id": f"F{fid:02d}", "compound": "(interferent)",
                 "neutral_mass": None, "adduct_category": "unassigned",
                 "is_isf": False, "gate": "n/a", "is_crux": False})
    fid += 1

pd.DataFrame(rows).sort_values("rt").to_csv("fused_input.csv", index=False)
pd.DataFrame(gold).to_csv("fused_gold.csv", index=False)
g = pd.DataFrame(gold)
print(g[g.adduct_category!="unassigned"].groupby("adduct_category").size())
print(f"\ncrux features: {g.is_crux.sum()} | gates: keep={sum(1 for f in FAMILIES if f[3]=='keep')} flag={sum(1 for f in FAMILIES if f[3]=='flag')}")
print(f"total features: {len(rows)}")
