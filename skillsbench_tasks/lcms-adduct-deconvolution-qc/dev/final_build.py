#!/usr/bin/env python3
"""Final dataset for lcms-adduct-deconvolution-qc.

Negative-mode HILIC feature list. Agent deconvolves (group co-eluting features ->
neutral mass), triages each upstream adduct_annotation ok/isf/dubious per the lab
taxonomy, and sets a per-compound keep/flag gate.

Crux (skill-critical, is_crux=True) concentrates the lab's non-derivable rule:
RAW-FORMULA-BRACKET adducts ([M+C2H3O2]- = acetate, [M+COOH]- = formate,
[M+C2F3O2]- = TFA) are chemically identical to OK shorthand but lab-flagged DUBIOUS.
A no-skill model reads them as valid adducts -> "ok", and keeps flag-bins built only
from them. Plus a few weird-ok / isotope-dimer cases for variety.

Verifier scores ONLY: adduct_category on crux features + keep/flag gate.
Deconvolution is required to answer the gate but is not itself scored.
"""
import pandas as pd

PROTON = 1.007276; E = 0.000549
H2O = 18.010565; CO2 = 43.989830
CL = 34.968853 + E; FORMATE = 44.998203; ACETATE = 59.013851
TFA = 112.985586; NEUTRON = 1.008665

def mz_of(M, ion):
    return {
        "MH":       M - PROTON,
        "acetate":  M + ACETATE,
        "formate":  M + FORMATE,
        "Cl":       M + CL,
        "TFA":      M + TFA,
        "dimer2i":  2*M - PROTON + 2*NEUTRON,
        "MH-H2O":   M - PROTON - H2O,
        "MH-CO2":   M - PROTON - CO2,
        "catdub":   M - PROTON - 2.0,
        "catisf":   M - PROTON - 176.0110,
    }[ion]

# templates of (annotation, ion, category, is_isf, is_crux)
def keep_std(extra_shorthand="acetate"):
    sh = {"acetate": "[M+CH3COO]-", "formate": "[M+HCOO]-"}[extra_shorthand]
    return [("[M-H]-", "MH", "ok", False, False),
            (sh, extra_shorthand, "ok", False, False),
            ("[M-H-H2O]-", "MH-H2O", "isf", True, False)]

def keep_lookalike(raw="acetate"):
    raws = {"acetate": "[M+C2H3O2]-", "formate": "[M+COOH]-", "TFA": "[M+C2F3O2]-"}
    return [("[M-H]-", "MH", "ok", False, False),
            (raws[raw], raw, "dubious", False, True),
            ("[M-H-CO2]-", "MH-CO2", "isf", True, False)]

def flag_lookalike(raw_a="acetate", raw_b="formate"):
    raws = {"acetate": "[M+C2H3O2]-", "formate": "[M+COOH]-", "TFA": "[M+C2F3O2]-"}
    return [(raws[raw_a], raw_a, "dubious", False, True),
            (raws[raw_b], raw_b, "dubious", False, True),
            ("[M-H-H2O]-", "MH-H2O", "isf", True, False)]

def keep_isotope():
    return [("[M-H]-", "MH", "ok", False, False),
            ("[2M-H+2i]-", "dimer2i", "ok", False, True)]

def keep_cat():
    return [("[M-H]-", "MH", "ok", False, False),
            ("[Cat-2H]-", "catdub", "dubious", False, True),
            ("[Cat-2H-C9H4O4]-", "catisf", "isf", True, False)]

# (compound, M, rt, gate, feats)
FAMILIES = [
    ("citric acid",     192.027003, 9.21, "keep", keep_std("acetate")),
    ("phenylalanine",   165.078979, 5.22, "keep", keep_std("acetate")),
    ("malic acid",      134.021524, 8.47, "keep", keep_lookalike("acetate")),
    ("glucose",         180.063388, 10.12,"keep", keep_lookalike("formate")),
    ("glutaric acid",   132.042259, 8.01, "keep", keep_lookalike("acetate")),
    ("tartaric acid",   150.016438, 9.03, "keep", keep_lookalike("TFA")),
    ("lactic acid",      90.031694, 6.51, "keep", keep_isotope()),
    ("gluconic acid",   196.058303, 11.20,"keep", keep_lookalike("formate")),
    ("succinic acid",   118.026609, 7.88, "flag", flag_lookalike("acetate", "formate")),
    ("fumaric acid",    116.010959, 7.05, "flag", flag_lookalike("formate", "TFA")),
    ("2-oxoglutarate",  146.021524, 8.90, "flag", flag_lookalike("acetate", "formate")),
    ("aconitic acid",   174.016438, 9.55, "flag", flag_lookalike("acetate", "TFA")),
]
INTERFERENTS = [("noise1", 250.1812, 4.10), ("noise2", 311.0827, 13.40),
                ("noise3", 288.9934, 3.05), ("noise4", 402.1521, 2.20)]

rows, gold = [], []
fid = 1
for comp, M, rt, gate, feats in FAMILIES:
    for j, (ann, ion, cat, isf, crux) in enumerate(feats):
        mz = mz_of(M, ion)
        rows.append({"feature_id": f"F{fid:02d}", "mz": round(mz, 4),
                     "rt": round(rt + (j-1)*0.012, 3), "intensity": 120000//(j+1),
                     "adduct_annotation": ann})
        gold.append({"feature_id": f"F{fid:02d}", "compound": comp,
                     "neutral_mass": round(M, 4), "adduct_category": cat,
                     "is_isf": isf, "gate": gate, "is_crux": crux})
        fid += 1
for name, mz, rt in INTERFERENTS:
    rows.append({"feature_id": f"F{fid:02d}", "mz": mz, "rt": rt,
                 "intensity": 60000, "adduct_annotation": "unknown"})
    gold.append({"feature_id": f"F{fid:02d}", "compound": "(interferent)",
                 "neutral_mass": None, "adduct_category": "unassigned",
                 "is_isf": False, "gate": "n/a", "is_crux": False})
    fid += 1

inp = pd.DataFrame(rows).sort_values("rt").reset_index(drop=True)
pd.DataFrame(gold).to_csv("final_gold.csv", index=False)
inp.to_csv("final_input.csv", index=False)
g = pd.DataFrame(gold)
print("category counts:\n", g[g.adduct_category!="unassigned"].groupby("adduct_category").size().to_string())
print(f"\ncrux features: {int(g.is_crux.sum())}  | gates: keep={sum(1 for f in FAMILIES if f[3]=='keep')} flag={sum(1 for f in FAMILIES if f[3]=='flag')}")
print(f"total features: {len(rows)}  (scored: {int(g.is_crux.sum())} crux + {len(FAMILIES)} gates = {int(g.is_crux.sum())+len(FAMILIES)} tests)")
