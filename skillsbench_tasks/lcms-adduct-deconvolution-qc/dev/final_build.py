#!/usr/bin/env python3
"""Final dataset for lcms-adduct-deconvolution-qc (rebalanced 2026-06-02).

Negative-mode HILIC feature list. The agent deconvolves (group co-eluting
features -> neutral mass), triages each upstream adduct_annotation
ok/isf/dubious per the lab taxonomy, and sets a per-compound keep/flag gate.

WHY THIS REBALANCE
------------------
The previous crux leaned entirely on "raw-formula-bracket -> dubious"
([M+C2H3O2]- = acetate, [M+COOH]- = formate, [M+C2F3O2]- = TFA). That points
in the GUESSABLE direction: a cautious frontier model already flags weird
notation on its own, so the no-skill baseline climbed and the with/without
delta collapsed.

The durable crux is the INVERSE direction -- false-positive suppression:
obscure-but-lab-validated adducts that a cautious model WRONGLY flags.
  - [2M-H+2i]-, [M-H+1i]-  : isotope-labeled internal-standard ions. This lab
                             spikes stable-isotope standards; "+Ni" = N heavy
                             isotopes. Non-derivable lab-method knowledge.
  - M-H1                   : historical deprotonation notation the legacy
                             pipeline still emits.
A no-skill model reads these as malformed/unusual -> "dubious", so it also
WRONGLY FLAGS any compound whose only intact-ion evidence is one of them.
This direction does NOT erode as models get stronger (more caution -> more
over-flagging), and the rule is a substantive, provenance-backed lab fact
rather than a punitive notation gotcha.

A minority of raw-formula-bracket cases is retained so the task still tests
both directions (under-flag as well as over-flag).

Verifier scores ONLY: adduct_category on crux features + keep/flag gate
(graded at each compound's principal/highest-intensity feature). Deconvolution
is required to answer the gate but is not itself scored.
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
        "TFA":      M + TFA,
        "dimer2i":  2*M - PROTON + 2*NEUTRON,
        "dimer4i":  2*M - PROTON + 4*NEUTRON,
        "mono1i":   M - PROTON + 1*NEUTRON,
        "MH-H2O":   M - PROTON - H2O,
        "MH-CO2":   M - PROTON - CO2,
    }[ion]

# Each template returns a list of (annotation, ion, category, is_isf, is_crux).
# The FIRST feature is the principal (highest intensity) -> where the gate is graded.

# --- trivial keep: has a plain [M-H]- ok; gate keeps regardless. crux=0 -----
def keep_std(shorthand="acetate"):
    sh = {"acetate": "[M+CH3COO]-", "formate": "[M+HCOO]-"}[shorthand]
    return [("[M-H]-", "MH", "ok", False, False),
            (sh, shorthand, "ok", False, False),
            ("[M-H-H2O]-", "MH-H2O", "isf", True, False)]

# --- keep with a raw-formula adduct alongside [M-H]-: tests Rule C crux, ----
#     but gate keeps trivially (decouples under-flag crux from the gate). -----
def keep_rawformula_withMH(raw="acetate"):
    raws = {"acetate": "[M+C2H3O2]-", "formate": "[M+COOH]-", "TFA": "[M+C2F3O2]-"}
    return [("[M-H]-", "MH", "ok", False, False),
            (raws[raw], raw, "dubious", False, True),
            ("[M-H-CO2]-", "MH-CO2", "isf", True, False)]

# --- DURABLE over-flag: the ONLY intact ion is an obscure-but-validated ----
#     adduct. gold=keep; a no-skill model over-flags it -> wrongly flags. -----
def keep_validated_only(ions):
    # ions: list of (annotation, ion-key) that are validated-OK (crux)
    feats = [(ann, key, "ok", False, True) for ann, key in ions]
    feats.append(("[M-H-H2O]-", "MH-H2O", "isf", True, False))
    return feats

# --- under-flag (eroding but real): only raw-formula dubious + isf, no MH. --
#     gold=flag; a no-skill model that reads raw-formula as ok wrongly keeps. -
def flag_rawformula(raw_a="acetate", raw_b="formate"):
    raws = {"acetate": "[M+C2H3O2]-", "formate": "[M+COOH]-", "TFA": "[M+C2F3O2]-"}
    return [(raws[raw_a], raw_a, "dubious", False, True),
            (raws[raw_b], raw_b, "dubious", False, True),
            ("[M-H-H2O]-", "MH-H2O", "isf", True, False)]

# (compound, M, rt, gate, feats)
FAMILIES = [
    # trivial keep controls
    ("citric acid",     192.027003, 9.21,  "keep", keep_std("acetate")),
    ("phenylalanine",   165.078979, 5.22,  "keep", keep_std("formate")),
    # raw-formula crux, gate-trivial keep (Rule C still represented)
    ("malic acid",      134.021524, 8.47,  "keep", keep_rawformula_withMH("acetate")),
    ("glucose",         180.063388, 10.12, "keep", keep_rawformula_withMH("formate")),
    # DURABLE over-flag keeps: only intact ion is a validated isotope/historical form
    ("lactic acid",      90.031694, 6.51,  "keep", keep_validated_only(
        [("[M-H+1i]-", "mono1i"), ("[2M-H+2i]-", "dimer2i")])),
    ("gluconic acid",   196.058303, 11.20, "keep", keep_validated_only(
        [("[2M-H+4i]-", "dimer4i")])),
    ("tartaric acid",   150.016438, 9.03,  "keep", keep_validated_only(
        [("[M-H+1i]-", "mono1i")])),
    ("glutaric acid",   132.042259, 8.01,  "keep", keep_validated_only(
        [("M-H1", "MH")])),
    ("aconitic acid",   174.016438, 9.55,  "keep", keep_validated_only(
        [("[2M-H+2i]-", "dimer2i")])),
    # under-flag flags: only raw-formula dubious + isf
    ("succinic acid",   118.026609, 7.88,  "flag", flag_rawformula("acetate", "formate")),
    ("fumaric acid",    116.010959, 7.05,  "flag", flag_rawformula("formate", "TFA")),
    ("2-oxoglutarate",  146.021524, 8.90,  "flag", flag_rawformula("acetate", "formate")),
]
INTERFERENTS = [("noise1", 250.1812, 4.10), ("noise2", 311.0827, 13.40),
                ("noise3", 288.9934, 3.05), ("noise4", 402.1521, 2.20)]

rows, gold = [], []
gate_principal, crux_map = {}, {}
fid = 1
for comp, M, rt, gate, feats in FAMILIES:
    principal_fid = f"F{fid:02d}"
    gate_principal[principal_fid] = gate
    for j, (ann, ion, cat, isf, crux) in enumerate(feats):
        f = f"F{fid:02d}"
        mz = mz_of(M, ion)
        rows.append({"feature_id": f, "mz": round(mz, 4),
                     "rt": round(rt + (j-1)*0.012, 3), "intensity": 120000//(j+1),
                     "adduct_annotation": ann})
        gold.append({"feature_id": f, "compound": comp,
                     "neutral_mass": round(M, 4), "adduct_category": cat,
                     "is_isf": isf, "gate": gate, "is_crux": crux})
        if crux:
            crux_map[f] = cat
        fid += 1
for name, mz, rt in INTERFERENTS:
    f = f"F{fid:02d}"
    rows.append({"feature_id": f, "mz": mz, "rt": rt,
                 "intensity": 60000, "adduct_annotation": "unknown"})
    gold.append({"feature_id": f, "compound": "(interferent)",
                 "neutral_mass": None, "adduct_category": "unassigned",
                 "is_isf": False, "gate": "n/a", "is_crux": False})
    fid += 1

inp = pd.DataFrame(rows).sort_values("rt").reset_index(drop=True)
pd.DataFrame(gold).to_csv("final_gold.csv", index=False)
inp.to_csv("final_input.csv", index=False)

g = pd.DataFrame(gold)
n_feat = len(rows)
print("category counts:\n", g[g.adduct_category != "unassigned"].groupby("adduct_category").size().to_string())
print(f"\ntotal features: {n_feat}  crux: {len(crux_map)}  gates: {len(gate_principal)} "
      f"(keep={sum(v=='keep' for v in gate_principal.values())} "
      f"flag={sum(v=='flag' for v in gate_principal.values())})")

# Emit the literals the verifier / eval need, sorted by feature id.
def _lit(d):
    return "{\n" + "".join(f'    "{k}": "{v}",\n' for k, v in sorted(d.items())) + "}"
print("\n# --- paste into eval.py and tests/test_outputs.py ---")
print(f"N_FEATURES = {n_feat}")
print("GOLD_CRUX =", _lit(crux_map))
print("GOLD_GATE =", _lit(gate_principal))
