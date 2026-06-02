#!/usr/bin/env python3
"""GC-MS derivatization gold engine + probe cases.

Convention (Fiehn-lab style MeOX/MSTFA GC-MS):
  - TMS (trimethylsilyl) replaces 1 active H on -OH / -COOH(OH) / -SH / N-H.
    Net monoisotopic increment per TMS group: +72.039853 Da  (add Si(CH3)3, remove H)
  - MeOX (methoximation) converts each aldehyde/ketone C=O to C=N-OCH3 FIRST,
    blocking it from TMS. Net increment per carbonyl: +29.026549 Da.
  - A carbonyl that is methoximated does NOT also get a TMS.
  - Protocol "tms": TMS only (carbonyls left underivatized, no increment).
  - Protocol "meox_tms": MeOX on all carbonyls, then TMS on all active-H groups.

We restrict probe compounds to OH/COOH/carbonyl chemistry (no amines/phosphates)
so active-H counting is unambiguous.
"""
TMS = 72.039853
MEOX = 29.026549

# name, underiv_mono_mass, n_OH, n_COOH, n_carbonyl(ald/ket), protocol
CASES = [
    ("lactic acid",     90.031694, 1, 1, 0, "meox_tms"),
    ("glucose",        180.063388, 5, 0, 1, "meox_tms"),
    ("fructose",       180.063388, 5, 0, 1, "meox_tms"),
    ("citric acid",    192.027003, 1, 3, 0, "meox_tms"),
    ("malic acid",     134.021524, 1, 2, 0, "meox_tms"),
    ("succinic acid",  118.026609, 0, 2, 0, "meox_tms"),
    ("fumaric acid",   116.010959, 0, 2, 0, "meox_tms"),
    ("oxalic acid",     89.995309, 0, 2, 0, "meox_tms"),
    ("glycerol",        92.047344, 3, 0, 0, "meox_tms"),
    ("palmitic acid",  256.240230, 0, 1, 0, "meox_tms"),
    ("cholesterol",    386.354866, 1, 0, 0, "meox_tms"),
    ("pyruvic acid",    88.016044, 0, 1, 1, "meox_tms"),
    # protocol contrasts: TMS-only leaves carbonyl bare (tests "carbonyl != TMS")
    ("glucose (TMS only)", 180.063388, 5, 0, 1, "tms"),
    ("pyruvic acid (TMS only)", 88.016044, 0, 1, 1, "tms"),
    ("ribose",         150.052823, 4, 0, 1, "meox_tms"),
]


def gold(n_oh, n_cooh, n_carbonyl, protocol):
    n_tms = n_oh + n_cooh  # each active-H heteroatom site -> 1 TMS
    if protocol == "meox_tms":
        n_meox = n_carbonyl
    else:  # tms only
        n_meox = 0
    return n_tms, n_meox


def deriv_mass(m, n_tms, n_meox):
    return m + n_tms * TMS + n_meox * MEOX


if __name__ == "__main__":
    import pandas as pd
    rows = []
    for name, m, oh, cooh, carb, proto in CASES:
        nt, nm = gold(oh, cooh, carb, proto)
        dm = deriv_mass(m, nt, nm)
        rows.append({"compound": name, "underiv_mass": round(m, 4),
                     "protocol": proto, "n_tms": nt, "n_meox": nm,
                     "deriv_mass": round(dm, 4)})
    df = pd.DataFrame(rows)
    print(df.to_string(index=False))
    df.to_csv("tms_cases_full.csv", index=False)
