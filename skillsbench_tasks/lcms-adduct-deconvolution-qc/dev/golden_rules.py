#!/usr/bin/env python3
"""Deterministic engine for the Seven Golden Rules (Kind & Fiehn 2007, BMC Bioinf 8:105).

This is the GROUND-TRUTH engine used to compute gold labels for the task.
Cutoffs follow the widely-cited values from the paper; the SKILL.md teaches
the same cutoffs. Gold = does a candidate formula PASS all applicable rules.

NOTE: the numeric cutoffs below should be reconciled against the paper's
Table 1 / Table 6 before final submission. They are the commonly-cited values.
"""
import re

# Monoisotopic masses (Da)
MONO = {
    "C": 12.0, "H": 1.00782503, "N": 14.00307401, "O": 15.99491462,
    "P": 30.97376151, "S": 31.97207069, "Si": 27.97692653,
    "F": 18.99840322, "Cl": 34.96885271, "Br": 78.91833760, "Na": 22.98976928,
}
# Common valences used for RDBE / SENIOR
VALENCE = {"C": 4, "H": 1, "N": 3, "O": 2, "P": 3, "S": 2, "Si": 4,
           "F": 1, "Cl": 1, "Br": 1, "Na": 1}


def parse_formula(f):
    """Parse 'C6H12O6' -> {'C':6,'H':12,'O':6}."""
    counts = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", f):
        if not el:
            continue
        counts[el] = counts.get(el, 0) + (int(n) if n else 1)
    return counts


def mono_mass(counts):
    return sum(MONO[e] * n for e, n in counts.items())


def rdbe(counts):
    """RDBE = 1 + 0.5 * sum(n_i * (v_i - 2)). Halogens/H lower it; C/N/P raise it."""
    return 1.0 + 0.5 * sum(n * (VALENCE[e] - 2) for e, n in counts.items())


def check_rules(counts, derivatization="none"):
    """Return (passes: bool, failed_rules: list[str]). Mass assumed < 500 Da tier."""
    failed = []
    C = counts.get("C", 0)
    H = counts.get("H", 0)
    N = counts.get("N", 0)
    O = counts.get("O", 0)
    P = counts.get("P", 0)
    S = counts.get("S", 0)
    Si = counts.get("Si", 0)

    tms = derivatization.lower() == "tms"

    # Rule 1: element-number maxima (<500 Da tier)
    maxima = {"C": 39, "H": 72, "N": 20, "O": 20, "P": 9, "S": 10,
              "Si": 14, "F": 16, "Cl": 10, "Br": 5}
    for e, n in counts.items():
        if e in maxima and n > maxima[e]:
            failed.append(f"R1:{e}>{maxima[e]}")

    # Rule 2: LEWIS/SENIOR + RDBE. Neutral molecule -> RDBE >= 0 and integer.
    r = rdbe(counts)
    if r < 0:
        failed.append("R2:RDBE<0")
    if abs(r - round(r)) > 1e-6:
        failed.append("R2:RDBE_noninteger")  # odd valence sum -> impossible neutral
    # SENIOR: sum of valences >= 2*(atoms-1)
    natoms = sum(counts.values())
    vsum = sum(VALENCE[e] * n for e, n in counts.items())
    if natoms > 1 and vsum < 2 * (natoms - 1):
        failed.append("R2:SENIOR_disconnected")

    if C == 0:
        # No carbon: ratio rules below are undefined; flag unless trivial.
        failed.append("R4:no_carbon")
        return (len(failed) == 0, failed)

    # Rule 4: H/C ratio. Common 0.2-3.1 (relaxed upper for TMS due to Si-CH3 groups).
    hc = H / C
    hc_hi = 3.7 if tms else 3.1
    if hc < 0.2 or hc > hc_hi:
        failed.append(f"R4:H/C={hc:.2f}")

    # Rule 5: heteroatom-to-carbon ratios.
    if N / C > 1.3:
        failed.append(f"R5:N/C={N/C:.2f}")
    if O / C > 1.2:
        failed.append(f"R5:O/C={O/C:.2f}")
    if P / C > 0.3:
        failed.append(f"R5:P/C={P/C:.2f}")
    if S / C > 0.8:
        failed.append(f"R5:S/C={S/C:.2f}")
    # Rule 7: TMS relaxes Si/C; non-TMS formulas with Si are unusual but not auto-fail.

    # Rule 6: combined heteroatom-probability heuristics (paper Table 6).
    if N > 1 and O > 1 and P > 1 and S > 1:
        if not (N <= 10 and O <= 20 and P <= 4 and S <= 3):
            failed.append("R6:NOPS_combo")
    if N > 3 and O > 3 and P > 3:
        if not (N <= 11 and O <= 22 and P <= 6):
            failed.append("R6:NOP_combo")
    if O > 1 and P > 1 and S > 1:
        if not (O <= 14 and P <= 3 and S <= 3):
            failed.append("R6:OPS_combo")
    if N > 6 and O > 6 and S > 6:
        if not (N <= 19 and O <= 14 and S <= 8):
            failed.append("R6:NOS_combo")

    return (len(failed) == 0, failed)


def theoretical_isotopes(counts):
    """Approx relative (%) M+1 and M+2 abundances vs monoisotopic peak."""
    C = counts.get("C", 0); N = counts.get("N", 0); O = counts.get("O", 0)
    S = counts.get("S", 0); Si = counts.get("Si", 0)
    Cl = counts.get("Cl", 0); Br = counts.get("Br", 0); H = counts.get("H", 0)
    # M+1: 13C 1.07%, 15N 0.37%, 2H 0.0115%, 29Si 4.68%, 33S 0.75%
    mp1 = 1.07 * C + 0.37 * N + 0.0115 * H + 4.68 * Si + 0.75 * S
    # M+2: 18O 0.205%, 34S 4.25%, 30Si 3.09%, 37Cl 31.96% each, 81Br 97.3% each,
    #       plus 13C2 combinatorial ~ (1.07*C)^2/200
    mp2 = (0.205 * O + 4.25 * S + 3.09 * Si + 31.96 * Cl + 97.3 * Br
           + (1.07 * C) ** 2 / 200.0)
    return mp1, mp2


if __name__ == "__main__":
    for f in ["C6H12O6", "C8H10N4O2", "CHCl3", "C2H7", "C50H2", "C9H11NO2"]:
        c = parse_formula(f)
        ok, fr = check_rules(c)
        i1, i2 = theoretical_isotopes(c)
        print(f"{f:12s} m={mono_mass(c):8.4f} RDBE={rdbe(c):4.1f} "
              f"pass={ok} M+1={i1:5.1f}% M+2={i2:5.1f}% {fr}")
