"""
build_features_v2.py — Hit-level feature table for Bayesian confidence scoring.

Loads curated Orbitrap HILIC neg+pos CSVs, matches annotations to MassWiki hits
by IK14 (with name fallback for null-SMILES hits), computes per-candidate features,
outputs data/feature_table_v2.csv (~65K rows, one per spectrum×candidate pair).

The confidence model scores the top-1 candidate per spectrum (ranked by
entropy_similarity). This script builds the full candidate table; the scoring
engine (bayesian_score_v2.py) selects the top-1 and scores it.

Key features per candidate: entropy_similarity, delta_mda, signed_delta_rt,
sim_gap, adduct evidence (ISF, ok_adduct, n_adducts), spectral_entropy, polarity,
peak-level MS² metrics (forward_cosine, reverse_cosine, cov_count, cov_int).

Usage:
    python code/build_features_v2.py
"""

import os, sys, json, re
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
RDLogger.DisableLog('rdApp.*')

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

NEG_CSV  = os.path.join(ROOT, 'data', 'Orbitrap_HILIC_negESI_curated_042126.csv')
POS_CSV  = os.path.join(ROOT, 'data', 'Orbitrap_HILIC_posESI_curated_042126.csv')
HITS_V2  = os.path.join(ROOT, 'data', 'orbitrap_hits_v2.csv')
ADDUCT_TAX = os.path.join(ROOT, 'data', 'adduct_taxonomy_oliver.csv')
INCHIKEY_CACHE = os.path.join(ROOT, 'data', 'inchikey_cache.json')

# Libraries trusted enough to keep their empty-SMILES deposits. Entries in these
# libraries with no parseable SMILES represent real compounds whose metadata is
# incomplete — not in silico predictions. For scoring, they only contribute
# entropy_similarity and sim_gap; signed_delta_rt is structurally NaN and
# passes through as 0 logLR. Curator annotations in these rows can still be
# recovered via the name-fallback in pick_scored_rows. Other libraries' empty-
# IK14 rows are dropped entirely (typically in silico / stripped metadata).
TRUSTED_EMPTY_IK14_DBS = {'NIST23'}
QUERY_PEAKS_CACHE   = os.path.join(ROOT, 'data', 'query_peaks_cache_v2.json')
LIBRARY_PEAKS_CACHE = os.path.join(ROOT, 'data', 'library_peaks_cache.json')
OUT_PATH = os.path.join(ROOT, 'data', 'feature_table_v2.csv')

MS2_MATCH_TOLERANCE = 0.01   # Da — m/z tolerance for peak matching, aligned with null-dist work

# ── Monoisotopic atomic masses ──────────────────────────────────────────────

ATOM_MASS = {
    'H':  1.00782503, 'C':  12.00000000, 'N':  14.00307401,
    'O':  15.99491462, 'P':  30.97376151, 'S':  31.97207069,
    'Cl': 34.96885268, 'Br': 78.91833710, 'F':  18.99840322,
    'Na': 22.98922190, 'K':  38.96370668, 'Li':  7.01600455,
    'I': 126.90447300, 'Si': 27.97692653,
}
ELECTRON_MASS = 0.00054858


def formula_mass(formula: str) -> float:
    """Compute monoisotopic mass from a molecular formula string like 'C2H4O2'."""
    tokens = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
    # Check that the regex consumed the entire formula
    reconstructed = ''.join(e + c for e, c in tokens)
    if reconstructed != formula:
        return np.nan
    mass = 0.0
    parsed_any = False
    for elem, count in tokens:
        if not elem:
            continue
        n = int(count) if count else 1
        if elem not in ATOM_MASS:
            return np.nan
        mass += ATOM_MASS[elem] * n
        parsed_any = True
    return mass if parsed_any else np.nan


# ── Adduct parsing ──────────────────────────────────────────────────────────

# Common adduct aliases that don't parse cleanly with regex
ADDUCT_ALIASES = {
    'FA':      'CH2O2',   # formic acid (HCOOH = CH2O2)
    'formate': 'CH2O2',
    'HCOO':    'CHO2',    # formate anion
    'CHO2':    'CHO2',
    'COOH':    'CHO2',
    'acetate': 'C2H4O2',
    'HAc':     'C2H4O2',
    'CH3COO':  'C2H3O2',
    'CH3COOH': 'C2H4O2',
    'TFA':     'C2HF3O2',
    'ACN':     'C2H3N',
    'DMSO':    'C2H6OS',
    'MeOH':    'CH4O',
    'IsoProp': 'C3H8O',
    'NH4':     'NH4',
    'NH3':     'NH3',
    'PO3H2':   'H2O3P',
    'H3PO4':   'H3O4P',
    'C2F3O2':  'C2F3O2',
    'C2H3O2':  'C2H3O2',
    'CHC3COO': 'C4H3ClO2',  # chloroacetate (approximate)
    'e':       '__electron__',
}


def _expand_alias(token: str) -> str:
    """Replace known aliases with molecular formulas."""
    return ADDUCT_ALIASES.get(token, token)


def parse_adduct(adduct_str: str):
    """
    Parse an adduct string into (multiplier, mass_shift, charge).

    Returns (n_molecules, mass_shift, abs_charge) or None if unparseable.
    mass_shift is the total mass added/removed per molecule (before dividing by charge).
    The theoretical m/z = (n * M + mass_shift) / abs_charge

    Examples:
        [M+H]+       → (1, +1.00728, 1)
        [M-H]-       → (1, -1.00728, 1)
        [M+Na]+      → (1, +22.9886, 1)
        [2M-H]-      → (2, -1.00728, 1)
        [M+2H]2+     → (1, +2*1.00728, 2)
        [M+H-H2O]+   → (1, +H - H2O, 1)
        [M+FA-H]-    → (1, +FA - H, 1)
    """
    if not isinstance(adduct_str, str):
        return None
    s = adduct_str.strip()

    # Handle special cases
    if s in ('Unknown', 'not given', 'Precursor ion scan', ''):
        return None

    # Strip enclosing brackets and extract charge
    # Patterns: [...]+ [...]+  [...]- [...]2+ [...]2- or bare M+H etc.
    charge_sign = None
    abs_charge = 1

    # Try to match bracketed form: [content]charge
    m = re.match(r'^\[(.+?)\](\d*)([\+\-]?)(\d*)$', s)
    if m:
        content = m.group(1)
        pre_sign_num = m.group(2)
        sign = m.group(3)
        post_sign_num = m.group(4)
        if sign:
            charge_sign = 1 if sign == '+' else -1
            abs_charge = int(pre_sign_num) if pre_sign_num else (int(post_sign_num) if post_sign_num else 1)
        else:
            # No explicit sign — try to infer
            charge_sign = None
    else:
        # Bare form: M+H, M-H, 2M-H, M+Na, Cat, M-e, etc.
        content = s
        # Try to find trailing +/-
        m2 = re.match(r'^(.+?)([\+\-])$', s)
        if m2:
            content = m2.group(1)
            charge_sign = 1 if m2.group(2) == '+' else -1

    # [Cat]+ / [Anion]- under the unified neutral-reference mass (see get_neutral_mass):
    # cation-form SMILES have already had 1 H subtracted to reach the reference, so [Cat]+
    # behaves identically to [M+H]+; [Anion]- identical to [M-H]-.
    if content == 'Cat':
        return (1, +ATOM_MASS['H'], 1)
    if content == 'Anion':
        return (1, -ATOM_MASS['H'], 1)
    # [Cat-X]-, [Cat+Y]+ etc. — composite cation-relative adducts. Leave unparseable for now.
    if content.startswith('Cat') or content.startswith('Anion'):
        return None

    # Parse multiplier: nM...
    m_mult = re.match(r'^(\d+)?M(.*)$', content)
    if not m_mult:
        return None

    n_molecules = int(m_mult.group(1)) if m_mult.group(1) else 1
    remainder = m_mult.group(2)  # everything after M

    # Parse the additive/subtractive terms: +H, -H2O, +Na, -CO2, etc.
    mass_shift = 0.0

    # Tokenize: split on +/- keeping the sign
    # e.g. "+H-H2O+Na" → [(+, H), (-, H2O), (+, Na)]
    terms = re.findall(r'([\+\-])([^\+\-]+)', remainder)

    for sign_str, token in terms:
        sign = 1 if sign_str == '+' else -1

        # Check for leading number: e.g. "2H" → 2 × H, "2H2O" → 2 × H2O
        num_match = re.match(r'^(\d+)([A-Z].*)$', token)
        if num_match:
            coeff = int(num_match.group(1))
            formula_str = num_match.group(2)
        else:
            coeff = 1
            formula_str = token

        # Handle "i" suffix (isotope label) — ignore
        formula_str = re.sub(r'i$', '', formula_str)
        if not formula_str:
            continue

        # Expand aliases
        formula_str = _expand_alias(formula_str)

        if formula_str == '__electron__':
            mass_shift += sign * coeff * ELECTRON_MASS
            continue

        # Compute mass of this formula
        fm = formula_mass(formula_str)
        if np.isnan(fm):
            return None  # can't parse
        mass_shift += sign * coeff * fm

    # Infer charge from added/removed protons for sloppy multi-charge notation like
    # `[M+2H]+` (user omitted the 2 in the charge). Only safe when the adduct has no
    # OTHER charge-contributing species — otherwise the H-count overshoots the net
    # charge. Example: [2M-2H+Na]- has 2H removed but Na added, so net = -1, not -2.
    # Before this guard, the inference clobbered abs_charge 1→2 on that case and
    # produced theoretical_mz off by ~190 Da for HippuricAcid and similar dimers.
    _CHARGE_CARRYING_NON_H = {
        'Na', 'K', 'NH4', 'Li',                                  # cations
        'Cl', 'Br', 'F', 'I', 'OH',                              # anions
        'FA', 'formate', 'HCOO', 'CHO2', 'COOH',                 # formate-family anions
        'acetate', 'CH3COO', 'HAc', 'CH3COOH', 'C2H3O2',         # acetate-family
        'TFA', 'C2HF3O2', 'C2F3O2',                              # trifluoroacetate
    }
    if abs_charge == 1:
        has_other_charge_species = any(
            token in _CHARGE_CARRYING_NON_H for _, token in terms
        )
        if not has_other_charge_species:
            for sign_str, token in terms:
                m_proton = re.match(r'^(\d+)H$', token)
                if m_proton:
                    n_protons = int(m_proton.group(1))
                    if n_protons > 1:
                        abs_charge = n_protons

    # Infer charge sign if not explicit
    if charge_sign is None:
        # Heuristic: if we added H or Na or K, likely positive; if removed H, likely negative
        if '+H' in remainder or '+Na' in remainder or '+K' in remainder or '+NH4' in remainder:
            charge_sign = 1
        elif '-H' in remainder:
            charge_sign = -1
        elif remainder == '' or remainder == '+':
            charge_sign = 1  # bare [M]+ or M
        else:
            charge_sign = 1  # default

    # Account for electrons: positive ion loses electrons, negative gains
    # m/z = (n*M + mass_shift - charge_sign * abs_charge * electron_mass) / abs_charge
    # But for high-res MS, the electron mass is negligible for most purposes.
    # Standard convention: m/z = (n*M + mass_shift) / abs_charge
    # where mass_shift already accounts for proton mass (not electron).

    return (n_molecules, mass_shift, abs_charge)


def compute_theoretical_mz(neutral_mass: float, adduct_str: str) -> float:
    """Compute theoretical m/z from neutral monoisotopic mass + adduct string."""
    parsed = parse_adduct(adduct_str)
    if parsed is None:
        return np.nan
    n_mol, mass_shift, abs_charge = parsed
    mz = (n_mol * neutral_mass + mass_shift) / abs_charge
    if mz <= 0:
        return np.nan
    return mz


# ── InChIKey cache ──────────────────────────────────────────────────────────

_ik_cache = {}

def _load_ik_cache():
    global _ik_cache
    if os.path.exists(INCHIKEY_CACHE):
        with open(INCHIKEY_CACHE) as f:
            _ik_cache = json.load(f)

def get_ik14(smiles: str) -> str:
    if not isinstance(smiles, str) or not smiles.strip():
        return ''
    if smiles in _ik_cache:
        return _ik_cache[smiles].get('ik14', '')
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ''
    try:
        ik = InchiToInchiKey(MolToInchi(mol))
    except Exception:
        return ''
    if ik:
        _ik_cache[smiles] = {'ik14': ik[:14], 'inchikey': ik}
        return ik[:14]
    return ''


def get_neutral_mass(smiles: str) -> float:
    """Get monoisotopic neutral reference mass from SMILES.

    Adjusts for SMILES that carry a formal charge (permanent cations, zwitterions)
    so that downstream adduct-mass math is consistent regardless of whether the
    library stored e.g. `C[N+](C)(C)R` (cation form) or `C[N+](C)(CR[O-])` (zwitterion).
    Without this adjustment, delta_mda was NaN for [Cat]+ adducts and would be
    systematically off by ~1 H for permanent-cation SMILES under any adduct.
    """
    if not isinstance(smiles, str) or not smiles.strip():
        return np.nan
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.nan
    mass = Descriptors.ExactMolWt(mol)
    fc = Chem.GetFormalCharge(mol)
    # Subtract fc×H_mass to normalize to a "neutral reference" mass. For a cation
    # SMILES (fc=+1), we subtract one H; the adduct math then adds it back via e.g.
    # [M+H]+ or [Cat]+ (which we treat as +H equivalent).
    return mass - fc * ATOM_MASS['H']


# ── MS² peak-level match metrics (forward/reverse cosine + library coverage) ─

_query_peaks_cache = None
_lib_peaks_cache = None


def _load_peak_caches():
    """Lazy-load peak caches into module globals (library cache is ~320 MB)."""
    global _query_peaks_cache, _lib_peaks_cache
    if _query_peaks_cache is None:
        if not os.path.exists(QUERY_PEAKS_CACHE):
            raise FileNotFoundError(f"Query peaks cache not found at {QUERY_PEAKS_CACHE}")
        with open(QUERY_PEAKS_CACHE) as f:
            _query_peaks_cache = json.load(f)
    if _lib_peaks_cache is None:
        if not os.path.exists(LIBRARY_PEAKS_CACHE):
            raise FileNotFoundError(f"Library peaks cache not found at {LIBRARY_PEAKS_CACHE}")
        with open(LIBRARY_PEAKS_CACHE) as f:
            _lib_peaks_cache = json.load(f)
    return _query_peaks_cache, _lib_peaks_cache


def compute_ms2_scores(q_pks, l_pks, tol: float = MS2_MATCH_TOLERANCE):
    """Compute 4 MS² match metrics between a query and library peak list.

    Returns (forward_cosine, reverse_cosine, cov_count, cov_int). Any returns
    NaN-quadruple if either peak list is empty / None.

    Metrics:
      * forward_cosine = Σ_matched(q_i × l_j) / (‖q_full‖ × ‖l_full‖)
            Classical spectral cosine. Query noise peaks penalize via ‖q_full‖.
      * reverse_cosine = Σ_matched(q_i × l_j) / (‖q_matched‖ × ‖l_full‖)
            NIST-style reverse match. Query noise dropped from normalizer, so
            reverse_cosine ≥ forward_cosine always. Designed as a complementary
            diagnostic — by itself under-discriminates wrong vs right hits
            because it ignores the evidence carried by unmatched query peaks.
      * cov_count = (# matched library peaks) / (# total library peaks)
      * cov_int   = (Σ matched library intensity) / (Σ total library intensity)

    Peak-matching: for each library peak, pick the single closest query peak
    within `tol` Da. A query peak can be matched by multiple library peaks
    (rare at tol=0.01 Da).

    Inputs expect the same layout as the cache files:
      q_pks, l_pks = list of [mz, intensity] pairs.
    """
    if not q_pks or not l_pks:
        return np.nan, np.nan, np.nan, np.nan

    q_mz = np.asarray([p[0] for p in q_pks], dtype=float)
    q_int = np.asarray([p[1] for p in q_pks], dtype=float)
    l_mz = np.asarray([p[0] for p in l_pks], dtype=float)
    l_int = np.asarray([p[1] for p in l_pks], dtype=float)

    q_full_norm = np.sqrt((q_int ** 2).sum()) + 1e-12
    l_full_norm = np.sqrt((l_int ** 2).sum()) + 1e-12

    # Sort query by mz so we can use binary search for each library peak
    sort_idx = np.argsort(q_mz)
    q_mz_s = q_mz[sort_idx]
    q_int_s = q_int[sort_idx]

    dot_raw = 0.0
    matched_count = 0
    matched_lib_int_raw = 0.0
    total_lib_int_raw = l_int.sum()
    q_matched_mask_sorted = np.zeros(len(q_mz_s), dtype=bool)

    for i in range(len(l_mz)):
        lmz = l_mz[i]
        j = np.searchsorted(q_mz_s, lmz)
        best_j = -1
        best_diff = tol
        if j < len(q_mz_s):
            d = abs(q_mz_s[j] - lmz)
            if d <= best_diff:
                best_j, best_diff = j, d
        if j > 0:
            d = abs(q_mz_s[j - 1] - lmz)
            if d <= best_diff:
                best_j, best_diff = j - 1, d
        if best_j >= 0:
            dot_raw += q_int_s[best_j] * l_int[i]
            matched_count += 1
            matched_lib_int_raw += l_int[i]
            q_matched_mask_sorted[best_j] = True

    forward_cosine = dot_raw / (q_full_norm * l_full_norm)
    q_matched_norm = np.sqrt((q_int_s[q_matched_mask_sorted] ** 2).sum()) + 1e-12
    reverse_cosine = dot_raw / (q_matched_norm * l_full_norm)
    cov_count = matched_count / len(l_pks)
    cov_int = matched_lib_int_raw / max(1e-12, total_lib_int_raw)
    return float(forward_cosine), float(reverse_cosine), float(cov_count), float(cov_int)


def add_ms2_peak_features(df: pd.DataFrame) -> pd.DataFrame:
    """Populate forward_cosine / reverse_cosine / cov_count / cov_int columns on `df`.

    Requires `wiki_id` and `library_wiki_id` columns. Rows whose wiki_id or
    library_wiki_id lacks a peak-cache entry get NaN for all four metrics.
    See `compute_ms2_scores` for metric definitions.
    """
    query_peaks, lib_peaks = _load_peak_caches()
    n = len(df)
    fc = np.full(n, np.nan)
    rc = np.full(n, np.nan)
    cc = np.full(n, np.nan)
    ci = np.full(n, np.nan)
    wids = df['wiki_id'].values
    lids = df['library_wiki_id'].values
    for i in range(n):
        q = query_peaks.get(wids[i])
        lid = lids[i]
        l = lib_peaks.get(lid) if isinstance(lid, str) else None
        if q is None or l is None:
            continue
        fc[i], rc[i], cc[i], ci[i] = compute_ms2_scores(q, l)
    out = df.copy()
    out['forward_cosine'] = fc
    out['reverse_cosine'] = rc
    out['cov_count'] = cc
    out['cov_int'] = ci
    return out


# ── Adduct features ─────────────────────────────────────────────────────────

def _load_adduct_taxonomy():
    tax = pd.read_csv(ADDUCT_TAX)
    lookup = dict(zip(tax['adduct'].str.strip(), tax['category']))
    # Also register normalized (bracket-stripped) form of every CSV entry so
    # classify_adduct (which norm_adducts first) matches both bracketed and bare forms.
    # Earlier, entries like `[M+HAc-H]-` never matched because the lookup key still
    # carried brackets but classify_adduct stripped them off.
    for adduct_raw, cat in list(lookup.items()):
        norm = norm_adduct(adduct_raw)
        if norm and norm not in lookup:
            lookup[norm] = cat
    # Overrides for bare forms + adducts missing from the CSV entirely
    for a in ['M+H', 'M-H', 'M+Na', 'M+NH4', 'M+K', '2M+H', '2M+Na', '2M+K', '2M+NH4']:
        lookup[a] = 'ok'
    # Permanent cation ionization — common for quaternary ammoniums (trigonelline, carnitines).
    lookup['Cat'] = 'ok'
    lookup['Anion'] = 'ok'
    # Water-loss in-source fragments written `M-H2O+H` (different ordering than `M+H-H2O`
    # which the classify_adduct fallback regex already catches).
    for a in ['M-H2O+H', 'M-2H2O+H', 'M-H2O-H', 'M-2H2O-H', 'M-H20-H']:  # last one: data-entry typo
        lookup.setdefault(a, 'isf')
    return lookup


def norm_adduct(s):
    if not isinstance(s, str):
        return ''
    s = s.strip()
    s = re.sub(r'^\[', '', s)
    s = re.sub(r'\][\+\-]?\d*[\+\-]?$', '', s)
    s = re.sub(r'[\+\-]$', '', s)
    return s.strip()


def classify_adduct(adduct_str, lookup):
    norm = norm_adduct(adduct_str)
    cat = lookup.get(norm)
    if cat:
        return cat
    if re.match(r'M\+H-', norm) or re.match(r'M-H-', norm):
        return 'isf'
    return 'unknown'


def apply_adduct_features(df, adduct_lookup):
    """Add adduct-derived binary/integer features."""
    df = df.copy()
    df['adduct_cat'] = df['adduct'].apply(lambda a: classify_adduct(a, adduct_lookup))
    df['is_isf_adduct'] = (df['adduct_cat'] == 'isf').astype(int)
    df['is_dubious_adduct'] = (df['adduct_cat'] == 'dubious').astype(int)

    name_lower = df['name'].fillna('').str.strip().str.lower()
    # For yy_ entries, strip the prefix for grouping
    name_lower = name_lower.str.replace(r'^yy_\s*', '', regex=True)
    df['_name_lower'] = name_lower

    # has_ok_adduct: does this compound have at least one ok adduct anywhere?
    ok_mask = df['adduct_cat'] == 'ok'
    has_ok = df[ok_mask].groupby('_name_lower')['wiki_id'].count().rename('_has_ok')
    df = df.merge(has_ok.reset_index(), on='_name_lower', how='left')
    df['has_ok_adduct'] = df['_has_ok'].fillna(0).clip(upper=1).astype(int)

    # isf_no_mh: ISF adduct with no ok adduct for the same compound
    df['isf_no_mh'] = ((df['is_isf_adduct'] == 1) & (df['has_ok_adduct'] == 0)).astype(int)

    # n_compound_adducts: number of distinct adducts for this compound
    n_add = df.groupby('_name_lower')['adduct'].nunique().rename('n_compound_adducts')
    df = df.merge(n_add.reset_index(), on='_name_lower', how='left')

    df = df.drop(columns=['_has_ok', '_name_lower'])
    return df


# ── Main pipeline ───────────────────────────────────────────────────────────

def load_curated_spectra():
    """Load neg+pos curated CSVs, assign spectrum-level labels, return annotated only."""
    neg = pd.read_csv(NEG_CSV)
    pos = pd.read_csv(POS_CSV)
    neg['polarity'] = 0  # neg
    pos['polarity'] = 1  # pos
    df = pd.concat([neg, pos], ignore_index=True)

    def assign_label(name):
        if pd.isna(name) or str(name).strip() == '':
            return 'blank'
        if str(name).startswith('yy_'):
            return 'FP'
        if str(name).startswith('zz_'):
            return 'zz'
        return 'TP'

    df['label'] = df['name'].apply(assign_label)
    # Include 'blank' rows (Oliver didn't annotate — identity_score typically < 0.7).
    # They get hit_label=0 by construction (spectrum_label != 'TP'). They're intended for
    # scoring / NoTA validation, NOT training — score_confidence_v2.py filters them
    # out of top1_train before fitting the Bayesian channels.
    df = df[df['label'].isin(['TP', 'FP', 'blank'])].copy().reset_index(drop=True)

    # Drop internal standards — not in libraries, can't be scored
    istd_mask = df['name'].str.contains(r'iSTD|ISTD|internal standard', case=False, na=False)
    n_istd = istd_mask.sum()
    if n_istd > 0:
        df = df[~istd_mask].copy().reset_index(drop=True)
        print(f'  Dropped {n_istd} internal standards')

    # Annotation SMILES: prefer corrected, fall back to raw
    if 'annotation-smiles' in df.columns:
        df['anno_smiles'] = df['annotation-smiles'].fillna(df['smiles'])
    else:
        df['anno_smiles'] = df['smiles']

    print('  Computing IK14 for annotations...')
    df['anno_ik14'] = df['anno_smiles'].apply(get_ik14)

    # Spectral entropy (spectrum-level, shared across all candidates)
    df['spectral_entropy'] = pd.to_numeric(df['entropy'], errors='coerce')

    # Normalized annotation name for fallback matching
    anno_name = df['name'].fillna('').str.strip().str.lower()
    anno_name = anno_name.str.replace(r'^yy_\s*', '', regex=True)
    df['anno_name_lower'] = anno_name

    cols_keep = ['wiki_id', 'name', 'adduct', 'anno_smiles', 'anno_ik14',
                 'anno_name_lower',
                 'precursor_mz', 'rt', 'label', 'polarity', 'spectral_entropy']
    return df[cols_keep].copy()


def build_hit_features(spectra, hits_raw, adduct_lookup):
    """
    Build hit-level feature table: one row per (spectrum, candidate compound).

    For each spectrum, deduplicate hits by IK14 (keep best per compound).
    Per candidate row, compute:
      - hit_label: 1 if this candidate's IK14 matches the annotation's IK14
                   (for TP spectra), 0 otherwise. For yy_ spectra, all = 0.
      - entropy_similarity: spectral match (from API)
      - delta_mda: |observed_mz - theoretical_mz(candidate SMILES + candidate adduct)|
      - signed_delta_rt: delta_predicted_rt from API (observed - predicted for this candidate)
      - sim_gap: this candidate's entropy_sim - next best candidate's entropy_sim
      - adduct features: is_isf_adduct, is_dubious_adduct for this candidate's adduct
      - rank: original rank from the library search
      - n_candidates: total candidates for this spectrum
    """
    print('  Computing IK14 for hits...')
    # 'ref_identity' is the original MassWiki API split (identity_search under reference_library);
    # 'reference' is the same thing under a newer API response (where the identity/neutral_loss
    # split is already applied upstream in the fetch script). Both are valid identity-search
    # hits for our scoring pipeline. 'annotation' (in-house lab annotations) and
    # 'ref_neutral_loss' (open search) are NOT used by the scorer.
    hits = hits_raw[hits_raw['hit_source'].isin(['ref_identity', 'reference'])].copy()
    hits['hit_ik14'] = hits['smiles'].fillna('').apply(get_ik14)
    hits['entropy_similarity'] = pd.to_numeric(hits['entropy_similarity'], errors='coerce')
    hits['delta_predicted_rt'] = pd.to_numeric(hits['delta_predicted_rt'], errors='coerce')

    # Compute neutral mass and theoretical m/z per hit
    print('  Computing theoretical m/z for all hits...')
    hits['hit_neutral_mass'] = hits['smiles'].apply(get_neutral_mass)
    hits['hit_theoretical_mz'] = hits.apply(
        lambda r: compute_theoretical_mz(r['hit_neutral_mass'], r['adduct'])
        if pd.notna(r['hit_neutral_mass']) else np.nan, axis=1)

    # Classify hit adducts
    hits['hit_adduct_cat'] = hits['adduct'].apply(lambda a: classify_adduct(a, adduct_lookup))
    hits['hit_is_isf'] = (hits['hit_adduct_cat'] == 'isf').astype(int)
    hits['hit_is_dubious'] = (hits['hit_adduct_cat'] == 'dubious').astype(int)

    # Build spectrum lookup: wiki_id → (anno_ik14, label, polarity, ...)
    spec_info = spectra.set_index('wiki_id')

    # Before dedup: check if each (wiki_id, hit_ik14) has any ok-adduct hit.
    # This answers: "for this candidate compound in this spectrum, does the molecular
    # ion (M+H, M-H, etc.) also appear, or is only an ISF fragment present?"
    print('  Computing per-compound adduct evidence...')
    hits['_is_ok_adduct'] = (hits['hit_adduct_cat'] == 'ok').astype(int)
    ok_adduct_per_compound = (
        hits.groupby(['wiki_id', 'hit_ik14'])['_is_ok_adduct']
        .max().rename('compound_has_ok_adduct').reset_index()
    )
    # Also count distinct adducts per (wiki_id, hit_ik14)
    n_adducts_per_compound = (
        hits.groupby(['wiki_id', 'hit_ik14'])['adduct']
        .nunique().rename('n_candidate_adducts').reset_index()
    )

    # --- Cross-bin, same-RT confirmation (Oliver 2026-06-09) ---------------------------
    # The within-bin flag above only sees adducts present IN THIS bin, so a real in-source
    # fragment scores as an "orphan ISF" even when the compound's clean [M+H] parent is
    # confirmed in another bin at the same RT. Validated (bench_rt_confirmation.py, all 3
    # pre-registered guards pass; non-circular, Δ0 OOF AUC): extend compound_has_ok_adduct
    # to "ok adduct in this bin OR same-compound clean adduct co-eluting in another bin".
    # LABEL-FREE: uses only observed adduct / rt / entropy_similarity, never hit_label.
    RT_CONFIRM_WIN = 10.0       # s — same-compound co-elution window (HILIC peaks are narrow)
    RT_CONFIRM_ESIM_MIN = 0.50  # quality gate on the confirming sibling (label-free)
    rt_by_wid = spec_info['rt']
    conf = hits[(hits['_is_ok_adduct'] == 1) &
                (pd.to_numeric(hits['entropy_similarity'], errors='coerce') >= RT_CONFIRM_ESIM_MIN)]
    conf = conf.assign(_rt=conf['wiki_id'].map(rt_by_wid))
    conf = conf[conf['_rt'].notna()]
    conf_bins = conf.groupby(['hit_ik14', 'wiki_id'])['_rt'].first().reset_index()
    conf_by_ik = {ik: g[['wiki_id', '_rt']].values for ik, g in conf_bins.groupby('hit_ik14')}

    def _rt_confirmed(wid, ik, rt):
        arr = conf_by_ik.get(ik)
        if arr is None or not np.isfinite(rt):
            return 0
        return int(any(w2 != wid and abs(rt - rt2) <= RT_CONFIRM_WIN for w2, rt2 in arr))

    oc = ok_adduct_per_compound
    oc_rt = oc['wiki_id'].map(rt_by_wid).values
    oc['compound_ok_rt_confirmed'] = [
        _rt_confirmed(w, ik, rt) for w, ik, rt in zip(oc['wiki_id'], oc['hit_ik14'], oc_rt)
    ]
    n_rescued = int(((oc['compound_has_ok_adduct'] == 0) & (oc['compound_ok_rt_confirmed'] == 1)).sum())
    oc['compound_has_ok_adduct'] = (
        (oc['compound_has_ok_adduct'] == 1) | (oc['compound_ok_rt_confirmed'] == 1)
    ).astype(int)
    print(f'  RT-confirmation: {n_rescued:,} (wiki_id,compound) pairs gain ok-adduct via co-eluting parent')

    # Dedup per (wiki_id, compound): IK14 first, then name fallback.
    # Policy (2026-04-22): empty-IK14 rows are kept ONLY when from a trusted
    # library (see TRUSTED_EMPTY_IK14_DBS) — those are real compounds with
    # incomplete library metadata. Untrusted-lib empty-IK14 rows (GNPS,
    # MassBank.us, etc.) are dropped as likely in silico / stripped entries.
    print('  Deduplicating hits by IK14 (trusted-lib empty-IK14 kept)...')
    hits['_name_lower'] = hits['name'].fillna('').str.strip().str.lower()
    hits_sorted = hits.sort_values('entropy_similarity', ascending=False)

    # Pass 1: dedup rows with a valid IK14 by (wiki_id, hit_ik14)
    has_ik = hits_sorted['hit_ik14'] != ''
    dedup_ik = hits_sorted[has_ik].drop_duplicates(subset=['wiki_id', 'hit_ik14'], keep='first')

    # Pass 2: empty-IK14 rows — keep only those from trusted libraries, dedup by name
    no_ik = ~has_ik
    null_ik_trusted = hits_sorted[no_ik & hits_sorted['db'].isin(TRUSTED_EMPTY_IK14_DBS)].copy()

    # Drop any trusted empty-IK14 row whose name duplicates an IK14-identified row
    # in the same spectrum (they're just stripped-metadata versions of the same compound).
    dedup_names_per_wid = dedup_ik.groupby('wiki_id')['_name_lower'].apply(set).to_dict()
    keep_mask = [
        r['_name_lower'] == '' or r['_name_lower'] not in dedup_names_per_wid.get(r['wiki_id'], set())
        for _, r in null_ik_trusted.iterrows()
    ]
    null_ik_kept = (null_ik_trusted[keep_mask]
                    .drop_duplicates(subset=['wiki_id', '_name_lower'], keep='first'))

    n_empty_dropped_untrusted = (no_ik & ~hits_sorted['db'].isin(TRUSTED_EMPTY_IK14_DBS)).sum()
    n_empty_dropped_duplicate = len(null_ik_trusted) - len(null_ik_kept)
    dedup = pd.concat([dedup_ik, null_ik_kept], ignore_index=True)
    print(f'  {len(hits):,} raw hits → {len(dedup):,} after dedup')
    print(f'    IK14-populated kept:       {len(dedup_ik):,}')
    print(f'    Trusted empty-IK14 kept:   {len(null_ik_kept):,} (from {sorted(TRUSTED_EMPTY_IK14_DBS)})')
    print(f'    Untrusted empty-IK14 dropped:  {n_empty_dropped_untrusted:,}')
    print(f'    Trusted empty-IK14 dropped (name-duplicate of IK14 row): {n_empty_dropped_duplicate:,}')

    # Merge compound-level adduct evidence onto deduped rows
    dedup = dedup.merge(ok_adduct_per_compound, on=['wiki_id', 'hit_ik14'], how='left')
    dedup = dedup.merge(n_adducts_per_compound, on=['wiki_id', 'hit_ik14'], how='left')

    # isf_no_ok: this candidate's best hit is ISF AND no ok-adduct hit exists for
    # this compound in this spectrum. Strongest ISF suspicion signal.
    dedup['hit_isf_no_ok'] = (
        (dedup['hit_is_isf'] == 1) &
        (dedup['compound_has_ok_adduct'] == 0)
    ).astype(int)

    # Only keep hits for spectra in our labeled set
    labeled_wids = set(spectra['wiki_id'])
    dedup = dedup[dedup['wiki_id'].isin(labeled_wids)].copy()
    print(f'  {len(dedup):,} hits for {dedup["wiki_id"].nunique()} labeled spectra')

    # Compute per-candidate features
    print('  Computing per-candidate features...')

    # observed precursor_mz from spectra
    obs_mz = spec_info['precursor_mz']
    dedup['obs_mz'] = dedup['wiki_id'].map(obs_mz)

    # delta_mda: |observed - theoretical| * 1000
    dedup['delta_mda'] = (dedup['obs_mz'] - dedup['hit_theoretical_mz']).abs() * 1000

    # signed_delta_rt: from API (observed RT - predicted RT for this candidate)
    dedup['signed_delta_rt'] = dedup['delta_predicted_rt']

    # sim_gap per candidate: this candidate's sim - next best candidate's sim (per spectrum).
    # Empty-IK14 rows are excluded from the competitor pool — a candidate we cannot
    # biologically identify should not count as a confidence-defeating competitor.
    # See project_bio_id_dedup.md (agreed interim 2026-04-21).
    print('  Computing sim_gap per candidate...')
    def compute_sim_gaps(group):
        competitor_sims = group.loc[group['hit_ik14'] != '', 'entropy_similarity'].values
        sims = group['entropy_similarity'].values
        group = group.copy()
        if len(competitor_sims) <= 1:
            # No other identifiable compound to compare against.
            # Top-identified candidate gets gap = its own sim (distance from zero baseline);
            # empty-IK14 rows and non-top rows get 0.0.
            if len(competitor_sims) == 1:
                top_sim = competitor_sims[0]
                group['sim_gap'] = [
                    s if (ik != '' and s == top_sim) else 0.0
                    for s, ik in zip(sims, group['hit_ik14'])
                ]
            else:
                group['sim_gap'] = 0.0
            return group
        sorted_comp = np.sort(competitor_sims)[::-1]  # descending
        top_comp = sorted_comp[0]
        runner_up_comp = sorted_comp[1]
        gaps = []
        for s, ik in zip(sims, group['hit_ik14']):
            if ik == '':
                # Phantom — not used as competitor, gap reports distance from top identifiable.
                gaps.append(s - top_comp)
            elif s == top_comp:
                gaps.append(s - runner_up_comp)
            else:
                gaps.append(s - top_comp)
        group['sim_gap'] = gaps
        return group

    dedup = dedup.groupby('wiki_id', group_keys=False).apply(compute_sim_gaps)

    # MS² peak-level match metrics: forward_cosine, reverse_cosine, cov_count, cov_int.
    # Computed from query + library peak caches — see compute_ms2_scores for definitions.
    # All four are correlated with entropy_similarity (r ≈ 0.73–0.84) so they are kept
    # as columns for GBM / diagnostic use, NOT wired into the 3-channel Bayesian.
    print('  Computing MS² peak-level match metrics (forward/reverse cosine, coverage)...')
    dedup = add_ms2_peak_features(dedup)

    # n_candidates per spectrum
    n_cand = dedup.groupby('wiki_id').size().rename('n_candidates')
    dedup = dedup.merge(n_cand.reset_index(), on='wiki_id', how='left')

    # hit_label: 1 if this candidate is the correct compound
    # For TP spectra: match by IK14 first, fall back to name if IK14 unavailable
    # For FP (yy_) spectra: all candidates = 0 (nothing is correct)
    anno_ik14s = spec_info['anno_ik14']
    anno_names = spec_info['anno_name_lower']
    labels = spec_info['label']
    dedup['spectrum_label'] = dedup['wiki_id'].map(labels)
    dedup['anno_ik14'] = dedup['wiki_id'].map(anno_ik14s)
    dedup['anno_name_lower'] = dedup['wiki_id'].map(anno_names)
    dedup['hit_name_lower'] = dedup['name'].fillna('').str.strip().str.lower()

    # Primary: IK14 match (when both sides have IK14)
    ik14_match = (
        (dedup['hit_ik14'] != '') &
        (dedup['anno_ik14'] != '') &
        (dedup['hit_ik14'] == dedup['anno_ik14'])
    )
    # Fallback: name match (when hit has no IK14, e.g. SMILES is null)
    name_match = (
        (~ik14_match) &  # didn't match by IK14
        (dedup['hit_ik14'] == '') &  # hit has no IK14 (SMILES missing)
        (dedup['anno_name_lower'] != '') &
        (dedup['hit_name_lower'] != '') &
        (dedup['anno_name_lower'] == dedup['hit_name_lower'])
    )

    dedup['hit_label'] = (
        (dedup['spectrum_label'] == 'TP') &
        (ik14_match | name_match)
    ).astype(int)

    # Spectrum-level context features
    dedup['polarity'] = dedup['wiki_id'].map(spec_info['polarity'])
    dedup['spectral_entropy'] = dedup['wiki_id'].map(spec_info['spectral_entropy'])
    # Observed retention time (per bin) — needed for cross-bin same-RT confirmation
    # (Oliver 2026-06-09: wipe the ISF penalty when the clean-adduct parent co-elutes).
    dedup['rt_obs'] = dedup['wiki_id'].map(spec_info['rt'])

    # Select output columns
    out_cols = [
        # Keys
        'wiki_id', 'hit_ik14', 'library_wiki_id',
        # Label
        'hit_label', 'spectrum_label',
        # Per-candidate features (channels)
        'entropy_similarity', 'delta_mda', 'signed_delta_rt', 'sim_gap',
        'hit_is_isf', 'hit_is_dubious', 'hit_isf_no_ok',
        'compound_has_ok_adduct', 'compound_ok_rt_confirmed', 'n_candidate_adducts',
        # MS² peak-level match features (see compute_ms2_scores docstring).
        # Kept for GBM / diagnostic use; not in 3-channel Bayesian scorer.
        'forward_cosine', 'reverse_cosine', 'cov_count', 'cov_int',
        # Per-candidate metadata
        'name', 'adduct', 'hit_adduct_cat', 'rank', 'db',
        'hit_theoretical_mz', 'precursor_mz',
        # Spectrum-level context
        'polarity', 'spectral_entropy', 'n_candidates', 'rt_obs',
        # For validation grouping + curator name-fallback match in pick_scored_rows
        'anno_ik14', 'anno_name_lower',
    ]
    result = dedup[[c for c in out_cols if c in dedup.columns]].copy()
    return result


def main():
    print('Phase 1: Building hit-level feature table v2')
    print('=' * 60)

    # Load caches
    print('Loading caches...')
    _load_ik_cache()
    adduct_lookup = _load_adduct_taxonomy()

    # Load curated spectra
    print('Loading curated spectra...')
    spectra = load_curated_spectra()
    print(f'  {len(spectra)} annotated spectra '
          f'(TP={sum(spectra["label"]=="TP")}, FP={sum(spectra["label"]=="FP")})')

    # Load hits
    print('Loading hits_v2...')
    hits_raw = pd.read_csv(HITS_V2)
    id_count = sum(hits_raw['hit_source'] == 'ref_identity')
    print(f'  {len(hits_raw):,} total hits, {id_count:,} identity-search')

    # Build hit-level feature table
    print('Building hit-level features...')
    result = build_hit_features(spectra, hits_raw, adduct_lookup)

    # Summary
    tp_hits = result['hit_label'].sum()
    fp_hits = len(result) - tp_hits
    print()
    print('Hit-level feature table summary:')
    print(f'  Rows: {len(result):,}')
    print(f'  Spectra: {result["wiki_id"].nunique():,}')
    print(f'  hit_label=1 (correct candidate): {tp_hits:,}')
    print(f'  hit_label=0 (wrong candidate): {fp_hits:,}')
    print(f'  Ratio: 1:{fp_hits/max(tp_hits,1):.0f}')
    print()
    for col in ['entropy_similarity', 'delta_mda', 'signed_delta_rt', 'sim_gap',
                'spectral_entropy',
                'forward_cosine', 'reverse_cosine', 'cov_count', 'cov_int']:
        if col not in result.columns:
            continue
        n = result[col].notna().sum()
        print(f'  {col}: {n:,}/{len(result):,} ({100*n/len(result):.1f}%)')

    # Save
    result.to_csv(OUT_PATH, index=False)
    print(f'\nSaved to {OUT_PATH}')

    # Save updated IK cache
    with open(INCHIKEY_CACHE, 'w') as f:
        json.dump(_ik_cache, f)
    print(f'Updated InChIKey cache ({len(_ik_cache)} entries)')

    return result


if __name__ == '__main__':
    main()
