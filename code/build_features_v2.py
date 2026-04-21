"""
build_features_v2.py — Hit-level feature table for Bayesian confidence scoring.

Loads curated Orbitrap HILIC neg+pos CSVs, matches annotations to MassWiki hits
by IK14 (with name fallback for null-SMILES hits), computes per-candidate features,
outputs data/feature_table_v2.csv (~65K rows, one per spectrum×candidate pair).

The confidence model scores the top-1 candidate per spectrum (ranked by
entropy_similarity). This script builds the full candidate table; the scoring
engine (bayesian_score_v2.py) selects the top-1 and scores it.

Key features per candidate: entropy_similarity, delta_mda, signed_delta_rt,
sim_gap, adduct evidence (ISF, ok_adduct, n_adducts), spectral_entropy, polarity.

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

NEG_CSV  = os.path.join(ROOT, 'data', 'Orbitrap_HILIC_negESI_curated_041326.csv')
POS_CSV  = os.path.join(ROOT, 'data', 'Orbitrap_HILIC_posESI_curated_041326.csv')
HITS_V2  = os.path.join(ROOT, 'data', 'orbitrap_hits_v2.csv')
ADDUCT_TAX = os.path.join(ROOT, 'data', 'adduct_taxonomy_oliver.csv')
INCHIKEY_CACHE = os.path.join(ROOT, 'data', 'inchikey_cache.json')
OUT_PATH = os.path.join(ROOT, 'data', 'feature_table_v2.csv')

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

    # Handle Cat (cation) forms — skip, can't compute theoretical mz
    if content.startswith('Cat') or content == 'Anion':
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

    # Infer charge from added protons for multiply-charged species
    # e.g. M+2H or [M+2H] → charge=2, [M-2H] → charge=2
    if abs_charge == 1:
        # Count net protons: look for +nH or -nH terms (not part of larger formula)
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
    """Get monoisotopic neutral mass from SMILES via RDKit."""
    if not isinstance(smiles, str) or not smiles.strip():
        return np.nan
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.nan
    return Descriptors.ExactMolWt(mol)


# ── Adduct features ─────────────────────────────────────────────────────────

def _load_adduct_taxonomy():
    tax = pd.read_csv(ADDUCT_TAX)
    lookup = dict(zip(tax['adduct'].str.strip(), tax['category']))
    # Overrides for bare forms
    for a in ['M+H', 'M-H', 'M+Na', 'M+NH4', 'M+K', '2M+H', '2M+Na', '2M+K', '2M+NH4']:
        lookup[a] = 'ok'
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
    df = df[df['label'].isin(['TP', 'FP'])].copy().reset_index(drop=True)

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
    hits = hits_raw[hits_raw['hit_source'] == 'ref_identity'].copy()
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

    # Deduplicate hits per (wiki_id, compound): IK14 first, then name fallback
    print('  Deduplicating hits by IK14...')
    hits['_name_lower'] = hits['name'].fillna('').str.strip().str.lower()
    hits_sorted = hits.sort_values('entropy_similarity', ascending=False)

    # Pass 1: dedup by IK14 (for hits that have IK14)
    has_ik = hits_sorted['hit_ik14'] != ''
    no_ik = ~has_ik
    dedup_ik = hits_sorted[has_ik].drop_duplicates(subset=['wiki_id', 'hit_ik14'], keep='first')

    # Pass 2: for hits with no IK14, check if their name matches any IK14-deduped
    # hit in the same spectrum. If so, drop them (they're duplicates we can't collapse by IK14).
    dedup_names_per_wid = dedup_ik.groupby('wiki_id')['_name_lower'].apply(set).to_dict()
    null_ik_hits = hits_sorted[no_ik].copy()

    keep_mask = []
    for _, r in null_ik_hits.iterrows():
        existing_names = dedup_names_per_wid.get(r['wiki_id'], set())
        if r['_name_lower'] and r['_name_lower'] in existing_names:
            keep_mask.append(False)  # duplicate of an IK14-identified hit
        else:
            keep_mask.append(True)

    null_ik_kept = null_ik_hits[keep_mask].drop_duplicates(
        subset=['wiki_id', '_name_lower'], keep='first')
    n_dropped = no_ik.sum() - len(null_ik_kept)

    dedup = pd.concat([dedup_ik, null_ik_kept], ignore_index=True)
    print(f'  {len(hits):,} raw hits → {len(dedup):,} after IK14+name dedup '
          f'({n_dropped} null-IK14 duplicates removed)')

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

    # sim_gap per candidate: this candidate's sim - next best candidate's sim (per spectrum)
    print('  Computing sim_gap per candidate...')
    def compute_sim_gaps(group):
        sims = group['entropy_similarity'].values
        if len(sims) <= 1:
            group = group.copy()
            group['sim_gap'] = sims[0] if len(sims) == 1 else 0.0
            return group
        # For each candidate, gap = its sim - max sim among OTHER candidates
        group = group.copy()
        sorted_sims = np.sort(sims)[::-1]  # descending
        gaps = []
        for s in sims:
            if s == sorted_sims[0]:
                # This is the top candidate; gap vs runner-up
                gaps.append(s - sorted_sims[1])
            else:
                # Not top; gap vs top (will be negative or zero)
                gaps.append(s - sorted_sims[0])
        group['sim_gap'] = gaps
        return group

    dedup = dedup.groupby('wiki_id', group_keys=False).apply(compute_sim_gaps)

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

    # Select output columns
    out_cols = [
        # Keys
        'wiki_id', 'hit_ik14', 'library_wiki_id',
        # Label
        'hit_label', 'spectrum_label',
        # Per-candidate features (channels)
        'entropy_similarity', 'delta_mda', 'signed_delta_rt', 'sim_gap',
        'hit_is_isf', 'hit_is_dubious', 'hit_isf_no_ok',
        'compound_has_ok_adduct', 'n_candidate_adducts',
        # Per-candidate metadata
        'name', 'adduct', 'hit_adduct_cat', 'rank', 'db',
        'hit_theoretical_mz', 'precursor_mz',
        # Spectrum-level context
        'polarity', 'spectral_entropy', 'n_candidates',
        # For validation grouping
        'anno_ik14',
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
                'spectral_entropy']:
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
