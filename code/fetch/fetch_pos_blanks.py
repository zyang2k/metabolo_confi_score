"""
fetch_pos_blanks.py — Fetch MassWiki reference-library hits for the pos-mode
blank spectra (identity_score < 0.7, no curator annotation).

Bypasses masswiki_pipeline_with_token's default filter (which keeps only
manually-annotated, non-yy_ rows). We want the OPPOSITE — unannotated rows
where we expect library hits to be sparse / weak, to validate the scorer's
NoTA behavior.

Token read from MASSWIKI_TOKEN env var. Never written to disk.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'code'))

import pandas as pd
# --- sibling-import path shim (code/ root) ---
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))

from masswiki_pipeline_with_token import (
    fetch_reference_hits_many,
    flatten_reference_hits,
)

POS_CUR  = os.path.join(ROOT, 'data', 'Orbitrap_HILIC_posESI_curated_042126.csv')
OUT      = os.path.join(ROOT, 'data', 'pos_blanks_hits.csv')
ERR_OUT  = os.path.join(ROOT, 'data', 'pos_blanks_hits_errors.csv')


def main():
    token = os.environ.get('MASSWIKI_TOKEN')
    if not token:
        sys.exit('Set MASSWIKI_TOKEN env var')

    pos = pd.read_csv(POS_CUR, low_memory=False)
    blanks = pos[
        (pos['identity_score'] < 0.7) &
        (pos['name'].isna() | (pos['name'].astype(str).str.strip() == ''))
    ]
    wids = blanks['wiki_id'].astype(str).tolist()
    print(f'Pos blanks to fetch: {len(wids):,}')

    primary = ('binbase', False)
    secondary = ('zyang2k', True)

    hits_dict, errors_df = fetch_reference_hits_many(
        wids, token, primary=primary, secondary=secondary,
        max_workers=8, rps=6.0,
    )

    df = flatten_reference_hits(hits_dict)
    df.to_csv(OUT, index=False)
    print(f'\nSaved {len(df):,} hit rows → {OUT}')

    if len(errors_df):
        errors_df.to_csv(ERR_OUT, index=False)
        print(f'Saved {len(errors_df):,} error rows → {ERR_OUT}')

    # Quick summary
    if len(df):
        print('\nPer-spectrum hit counts:')
        counts = df.groupby('wiki_id').size()
        print(f'  spectra with ≥1 hit: {len(counts):,}')
        print(f'  mean hits per spectrum: {counts.mean():.1f}')
        print(f'  median: {int(counts.median())}')


if __name__ == '__main__':
    main()
