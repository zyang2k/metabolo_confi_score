"""
fetch_blanks.py — Fetch MassWiki reference-library hits for blank spectra
(identity_score < 0.7, no curator annotation) in either polarity.

Bypasses masswiki_pipeline_with_token's default filter (which keeps only
manually-annotated, non-yy_ rows). We want the OPPOSITE — unannotated rows
for NoTA validation of the confidence scorer.

Token read from MASSWIKI_TOKEN env var. Never written to disk.

Usage:
    MASSWIKI_TOKEN=... python code/fetch_blanks.py --polarity pos
    MASSWIKI_TOKEN=... python code/fetch_blanks.py --polarity neg
"""
import argparse
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'code'))

import pandas as pd
from masswiki_pipeline_with_token import (
    fetch_reference_hits_many,
    flatten_reference_hits,
)

CUR_TEMPLATE = os.path.join(
    ROOT, 'data', 'Orbitrap_HILIC_{polarity}ESI_curated_042126.csv'
)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--polarity', required=True, choices=['pos', 'neg'])
    args = p.parse_args()

    token = os.environ.get('MASSWIKI_TOKEN')
    if not token:
        sys.exit('Set MASSWIKI_TOKEN env var')

    cur_path = CUR_TEMPLATE.format(polarity=args.polarity)
    out_path = os.path.join(ROOT, 'data', f'{args.polarity}_blanks_hits.csv')
    err_path = os.path.join(ROOT, 'data', f'{args.polarity}_blanks_hits_errors.csv')

    df = pd.read_csv(cur_path, low_memory=False)
    blanks = df[
        (df['identity_score'] < 0.7) &
        (df['name'].isna() | (df['name'].astype(str).str.strip() == ''))
    ]
    wids = blanks['wiki_id'].astype(str).tolist()
    print(f'[{args.polarity}] blanks to fetch: {len(wids):,}')

    hits_dict, errors_df = fetch_reference_hits_many(
        wids, token,
        primary=('binbase', False), secondary=('zyang2k', True),
        max_workers=8, rps=6.0,
    )
    out_df = flatten_reference_hits(hits_dict)
    out_df.to_csv(out_path, index=False)
    print(f'\nSaved {len(out_df):,} hit rows → {out_path}')

    if len(errors_df):
        errors_df.to_csv(err_path, index=False)
        print(f'Saved {len(errors_df):,} error rows → {err_path}')

    if len(out_df):
        counts = out_df.groupby('wiki_id').size()
        print(f'\nSpectra with ≥1 hit: {len(counts):,}')
        print(f'Mean hits per spectrum: {counts.mean():.1f}   median: {int(counts.median())}')


if __name__ == '__main__':
    main()
