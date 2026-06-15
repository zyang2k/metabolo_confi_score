"""check_deliverable_gate.py — Pre/post-deliverable quality gate.

Modes
  --pre    Snapshot data/deliverable_scores_v2.csv to a timestamped .bak
           before a scoring script overwrites it. Exit 0 (informational).
  --post   Read the freshly-written deliverable and verify it passes the
           shipping thresholds. Exit 1 on regression, 2 on structural issue,
           0 on pass. Block hooks should treat non-zero as "do not commit /
           do not ship".
  --audit  Same checks as --post but never exits non-zero. Pure report.

Thresholds (sourced from MEMORY — edit here if memory changes):
  OOF AUC          >= 0.913    project_gbm_primary_20260423
  Brier_cal        <= 0.1011   project_calibration_20260422
  ECE_cal          <= 0.010    project_calibration_20260422
  FDR @ conf >=0.9 <= 0.050    project_confidence_score_mvp

Where metrics come from
  Preferred:    data/deliverable_metrics_v2.json (sidecar written by the
                scoring script with OOF AUC / Brier_cal / ECE_cal).
  Fallback:     in-sample metrics computed on the deliverable rows whose
                spectrum_label is TP or FP. These overestimate (the model
                was trained on the same rows) — gate uses tighter in-sample
                thresholds and prints a loud note.

Structural contract checked in every --post / --audit run:
  required columns present, no NaN in confidence_pct, hit_label preserved
  on labeled rows, spectrum_label values within the known vocabulary.

Override
  DELIVERABLE_GATE_FORCE=1   skip all checks, exit 0. Use for intentional
                              regressions you've documented elsewhere.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import time

import numpy as np
import pandas as pd
from sklearn.metrics import brier_score_loss, roc_auc_score

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DELIVERABLE = os.path.join(ROOT, 'data', 'deliverable_scores_v2.csv')
METRICS_SIDECAR = os.path.join(ROOT, 'data', 'deliverable_metrics_v2.json')
BAK_DIR = os.path.join(ROOT, 'data')

REQUIRED_COLS = ['wiki_id', 'spectrum_label', 'hit_label',
                 'annotation', 'confidence', 'confidence_pct']
LABEL_VOCAB = {'tp', 'fp', 'blank', 'uncertain', ''}  # compared case-insensitively

# OOF thresholds (preferred path)
OOF_AUC_MIN = 0.913
BRIER_CAL_MAX = 0.1011
ECE_CAL_MAX = 0.010
FDR_AT_09_MAX = 0.050

# In-sample fallback: the deliverable's confidence is final-fit predictions
# (not OOF), and isotonic is fit on OOF — so in-sample Brier/ECE don't relate
# cleanly to the OOF thresholds. Without the sidecar we only enforce a sanity
# AUC floor; Brier/ECE are reported as WARN-only context, not blockers.
IS_AUC_MIN = 0.94


def ece(p: np.ndarray, y: np.ndarray, nbins: int = 10) -> float:
    edges = np.linspace(0, 1, nbins + 1)
    bi = np.clip(np.digitize(p, edges[1:-1]), 0, nbins - 1)
    total, n = 0.0, 0
    for b in range(nbins):
        m = bi == b
        if m.sum() == 0:
            continue
        total += m.sum() * abs(p[m].mean() - y[m].mean())
        n += m.sum()
    return total / max(n, 1)


def fdr_at(confidence: np.ndarray, hit_label: np.ndarray, threshold: float) -> float:
    mask = confidence >= threshold
    if mask.sum() == 0:
        return 0.0
    tp = (hit_label[mask] == 1).sum()
    fp = (hit_label[mask] == 0).sum()
    return fp / max(tp + fp, 1)


def cmd_pre() -> int:
    if not os.path.exists(DELIVERABLE):
        print(f'[gate --pre] no current {DELIVERABLE}; nothing to snapshot')
        return 0
    ts = time.strftime('%Y%m%d_%H%M%S')
    bak = os.path.join(BAK_DIR, f'deliverable_scores_v2.bak.{ts}.csv')
    shutil.copy2(DELIVERABLE, bak)
    print(f'[gate --pre] snapshot -> {bak}')
    return 0


def check_structure(df: pd.DataFrame, issues: list) -> None:
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        issues.append(('FAIL', f'missing required columns: {missing}'))
        return
    n_bad_conf = df['confidence_pct'].isna().sum()
    if n_bad_conf:
        issues.append(('FAIL', f'{n_bad_conf} rows with NaN confidence_pct'))
    seen = {s.lower() for s in df['spectrum_label'].fillna('').astype(str).unique()}
    bad_labels = seen - LABEL_VOCAB
    if bad_labels:
        issues.append(('WARN', f'unexpected spectrum_label values: {bad_labels}'))
    labeled = df[df['spectrum_label'].astype(str).str.lower().isin(['tp', 'fp'])]
    if labeled.empty:
        issues.append(('FAIL', 'no TP/FP labeled rows in deliverable'))
    if 'hit_label' in df.columns and df['hit_label'].isna().any():
        n = df['hit_label'].isna().sum()
        issues.append(('WARN', f'{n} rows with NaN hit_label'))


def check_metrics_sidecar(metrics: dict, issues: list) -> None:
    auc = metrics.get('oof_auc')
    brier = metrics.get('brier_cal')
    e_cal = metrics.get('ece_cal')
    fdr = metrics.get('fdr_at_09')
    if auc is None or auc < OOF_AUC_MIN:
        issues.append(('FAIL', f'OOF AUC {auc} < {OOF_AUC_MIN}'))
    if brier is None or brier > BRIER_CAL_MAX:
        issues.append(('FAIL', f'Brier_cal {brier} > {BRIER_CAL_MAX}'))
    if e_cal is None or e_cal > ECE_CAL_MAX:
        issues.append(('FAIL', f'ECE_cal {e_cal} > {ECE_CAL_MAX}'))
    if fdr is not None and fdr > FDR_AT_09_MAX:
        issues.append(('FAIL', f'FDR@0.9 {fdr} > {FDR_AT_09_MAX}'))
    print(f'  [sidecar] AUC={auc}  Brier_cal={brier}  ECE_cal={e_cal}  FDR@0.9={fdr}')


def check_metrics_insample(df: pd.DataFrame, issues: list) -> None:
    labeled = df[df['spectrum_label'].astype(str).str.lower().isin(['tp', 'fp'])].copy()
    if labeled.empty:
        return
    y = labeled['hit_label'].astype(int).values
    p = labeled['confidence'].astype(float).values
    auc = roc_auc_score(y, p) if len(set(y)) > 1 else float('nan')
    brier = brier_score_loss(y, p)
    e = ece(p, y)
    fdr = fdr_at(p, y, 0.9)
    print('  [in-sample, no OOF sidecar — AUC is the only blocker; '
          'Brier/ECE shown for context, not enforced]')
    print(f'  AUC={auc:.4f}  Brier={brier:.4f}  ECE={e:.4f}  FDR@0.9={fdr:.4f}')
    if auc < IS_AUC_MIN:
        issues.append(('FAIL', f'in-sample AUC {auc:.4f} < {IS_AUC_MIN}'))
    if fdr > FDR_AT_09_MAX:
        issues.append(('FAIL', f'FDR@0.9 {fdr:.4f} > {FDR_AT_09_MAX}'))


def cmd_check(audit_only: bool) -> int:
    if os.environ.get('DELIVERABLE_GATE_FORCE') == '1':
        print('[gate] DELIVERABLE_GATE_FORCE=1 — skipping all checks')
        return 0
    if not os.path.exists(DELIVERABLE):
        print(f'[gate] no {DELIVERABLE} to check')
        return 0 if audit_only else 2

    df = pd.read_csv(DELIVERABLE, low_memory=False)
    print(f'[gate] {DELIVERABLE}: {len(df):,} rows')

    issues: list[tuple[str, str]] = []
    check_structure(df, issues)
    if os.path.exists(METRICS_SIDECAR):
        with open(METRICS_SIDECAR) as f:
            check_metrics_sidecar(json.load(f), issues)
    else:
        check_metrics_insample(df, issues)
        print(f'  hint: have score_gbm_v2.py dump OOF AUC/Brier/ECE to '
              f'{METRICS_SIDECAR} for stricter checks.')

    fails = [m for lvl, m in issues if lvl == 'FAIL']
    warns = [m for lvl, m in issues if lvl == 'WARN']
    print(f'\n[gate] {len(fails)} FAIL, {len(warns)} WARN')
    for m in fails:
        print(f'  FAIL: {m}')
    for m in warns:
        print(f'  WARN: {m}')

    if audit_only:
        return 0
    if fails:
        print('\n[gate] regression — restore from a .bak in data/ or set '
              'DELIVERABLE_GATE_FORCE=1 to override.')
        return 1
    return 0


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument('--pre', action='store_true', help='snapshot before run')
    p.add_argument('--post', action='store_true', help='check after run, exit non-zero on regress')
    p.add_argument('--audit', action='store_true', help='check without blocking')
    args = p.parse_args()
    if args.pre:
        return cmd_pre()
    if args.post:
        return cmd_check(audit_only=False)
    if args.audit:
        return cmd_check(audit_only=True)
    p.print_help()
    return 0


if __name__ == '__main__':
    sys.exit(main())
