"""5-fold GroupKFold OOF training of EvidenceSetTransformer.

Mirrors the GBM's evaluation protocol (score_gbm_v2.py):
  - GroupKFold(n_splits=5) on `anno_ik14` (empty → unique placeholder)
  - Per-fold: train with early stopping on a 10% inner-validation slice of the
    training fold, then predict on the held-out test fold
  - Concatenate test-fold predictions into OOF; isotonic-calibrate the top-1
    score against hit_label
  - Report OOF AUC (top-1 + all-candidates), Brier, ECE; compare vs GBM if
    deliverable_scores_v2.csv exists.

Outputs:
  data/candidate_scores_set_tf.csv     (wiki_id, rank, st_raw, st_cal, hit_label)
  data/deliverable_scores_set_tf.csv   (one row per spectrum, top-1 confidence)

Usage:
  python -m code.set_transformer.train          # full run
  python -m code.set_transformer.train --smoke  # 1 fold, 2 epochs
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from typing import Dict

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import GroupKFold
from torch.utils.data import DataLoader, Subset

# Allow `python code/set_transformer/train.py` direct execution
if __name__ == '__main__' and __package__ is None:
    HERE = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.dirname(HERE))
    from set_transformer.dataset import build_dataset, collate_passthrough
    from set_transformer.model import (EvidenceSetTransformer, STConfig,
                                        masked_bce_per_hit)
else:
    from .dataset import build_dataset, collate_passthrough
    from .model import EvidenceSetTransformer, STConfig, masked_bce_per_hit


# ── Metric helpers ─────────────────────────────────────────────────────────────

def ece(p: np.ndarray, y: np.ndarray, nbins: int = 10) -> float:
    """Expected calibration error, weighted by bin count (matches GBM script)."""
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


def sigmoid(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(x, -50, 50)))


# ── Training loop ──────────────────────────────────────────────────────────────

def train_one_fold(model: EvidenceSetTransformer,
                   train_loader: DataLoader,
                   val_loader: DataLoader,
                   *,
                   n_epochs: int,
                   pos_weight: float,
                   lambda_set: float,
                   lr: float,
                   wd: float,
                   patience: int,
                   verbose: bool = True) -> Dict:
    device = next(model.parameters()).device
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=wd)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=n_epochs)

    best_auc = -np.inf
    best_state = None
    best_epoch = -1
    bad_epochs = 0

    for epoch in range(n_epochs):
        model.train()
        loss_sum, n_seen = 0.0, 0
        for batch in train_loader:
            X = batch['X'].to(device)
            C = batch['C'].to(device)
            M = batch['M'].to(device)
            y_hit = batch['y_hit'].to(device)
            y_set = batch['y_set'].to(device)

            optimizer.zero_grad()
            out = model(X, C, M)
            loss_hit = masked_bce_per_hit(out.equiv_logits, y_hit, M, pos_weight=pos_weight)
            loss_set = nn.functional.binary_cross_entropy_with_logits(out.set_logit, y_set)
            loss = loss_hit + lambda_set * loss_set
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            loss_sum += loss.item() * X.shape[0]
            n_seen += X.shape[0]
        scheduler.step()
        train_loss = loss_sum / max(n_seen, 1)

        # Validation: top-1-by-esim AUC (matches GBM deliverable shape).
        # Top-1 is rank 0 because dataset is sorted by esim desc.
        model.eval()
        v_logits, v_labels = [], []
        with torch.no_grad():
            for batch in val_loader:
                X = batch['X'].to(device)
                C = batch['C'].to(device)
                M = batch['M'].to(device)
                out = model(X, C, M)
                v_logits.append(out.equiv_logits[:, 0].cpu().numpy())
                v_labels.append(batch['y_hit'][:, 0].numpy())
        v_logits = np.concatenate(v_logits)
        v_labels = np.concatenate(v_labels)
        valid = v_labels >= 0
        if valid.sum() < 2 or len(set(v_labels[valid].tolist())) < 2:
            val_auc = float('nan')
        else:
            val_auc = roc_auc_score(v_labels[valid], v_logits[valid])

        if verbose:
            print(f'    epoch {epoch:3d}  train_loss={train_loss:.4f}  val_top1_AUC={val_auc:.4f}')

        if np.isfinite(val_auc) and val_auc > best_auc:
            best_auc = val_auc
            best_epoch = epoch
            best_state = {k: v.detach().clone().cpu() for k, v in model.state_dict().items()}
            bad_epochs = 0
        else:
            bad_epochs += 1
            if bad_epochs >= patience:
                if verbose:
                    print(f'    early stop at epoch {epoch}, best_AUC={best_auc:.4f} @ epoch {best_epoch}')
                break

    if best_state is not None:
        model.load_state_dict(best_state)
    return {'best_val_auc': float(best_auc), 'best_epoch': int(best_epoch)}


# ── Main ───────────────────────────────────────────────────────────────────────

def main(*,
         seed: int = 42,
         K_max: int = 30,
         n_epochs: int = 50,
         batch_size: int = 32,
         lr: float = 3e-4,
         wd: float = 1e-4,
         pos_weight: float = 12.0,
         lambda_set: float = 0.3,
         patience: int = 8,
         d_model: int = 64,
         n_heads: int = 4,
         n_inducing: int = 8,
         n_isab: int = 2,
         dropout: float = 0.1,
         smoke: bool = False) -> None:

    if smoke:
        n_epochs = 2
        patience = 99

    torch.manual_seed(seed)
    np.random.seed(seed)

    HERE = os.path.dirname(os.path.abspath(__file__))
    ROOT = os.path.dirname(os.path.dirname(HERE))
    FEATURE_TABLE = os.path.join(ROOT, 'data', 'feature_table_v2.csv')
    OUT_CAND = os.path.join(ROOT, 'data', 'candidate_scores_set_tf.csv')
    OUT_DELIV = os.path.join(ROOT, 'data', 'deliverable_scores_set_tf.csv')
    OUT_EMB = os.path.join(ROOT, 'data', 'st_oof_embeddings.npz')
    GBM_DELIV = os.path.join(ROOT, 'data', 'deliverable_scores_v2.csv')

    print(f'Building dataset (K_max={K_max})...')
    t0 = time.time()
    ds = build_dataset(FEATURE_TABLE, K_max=K_max, labeled_only=True)
    n_spec = len(ds)
    print(f'  built in {time.time()-t0:.1f}s')
    print(f'  n_spectra: {n_spec:,}   total candidates: {ds.M.sum().item():,}')
    print(f'  set-level prior P(any TP): {ds.y_set.mean().item():.3f}')
    n_valid_hits = int(ds.M.sum().item())
    n_tp_hits = int((ds.y_hit == 1).sum().item())
    print(f'  hit-level prior P(TP):     {n_tp_hits / max(n_valid_hits,1):.4f}')

    cfg = STConfig(
        d_num=ds.d_num,
        cat_cardinalities=ds.cat_cardinalities,
        cat_emb_dim=8,
        d_model=d_model,
        n_heads=n_heads,
        n_inducing=n_inducing,
        n_isab=n_isab,
        dropout=dropout,
    )

    gkf = GroupKFold(n_splits=5)
    splits = list(gkf.split(np.arange(n_spec), ds.y_set.numpy(), ds.groups))
    if smoke:
        splits = splits[:1]

    oof_logits_per_hit = np.full((n_spec, K_max), np.nan, dtype=np.float32)
    oof_logits_set = np.full(n_spec, np.nan, dtype=np.float32)
    oof_z_set = np.full((n_spec, d_model), np.nan, dtype=np.float32)
    oof_z_top = np.full((n_spec, d_model), np.nan, dtype=np.float32)
    fold_records = []

    rng = np.random.default_rng(seed)

    for fold, (tr_idx, te_idx) in enumerate(splits):
        print(f'\n=== Fold {fold} ===  train={len(tr_idx):,}  test={len(te_idx):,}')

        # 10% inner validation for early stopping, sampled once per fold.
        # Note: we don't re-group by anno_ik14 here (would need GroupShuffleSplit);
        # the OOF test set is already group-correct, so the inner val is just for
        # early-stopping signal and small leakage there is acceptable.
        perm = rng.permutation(len(tr_idx))
        n_val = max(int(0.1 * len(tr_idx)), 50)
        val_inner = tr_idx[perm[:n_val]]
        tr_inner = tr_idx[perm[n_val:]]

        train_loader = DataLoader(Subset(ds, tr_inner.tolist()),
                                  batch_size=batch_size, shuffle=True,
                                  collate_fn=collate_passthrough)
        val_loader = DataLoader(Subset(ds, val_inner.tolist()),
                                batch_size=batch_size, shuffle=False,
                                collate_fn=collate_passthrough)
        test_loader = DataLoader(Subset(ds, te_idx.tolist()),
                                 batch_size=batch_size, shuffle=False,
                                 collate_fn=collate_passthrough)

        model = EvidenceSetTransformer(cfg)
        info = train_one_fold(model, train_loader, val_loader,
                              n_epochs=n_epochs, pos_weight=pos_weight,
                              lambda_set=lambda_set, lr=lr, wd=wd,
                              patience=patience)

        # Predict on test fold
        model.eval()
        cursor = 0
        with torch.no_grad():
            for batch in test_loader:
                X = batch['X']; C = batch['C']; M = batch['M']
                out = model(X, C, M)
                bs = X.shape[0]
                gidx = te_idx[cursor:cursor + bs]
                # Mask padded positions back to NaN for clarity in the OOF dump
                logit_arr = out.equiv_logits.cpu().numpy()
                logit_arr[~M.cpu().numpy()] = np.nan
                oof_logits_per_hit[gidx] = logit_arr
                oof_logits_set[gidx] = out.set_logit.cpu().numpy()
                oof_z_set[gidx] = out.z_set.cpu().numpy()
                oof_z_top[gidx] = out.z_top.cpu().numpy()
                cursor += bs

        # Fold AUC on top-1
        top1_logit = oof_logits_per_hit[te_idx, 0]
        top1_label = ds.y_hit[te_idx, 0].numpy()
        valid = (top1_label >= 0) & np.isfinite(top1_logit)
        if valid.sum() >= 2 and len(set(top1_label[valid])) > 1:
            fold_auc = roc_auc_score(top1_label[valid], top1_logit[valid])
        else:
            fold_auc = float('nan')
        print(f'  best_val_auc={info["best_val_auc"]:.4f} @ epoch {info["best_epoch"]}'
              f'   test_top1_AUC={fold_auc:.4f}')
        fold_records.append({'fold': fold, 'test_top1_auc': fold_auc, **info})

    # ── Aggregate OOF metrics ──
    print(f'\n=== Aggregate OOF ===')
    yh = ds.y_hit.numpy()
    M_np = ds.M.numpy()

    # Top-1 (matches GBM deliverable)
    top1_logit_oof = oof_logits_per_hit[:, 0]
    top1_label_oof = yh[:, 0]
    valid_top1 = (top1_label_oof >= 0) & np.isfinite(top1_logit_oof)
    auc_top1 = roc_auc_score(top1_label_oof[valid_top1], top1_logit_oof[valid_top1])

    # All-candidate
    valid_all = M_np & (yh >= 0) & np.isfinite(oof_logits_per_hit)
    auc_all = roc_auc_score(yh[valid_all], oof_logits_per_hit[valid_all])

    print(f'OOF top-1 AUC    : {auc_top1:.4f}      (vs GBM 0.913)')
    print(f'OOF all-cand AUC : {auc_all:.4f}')
    fold_aucs = [r['test_top1_auc'] for r in fold_records if np.isfinite(r['test_top1_auc'])]
    if fold_aucs:
        print(f'Fold top-1 AUCs  : [{", ".join(f"{a:.4f}" for a in fold_aucs)}]')

    # Calibrate top-1 sigmoid → isotonic on hit_label
    p_top1 = sigmoid(top1_logit_oof[valid_top1])
    iso = IsotonicRegression(out_of_bounds='clip')
    iso.fit(p_top1, top1_label_oof[valid_top1])
    p_top1_cal = iso.transform(p_top1)

    brier_raw = brier_score_loss(top1_label_oof[valid_top1], p_top1)
    brier_cal = brier_score_loss(top1_label_oof[valid_top1], p_top1_cal)
    ece_raw = ece(p_top1, top1_label_oof[valid_top1])
    ece_cal = ece(p_top1_cal, top1_label_oof[valid_top1])
    print(f'Top-1 Brier      : raw={brier_raw:.4f}  cal={brier_cal:.4f}')
    print(f'Top-1 ECE        : raw={ece_raw:.4f}  cal={ece_cal:.4f}')

    # ── Write candidate-level + deliverable outputs ──
    all_p = sigmoid(oof_logits_per_hit)
    flat = all_p.reshape(-1)
    finite = np.isfinite(flat)
    flat_cal = np.full_like(flat, np.nan)
    flat_cal[finite] = iso.transform(flat[finite])
    all_p_cal = flat_cal.reshape(all_p.shape)

    # Per-(spectrum, rank) rows for valid candidates only
    rows = []
    for i, wid in enumerate(ds.wiki_id):
        for k in range(K_max):
            if not M_np[i, k]:
                continue
            rows.append({
                'wiki_id': wid,
                'rank': k,
                'st_raw': float(all_p[i, k]) if np.isfinite(all_p[i, k]) else np.nan,
                'st_cal': float(all_p_cal[i, k]) if np.isfinite(all_p_cal[i, k]) else np.nan,
                'hit_label': int(yh[i, k]),
            })
    cand_df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(OUT_CAND), exist_ok=True)
    cand_df.to_csv(OUT_CAND, index=False)
    print(f'\nWrote {OUT_CAND}: {len(cand_df):,} rows')

    deliv = cand_df[cand_df['rank'] == 0].copy()
    deliv['confidence'] = deliv['st_cal']
    deliv['confidence_pct'] = (deliv['confidence'] * 100).round(1)
    deliv['confidence_raw'] = deliv['st_raw']
    deliv['set_logit'] = oof_logits_set[:len(deliv)]
    deliv.to_csv(OUT_DELIV, index=False)
    print(f'Wrote {OUT_DELIV}: {len(deliv):,} rows')

    # Save OOF embeddings for the hybrid (head C) stack
    np.savez_compressed(OUT_EMB,
                        wiki_id=ds.wiki_id,
                        z_set=oof_z_set,
                        z_top=oof_z_top,
                        set_logit=oof_logits_set,
                        top1_logit=oof_logits_per_hit[:, 0])
    print(f'Wrote {OUT_EMB}: z_set {oof_z_set.shape}, z_top {oof_z_top.shape}')

    # ── Compare to GBM if available ──
    if os.path.exists(GBM_DELIV):
        try:
            gbm = pd.read_csv(GBM_DELIV)[['wiki_id', 'confidence', 'hit_label']].rename(
                columns={'confidence': 'gbm_conf', 'hit_label': 'gbm_hit_label'})
            st_side = deliv[['wiki_id', 'confidence', 'hit_label']].rename(
                columns={'confidence': 'st_conf', 'hit_label': 'st_hit_label'})
            cmp = st_side.merge(gbm, on='wiki_id', how='inner')
            print(f'\n=== ST vs GBM (n={len(cmp):,} shared spectra) ===')
            for thr in [0.5, 0.7, 0.9]:
                g = (cmp['gbm_conf'] >= thr).mean()
                s = (cmp['st_conf'] >= thr).mean()
                print(f'  ≥ {thr:.1f}:  GBM {100*g:5.1f}%   ST {100*s:5.1f}%')
            # Use ST's hit_label (curator-aligned via top-1-by-esim of feature_table_v2);
            # GBM's hit_label refers to the curator pick which may be name_fallback rather
            # than the top-1 row, so the two label vectors aren't always identical.
            valid = cmp['st_hit_label'] >= 0
            if valid.sum() > 1:
                auc_g = roc_auc_score(cmp.loc[valid, 'st_hit_label'], cmp.loc[valid, 'gbm_conf'])
                auc_s = roc_auc_score(cmp.loc[valid, 'st_hit_label'], cmp.loc[valid, 'st_conf'])
                print(f'  AUC on shared (ST top-1 label):  GBM {auc_g:.4f}   ST {auc_s:.4f}')
        except Exception as e:
            print(f'  GBM comparison failed: {e}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--smoke', action='store_true', help='1 fold, 2 epochs')
    parser.add_argument('--epochs', type=int, default=50)
    parser.add_argument('--K_max', type=int, default=30)
    parser.add_argument('--batch_size', type=int, default=32)
    parser.add_argument('--lr', type=float, default=3e-4)
    parser.add_argument('--pos_weight', type=float, default=12.0)
    parser.add_argument('--lambda_set', type=float, default=0.3)
    parser.add_argument('--patience', type=int, default=8)
    parser.add_argument('--d_model', type=int, default=64)
    parser.add_argument('--n_isab', type=int, default=2)
    parser.add_argument('--dropout', type=float, default=0.1)
    args = parser.parse_args()
    main(smoke=args.smoke, n_epochs=args.epochs, K_max=args.K_max,
         batch_size=args.batch_size, lr=args.lr, pos_weight=args.pos_weight,
         lambda_set=args.lambda_set, patience=args.patience,
         d_model=args.d_model, n_isab=args.n_isab, dropout=args.dropout)
