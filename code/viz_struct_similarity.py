"""viz_struct_similarity.py — Four-panel visualization for the structural-similarity
closure narrative.

Panels:
  A. Within-bin similarity distributions — Morgan-2 Tanimoto vs MolRex cosine
     (the bimodal-vs-well-shaped finding)
  B. TP vs FP distribution of cos_to_gbm_top1 — signal IS there at AUC 0.72
  C. Ablation AUC across 6 arms with per-fold error bars
  D. Recovery: how much each removal costs vs how much struct recovers

Saves figures/struct_similarity_summary.png (multi-panel) plus per-panel PNGs.
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from sklearn.metrics import roc_auc_score, roc_curve

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / 'figures'
OUT_DIR.mkdir(exist_ok=True)

# Match meeting deck style — set per-figure dpi to avoid rcParams interactions
plt.rcParams.update({
    'font.size': 10, 'axes.titlesize': 11, 'axes.labelsize': 10,
    'axes.spines.right': False, 'axes.spines.top': False,
    'legend.frameon': False, 'legend.fontsize': 9,
    'figure.facecolor': 'white',
})

# Colors
C_MORGAN = '#7A4FBE'      # purple — discrete fingerprint
C_MOLREX = '#2E86AB'      # blue   — neural embedding
C_TP     = '#1B998B'      # teal
C_FP     = '#D7263D'      # red
C_ARM    = '#5C677D'      # slate (default bars)
C_DROP   = '#D7263D'      # red (removal/loss)
C_REC    = '#1B998B'      # teal (recovery)


def gather_morgan_tanimoto(sample_size=100, seed=0):
    """Sample bins with ≥5 candidates, return within-bin pairwise Tanimotos."""
    import json
    cand = pd.read_csv(ROOT / 'data' / 'candidate_scores_v2.csv',
                       usecols=['wiki_id', 'hit_ik14', 'library_wiki_id'])
    # Need SMILES per library hit — pull from orbitrap_hits_v2 if present, else feature_table
    hits_path = ROOT / 'data' / 'orbitrap_hits_v2.csv'
    if hits_path.exists():
        h = pd.read_csv(hits_path, low_memory=False, usecols=['library_wiki_id', 'smiles'])
        h = h.drop_duplicates('library_wiki_id').set_index('library_wiki_id')['smiles']
        cand['smiles'] = cand['library_wiki_id'].map(h)
    else:
        raise SystemExit('orbitrap_hits_v2.csv not found — needed for SMILES per candidate')
    cand = cand[cand['smiles'].notna() & (cand['smiles'] != '')].copy()
    rng = np.random.RandomState(seed)
    bin_sizes = cand.groupby('wiki_id').size()
    big_bins = bin_sizes[bin_sizes >= 5].index.tolist()
    sample = rng.choice(big_bins, size=min(sample_size, len(big_bins)), replace=False)
    pairs = []
    for wid in sample:
        smis = cand.loc[cand['wiki_id'] == wid, 'smiles'].tolist()
        fps = []
        for s in smis:
            m = Chem.MolFromSmiles(s)
            if m is None: continue
            fps.append(AllChem.GetMorganFingerprintAsBitVect(m, radius=2, nBits=2048))
        if len(fps) < 2: continue
        for i in range(len(fps)):
            for j in range(i + 1, len(fps)):
                pairs.append(DataStructs.TanimotoSimilarity(fps[i], fps[j]))
    return np.asarray(pairs)


def gather_molrex_cosines(sample_size=100, seed=0):
    e = np.load(ROOT / 'data' / 'molrex_embeddings.npz', allow_pickle=True)
    iks = e['ik14']
    V = e['embedding']
    Vn = V / (np.linalg.norm(V, axis=1, keepdims=True) + 1e-12)
    ik2idx = {ik: i for i, ik in enumerate(iks)}
    cand = pd.read_csv(ROOT / 'data' / 'candidate_scores_v2.csv',
                       usecols=['wiki_id', 'hit_ik14'])
    rng = np.random.RandomState(seed)
    bin_sizes = cand.groupby('wiki_id').size()
    big_bins = bin_sizes[bin_sizes >= 5].index.tolist()
    sample = rng.choice(big_bins, size=min(sample_size, len(big_bins)), replace=False)
    pairs = []
    for wid in sample:
        iks_bin = cand.loc[cand['wiki_id'] == wid, 'hit_ik14'].dropna().tolist()
        iks_bin = [ik for ik in iks_bin if ik in ik2idx]
        if len(iks_bin) < 2: continue
        V_g = Vn[[ik2idx[ik] for ik in iks_bin]]
        C = V_g @ V_g.T
        for i in range(len(iks_bin)):
            for j in range(i + 1, len(iks_bin)):
                pairs.append(float(C[i, j]))
    return np.asarray(pairs)


def gather_tpfp_struct_cosine():
    """Per-(bin, candidate) cos_to_gbm_top1 vs hit_label, full labeled set."""
    e = np.load(ROOT / 'data' / 'molrex_embeddings.npz', allow_pickle=True)
    iks = e['ik14']
    V = e['embedding']
    Vn = V / (np.linalg.norm(V, axis=1, keepdims=True) + 1e-12)
    ik2idx = {ik: i for i, ik in enumerate(iks)}
    cand = pd.read_csv(ROOT / 'data' / 'candidate_scores_v2.csv',
                       usecols=['wiki_id', 'hit_ik14', 'gbm_cal', 'hit_label'])
    cand = cand[cand['hit_ik14'].notna() & (cand['hit_ik14'] != '')].copy()
    cand['mol_idx'] = cand['hit_ik14'].map(ik2idx)
    cand = cand[cand['mol_idx'].notna()].copy()
    cand['mol_idx'] = cand['mol_idx'].astype(int)
    rows = []
    for wid, grp in cand.groupby('wiki_id'):
        if len(grp) < 2: continue
        idx = grp['mol_idx'].values
        V_g = Vn[idx]
        C = V_g @ V_g.T
        np.fill_diagonal(C, np.nan)
        gbm_vals = grp['gbm_cal'].values.copy().astype(float)
        labels = grp['hit_label'].values
        for k in range(len(grp)):
            ge = gbm_vals.copy()
            ge[k] = -np.inf
            top1_local = int(np.argmax(ge))
            rows.append((labels[k], float(C[k, top1_local])))
    df = pd.DataFrame(rows, columns=['hit_label', 'cos_to_gbm_top1'])
    df = df.dropna()
    return df


def panel_A(ax, morgan_t, molrex_c):
    bins = np.linspace(0, 1, 41)
    ax.hist(morgan_t, bins=bins, alpha=0.6, color=C_MORGAN, density=True,
            label=f'Morgan-2 Tanimoto  (median {np.median(morgan_t):.2f})',
            edgecolor='none')
    ax.hist(molrex_c, bins=bins, alpha=0.6, color=C_MOLREX, density=True,
            label=f'MolRex cosine  (median {np.median(molrex_c):.2f})',
            edgecolor='none')
    ax.axvspan(0.3, 0.8, alpha=0.08, color='black', zorder=0)
    ax.text(0.55, ax.get_ylim()[1] * 0.92, 'productive range\nfor competition',
            ha='center', va='top', fontsize=8, color='#444', style='italic',
            transform=ax.get_xaxis_transform())
    morgan_mid = ((morgan_t >= 0.3) & (morgan_t <= 0.8)).mean()
    molrex_mid = ((molrex_c >= 0.3) & (molrex_c <= 0.8)).mean()
    ax.set_xlabel('Within-bin pairwise similarity')
    ax.set_ylabel('Density')
    ax.set_title('A. Within-bin similarity: Morgan-2 vs MolRex\n'
                 f'Morgan: {100*morgan_mid:.0f}% in productive range  ·  '
                 f'MolRex: {100*molrex_mid:.0f}%')
    ax.legend(loc='upper center')
    ax.set_xlim(0, 1.02)


def panel_B(ax, tpfp_df):
    tp = tpfp_df[tpfp_df['hit_label'] == 1]['cos_to_gbm_top1'].values
    fp = tpfp_df[tpfp_df['hit_label'] == 0]['cos_to_gbm_top1'].values
    bins = np.linspace(-0.2, 1.0, 36)
    ax.hist(fp, bins=bins, alpha=0.55, color=C_FP, density=True,
            label=f'FP  (n={len(fp):,}, mean {fp.mean():.2f})', edgecolor='none')
    ax.hist(tp, bins=bins, alpha=0.55, color=C_TP, density=True,
            label=f'TP  (n={len(tp):,}, mean {tp.mean():.2f})', edgecolor='none')
    auc = roc_auc_score(tpfp_df['hit_label'], tpfp_df['cos_to_gbm_top1'])
    ax.set_xlabel('cos(this candidate, bin\'s GBM top-1)  — MolRex space')
    ax.set_ylabel('Density')
    ax.set_title(f'B. TP vs FP — within-bin structural similarity to top-1\n'
                 f'Univariate AUC = {auc:.3f}  ·  Cohen d ≈ +0.78')
    ax.legend(loc='upper left')
    ax.set_xlim(-0.1, 1.02)


def panel_C(ax, summary):
    fold_aucs = {row['tag']: eval(row['fold_aucs']) if isinstance(row['fold_aucs'], str)
                 else row['fold_aucs'] for _, row in summary.iterrows()}
    arms = list(summary['tag'])
    auc = [summary[summary['tag'] == a]['auc_oof'].iloc[0] for a in arms]
    ses = [np.std(fold_aucs[a]) / np.sqrt(len(fold_aucs[a])) for a in arms]

    short_names = {
        'A1 baseline_full': 'A1\nbaseline\n(18 feat)',
        'A2 baseline_full + struct': 'A2\n+struct',
        'A3 minus_entropy_sim': 'A3\n−entropy',
        'A4 minus_entropy_sim + struct': 'A4\n−entropy\n+struct',
        'A5 minus_entropy_sim_AND_sim_gap': 'A5\n−entropy\n−sim_gap',
        'A6 minus_entropy_sim_AND_sim_gap + struct': 'A6\n−both\n+struct',
    }
    labels = [short_names.get(a, a) for a in arms]
    has_struct = ['struct' in a or '+ struct' in a for a in arms]
    colors = [C_MOLREX if hs else C_ARM for hs in has_struct]

    x = np.arange(len(arms))
    ax.bar(x, auc, yerr=ses, color=colors, alpha=0.85, capsize=4, edgecolor='white')
    for i, (a, s) in enumerate(zip(auc, ses)):
        ax.text(i, a + s + 0.0015, f'{a:.4f}', ha='center', va='bottom', fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylim(0.86, 0.92)
    ax.set_ylabel('OOF AUC (5-fold IK14 GroupKFold)')
    ax.set_title('C. Ablation — does struct recover the lost AUC?')
    ax.axhline(auc[0], color='#888', linestyle='--', linewidth=0.6)
    ax.text(0.02, auc[0] + 0.001, 'baseline', fontsize=7, color='#666',
            transform=ax.get_yaxis_transform())


def panel_D(ax, summary):
    fold_aucs = {row['tag']: eval(row['fold_aucs']) if isinstance(row['fold_aucs'], str)
                 else row['fold_aucs'] for _, row in summary.iterrows()}

    def get(tag):
        return np.array(fold_aucs[tag])

    a1 = get('A1 baseline_full')
    a3 = get('A3 minus_entropy_sim')
    a4 = get('A4 minus_entropy_sim + struct')
    a5 = get('A5 minus_entropy_sim_AND_sim_gap')
    a6 = get('A6 minus_entropy_sim_AND_sim_gap + struct')

    loss_e   = (a1 - a3).mean()
    loss_es  = (a1 - a5).mean()
    rec_e    = (a4 - a3).mean()
    rec_es   = (a6 - a5).mean()
    se_loss_e  = (a1 - a3).std() / np.sqrt(len(a1))
    se_loss_es = (a1 - a5).std() / np.sqrt(len(a1))
    se_rec_e   = (a4 - a3).std() / np.sqrt(len(a1))
    se_rec_es  = (a6 - a5).std() / np.sqrt(len(a1))

    cats = ['Remove\nentropy_sim', 'Remove entropy_sim\n+ sim_gap']
    losses    = [loss_e,  loss_es]
    recoveries= [rec_e,   rec_es]
    se_losses = [se_loss_e, se_loss_es]
    se_recs   = [se_rec_e,  se_rec_es]

    x = np.arange(len(cats))
    width = 0.35
    ax.bar(x - width/2, losses, width, yerr=se_losses, color=C_DROP, alpha=0.85,
           capsize=4, edgecolor='white', label='AUC lost by removal')
    ax.bar(x + width/2, recoveries, width, yerr=se_recs, color=C_REC, alpha=0.85,
           capsize=4, edgecolor='white', label='AUC recovered by adding struct')
    for i, v in enumerate(losses):
        ax.text(x[i] - width/2, v + se_losses[i] + 0.0008, f'{v:+.4f}',
                ha='center', va='bottom', fontsize=8, color=C_DROP)
    for i, v in enumerate(recoveries):
        ax.text(x[i] + width/2, v + se_recs[i] + 0.0008, f'{v:+.4f}',
                ha='center', va='bottom', fontsize=8, color=C_REC)

    ax.set_xticks(x)
    ax.set_xticklabels(cats)
    ax.set_ylim(0, max(losses) * 1.15)
    ax.set_ylabel('Δ AUC')
    ax.set_title('D. Struct cannot substitute for MS²-similarity\n'
                 '(red bars = lost; green = recovered with struct)')
    ax.legend(loc='upper left')


def main():
    print('Gathering Morgan-2 Tanimotos...')
    morgan_t = gather_morgan_tanimoto(sample_size=100)
    print(f'  n pairs: {len(morgan_t):,}')

    print('Gathering MolRex cosines...')
    molrex_c = gather_molrex_cosines(sample_size=100)
    print(f'  n pairs: {len(molrex_c):,}')

    print('Gathering TP/FP cos_to_gbm_top1 (full set)...')
    tpfp_df = gather_tpfp_struct_cosine()
    print(f'  n labeled rows: {len(tpfp_df):,}')

    print('Loading ablation summary...')
    summary = pd.read_csv(ROOT / 'data' / 'bench_struct_redundancy_summary.csv')

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=120,
                              constrained_layout=True)
    panel_A(axes[0, 0], morgan_t, molrex_c)
    panel_B(axes[0, 1], tpfp_df)
    panel_C(axes[1, 0], summary)
    panel_D(axes[1, 1], summary)
    fig.suptitle('Within-spectrum structural similarity — empirically closed across three encodings',
                 fontsize=13, fontweight='bold')
    out = OUT_DIR / 'struct_similarity_summary.png'
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f'\nWrote {out}')

    # Per-panel for slide-by-slide use
    for letter, panel_fn, data in [
        ('A', panel_A, (morgan_t, molrex_c)),
        ('B', panel_B, (tpfp_df,)),
        ('C', panel_C, (summary,)),
        ('D', panel_D, (summary,)),
    ]:
        f2, a2 = plt.subplots(figsize=(7, 5), dpi=120, constrained_layout=True)
        panel_fn(a2, *data)
        p = OUT_DIR / f'struct_similarity_{letter}.png'
        f2.savefig(p, dpi=150)
        plt.close(f2)
        print(f'  + {p}')


if __name__ == '__main__':
    main()
