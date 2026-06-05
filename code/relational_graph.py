"""relational_graph.py — frozen-graph relational confusability features (gate 2).

Ships the validated rel_nn_sim / rel_n_confirmed_nbr features (project_relational_kg
+0.0070 AUC, TTOF-passed) into production scoring the leakage-safe way.

A bin's relational features search its observed MS2 against a GRAPH of reference bins,
blocked by polarity + precursor m/z (±10 mDa, lab convention), EXCLUDING graph nodes that
share the query's compound block (anno_ik14, else name) so a bin's own same-compound
replicates can't contribute their labels. Identical neighbourhood semantics to
build_relational_features.py — the difference is the graph is an explicit, freezable object:

  • OOF eval        graph = TRAIN fold only  → a test bin never sees its own fold's labels
                    (GroupKFold(anno_ik14) already keeps same-compound bins together, so the
                    train-fold graph contains zero same-compound nodes — honest by construction).
  • Final scoring   graph = all labeled training bins (frozen). In-sample labeled bins still
                    exclude self + same-compound; out-of-sample production bins are scored
                    against the same frozen graph so the feature distribution matches training.

  rel_nn_sim          max entropy sim to any different-compound graph bin (confusability)
  rel_n_confirmed_nbr # different-compound TP graph bins with sim >= SIM_THRESH
"""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
import ms_entropy as me

MS1_TOL = 0.01       # Da, precursor block half-width (10 mDa)
MS2_TOL = 0.02       # Da, entropy fragment tolerance
SIM_THRESH = 0.70    # entropy sim to count a confirmed spectral edge
REL_COLS = ['rel_nn_sim', 'rel_n_confirmed_nbr']


def _blocks(df: pd.DataFrame) -> np.ndarray:
    """Compound block key: anno_ik14 if present, else name fallback (mirrors builder)."""
    anno = df['anno_ik14'].fillna('').astype(str).values
    name = (df['anno_name_lower'] if 'anno_name_lower' in df.columns
            else pd.Series('', index=df.index)).fillna('').astype(str).values
    return np.where(anno != '', anno, np.char.add('name::', name))


class Graph:
    """Precursor-sorted reference graph, one entry per polarity, for blocked NN search."""

    def __init__(self, df: pd.DataFrame, peaks: dict):
        df = df.drop_duplicates('wiki_id')
        mz = pd.to_numeric(df['precursor_mz'], errors='coerce').values
        keep = ~np.isnan(mz) & df['wiki_id'].isin(peaks).values
        df = df[keep].reset_index(drop=True)
        block = _blocks(df)
        label = pd.to_numeric(df['hit_label'], errors='coerce').fillna(0).values
        self.peaks = peaks
        self.by_pol = {}
        for pol in pd.unique(df['polarity']):
            m = (df['polarity'] == pol).values
            mzp = pd.to_numeric(df.loc[m, 'precursor_mz'], errors='coerce').values
            order = np.argsort(mzp)
            self.by_pol[pol] = {
                'wid': df.loc[m, 'wiki_id'].values[order],
                'mz': mzp[order],
                'block': block[m][order],
                'label': label[m][order],
            }

    def features(self, query_df: pd.DataFrame) -> pd.DataFrame:
        """Per query bin → (wiki_id, rel_nn_sim, rel_n_confirmed_nbr)."""
        q = query_df.drop_duplicates('wiki_id')
        qmz = pd.to_numeric(q['precursor_mz'], errors='coerce').values
        qblock = _blocks(q)
        qwid = q['wiki_id'].values
        qpol = q['polarity'].values
        peaks = self.peaks
        rows = []
        for i in range(len(q)):
            wid = qwid[i]
            g = self.by_pol.get(qpol[i])
            if wid not in peaks or np.isnan(qmz[i]) or g is None:
                rows.append((wid, 0.0, 0))
                continue
            lo = np.searchsorted(g['mz'], qmz[i] - MS1_TOL, 'left')
            hi = np.searchsorted(g['mz'], qmz[i] + MS1_TOL, 'right')
            qp = np.asarray(peaks[wid], dtype=np.float64)
            best, n_conf = 0.0, 0
            for p in range(lo, hi):
                if g['wid'][p] == wid or g['block'][p] == qblock[i]:
                    continue
                s = me.calculate_entropy_similarity(
                    qp, np.asarray(peaks[g['wid'][p]], dtype=np.float64),
                    ms2_tolerance_in_da=MS2_TOL, clean_spectra=True)
                if s > best:
                    best = s
                if s >= SIM_THRESH and g['label'][p] == 1:
                    n_conf += 1
            rows.append((wid, float(best), int(n_conf)))
        return pd.DataFrame(rows, columns=['wiki_id'] + REL_COLS)

    def freeze(self, out_path: str):
        """Serialize the graph (meta + the peaks of its nodes) for out-of-sample scoring."""
        node_wids = set()
        meta = {}
        for pol, g in self.by_pol.items():
            meta[pol] = {k: (v.tolist() if isinstance(v, np.ndarray) else v)
                         for k, v in g.items()}
            node_wids.update(g['wid'].tolist())
        payload = {
            'ms1_tol': MS1_TOL, 'ms2_tol': MS2_TOL, 'sim_thresh': SIM_THRESH,
            'by_pol': meta,
            'peaks': {w: self.peaks[w] for w in node_wids if w in self.peaks},
        }
        with open(out_path, 'w') as f:
            json.dump(payload, f)
        return len(node_wids)


def load_graph(path: str) -> 'Graph':
    """Reconstruct a frozen Graph (for the deploy inference module)."""
    payload = json.load(open(path))
    g = Graph.__new__(Graph)
    g.peaks = {k: np.asarray(v, dtype=np.float64) for k, v in payload['peaks'].items()}
    g.by_pol = {}
    for pol, d in payload['by_pol'].items():
        g.by_pol[pol] = {
            'wid': np.asarray(d['wid'], dtype=object),
            'mz': np.asarray(d['mz'], dtype=float),
            'block': np.asarray(d['block'], dtype=object),
            'label': np.asarray(d['label'], dtype=float),
        }
    return g
