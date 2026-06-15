"""Self-contained SMILES -> Graphormer3D embedding.

Single-file utility: given a Graphormer3D checkpoint and a SMILES string, run
the encoder forward pass and return the pooled graph embedding as a numpy
array. The featurizer (atom features, 3D conformer) and the encoder model
definition are both inlined so this file has no project-internal imports.

External deps: numpy, torch, rdkit.

Checkpoint formats accepted:
  * encoder-only ``encoder_best.pt``: ``{"encoder": state_dict,
    "encoder_cfg": {...}, ...}`` produced by pretrain_graphomer3d_pt.py.
  * full regressor ``best.pt``: ``{"model": state_dict, "config": cfg, ...}``
    produced by src/train.py. ``model.encoder.*`` keys are extracted and
    encoder config is read from ``config["model"]``.

Usage:
    from scripts.embed_graphormer3d import embed_smiles
    emb = embed_smiles("/path/to/checkpoint.pt", "CCO")
    print(emb.shape)  # (d_model,)

CLI:
    python scripts/embed_graphormer3d.py <checkpoint.pt> <SMILES>
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")


# =============================================================================
# Featurizer (OGB-style; mirror of src/data/featurizer.py)
# =============================================================================

ATOM_FEATURES = {
    "atomic_num": list(range(1, 119)) + ["misc"],
    "chirality": [
        "CHI_UNSPECIFIED", "CHI_TETRAHEDRAL_CW", "CHI_TETRAHEDRAL_CCW",
        "CHI_OTHER", "misc",
    ],
    "degree": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, "misc"],
    "formal_charge": [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, "misc"],
    "num_h": [0, 1, 2, 3, 4, 5, 6, 7, 8, "misc"],
    "num_radical": [0, 1, 2, 3, 4, "misc"],
    "hybridization": ["SP", "SP2", "SP3", "SP3D", "SP3D2", "misc"],
    "is_aromatic": [False, True],
    "is_in_ring": [False, True],
}
ATOM_FEATURE_DIMS: List[int] = [len(v) for v in ATOM_FEATURES.values()]


def _safe_index(seq, item) -> int:
    try:
        return seq.index(item)
    except ValueError:
        return len(seq) - 1


def _atom_feature_row(atom: Chem.Atom) -> List[int]:
    return [
        _safe_index(ATOM_FEATURES["atomic_num"], atom.GetAtomicNum()),
        _safe_index(ATOM_FEATURES["chirality"], str(atom.GetChiralTag())),
        _safe_index(ATOM_FEATURES["degree"], atom.GetTotalDegree()),
        _safe_index(ATOM_FEATURES["formal_charge"], atom.GetFormalCharge()),
        _safe_index(ATOM_FEATURES["num_h"], atom.GetTotalNumHs()),
        _safe_index(ATOM_FEATURES["num_radical"], atom.GetNumRadicalElectrons()),
        _safe_index(ATOM_FEATURES["hybridization"], str(atom.GetHybridization())),
        int(atom.GetIsAromatic()),
        int(atom.IsInRing()),
    ]


def _embed_3d_best(
    mol: Chem.Mol,
    n_confs: int = 10,
    mmff: bool = True,
    max_iters: int = 200,
    random_seed: int = 0xC0FFEE,
) -> Optional[np.ndarray]:
    """ETKDGv3 -> K conformers -> MMFF -> pick lowest-energy heavy-atom coords."""
    try:
        n_heavy = mol.GetNumAtoms()
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = int(random_seed)
        params.numThreads = 1
        cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=int(n_confs),
                                               params=params))
        if not cids:
            params.useRandomCoords = True
            cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=int(n_confs),
                                                   params=params))
            if not cids:
                return None
        best_cid = int(cids[0])
        if mmff:
            try:
                results = AllChem.MMFFOptimizeMoleculeConfs(
                    mol, maxIters=int(max_iters), numThreads=1,
                )
            except Exception:
                results = None
            if results is not None:
                best_energy = float("inf")
                for cid, item in zip(cids, results):
                    try:
                        _converged, energy = item
                    except (TypeError, ValueError):
                        continue
                    if energy is None or not np.isfinite(energy):
                        continue
                    if energy < best_energy:
                        best_energy = float(energy)
                        best_cid = int(cid)
        conf = mol.GetConformer(best_cid)
        return np.asarray(
            [list(conf.GetAtomPosition(i)) for i in range(n_heavy)],
            dtype=np.float32,
        )
    except Exception:
        return None


def featurize_smiles(
    smiles: str,
    n_confs: int = 10,
    mmff: bool = True,
    max_iters: int = 200,
) -> Tuple[np.ndarray, np.ndarray]:
    """SMILES -> (atom_features (N, F_atom) int64 raw, coords (N, 3) float32).

    Raises ValueError on parse or embedding failure.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    n_atoms = mol.GetNumAtoms()
    if n_atoms == 0:
        raise ValueError(f"SMILES has zero atoms: {smiles!r}")
    atom_feats = np.asarray(
        [_atom_feature_row(a) for a in mol.GetAtoms()], dtype=np.int64,
    )
    coords = _embed_3d_best(Chem.Mol(mol), n_confs=n_confs, mmff=mmff,
                            max_iters=max_iters)
    if coords is None:
        raise ValueError(
            f"3D conformer embedding failed for SMILES: {smiles!r}. "
            "Graphormer3D requires valid 3D coordinates."
        )
    return atom_feats, coords


# =============================================================================
# Minimal batch container -- only the fields the encoder reads.
# =============================================================================

@dataclass
class _Batch:
    atom_features: torch.Tensor   # (B, N, F_atom) long, +1-shifted (PAD = 0)
    atom_mask: torch.Tensor       # (B, N) bool
    coords: torch.Tensor          # (B, N, 3) float
    has_3d: torch.Tensor          # (B,) bool
    pos_enc: torch.Tensor         # (B, N, pe_in_dim) float (zeros unless fuse)


# =============================================================================
# Graphormer3D model (faithful copy of src/models/graphormer3d.py, only what's
# needed for forward inference).
# =============================================================================

class _MultiEmbedding(nn.Module):
    def __init__(self, dims: List[int], d_model: int, init_std: float = 0.02):
        super().__init__()
        self.embs = nn.ModuleList(
            [nn.Embedding(d + 2, d_model, padding_idx=0) for d in dims]
        )
        for emb in self.embs:
            nn.init.normal_(emb.weight, mean=0.0, std=init_std)
            with torch.no_grad():
                emb.weight[0].zero_()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        out = 0
        for i, emb in enumerate(self.embs):
            out = out + emb(x[..., i])
        return out


class _AtomicNumberEmbedding(nn.Module):
    def __init__(self, d_model: int, init_std: float = 0.02):
        super().__init__()
        self.emb = nn.Embedding(ATOM_FEATURE_DIMS[0] + 2, d_model, padding_idx=0)
        nn.init.normal_(self.emb.weight, mean=0.0, std=init_std)
        with torch.no_grad():
            self.emb.weight[0].zero_()

    def forward(self, atom_features: torch.Tensor) -> torch.Tensor:
        return self.emb(atom_features[..., 0])


class _FuseEncoder(nn.Module):
    def __init__(self, dim_emb: int, dim_pe: int, pe_in_dim: int,
                 pe_model: str = "mlp", pe_layers: int = 2,
                 pe_norm: str = "none"):
        super().__init__()
        self.dim_emb = dim_emb
        self.dim_pe = dim_pe
        self.pe_in_dim = pe_in_dim
        self.atom_dim = dim_emb - dim_pe
        self.pe_norm_kind = pe_norm
        self.atom_emb = (
            _MultiEmbedding(ATOM_FEATURE_DIMS, self.atom_dim)
            if self.atom_dim > 0 else None
        )
        if pe_norm == "batch":
            self.raw_norm = nn.BatchNorm1d(pe_in_dim)
        elif pe_norm == "layer":
            self.raw_norm = nn.LayerNorm(pe_in_dim)
        else:
            self.raw_norm = None

        if pe_model == "linear":
            self.pe_encoder = nn.Linear(pe_in_dim, dim_pe)
        elif pe_model == "neo":
            inner = 2 * dim_pe
            self.pe_encoder = nn.Sequential(
                nn.Linear(pe_in_dim, inner), nn.LayerNorm(inner), nn.GELU(),
                nn.Linear(inner, dim_pe), nn.LayerNorm(dim_pe),
            )
        else:  # mlp
            layers: list[nn.Module] = []
            if pe_layers == 1:
                layers.append(nn.Linear(pe_in_dim, dim_pe))
                layers.append(nn.ReLU())
            else:
                layers.append(nn.Linear(pe_in_dim, 2 * dim_pe))
                layers.append(nn.ReLU())
                for _ in range(pe_layers - 2):
                    layers.append(nn.Linear(2 * dim_pe, 2 * dim_pe))
                    layers.append(nn.ReLU())
                layers.append(nn.Linear(2 * dim_pe, dim_pe))
                layers.append(nn.ReLU())
            self.pe_encoder = nn.Sequential(*layers)

    def forward(self, atom_features: torch.Tensor,
                pos_enc: torch.Tensor,
                atom_mask: torch.Tensor) -> torch.Tensor:
        if self.raw_norm is not None:
            if isinstance(self.raw_norm, nn.BatchNorm1d):
                B, K, D = pos_enc.shape
                flat = pos_enc.reshape(-1, D)
                fm = atom_mask.reshape(-1)
                real = self.raw_norm(flat[fm])
                out = torch.zeros_like(flat)
                out[fm] = real.to(out.dtype)
                pos_enc = out.reshape(B, K, D)
            else:
                pos_enc = self.raw_norm(pos_enc) * atom_mask.unsqueeze(-1).to(pos_enc.dtype)
        pe = self.pe_encoder(pos_enc)
        if self.atom_emb is None:
            return pe
        h = self.atom_emb(atom_features)
        return torch.cat([h, pe], dim=-1)


def _build_node_encoder(kind: str, dim_emb: int, pe_dim: int, pe_in_dim: int,
                        pe_model: str, pe_layers: int, pe_norm: str) -> nn.Module:
    if kind == "atomic":
        return _MultiEmbedding(ATOM_FEATURE_DIMS, dim_emb)
    if kind == "atomic_num_only":
        return _AtomicNumberEmbedding(dim_emb)
    if kind == "fuse":
        return _FuseEncoder(dim_emb=dim_emb, dim_pe=pe_dim, pe_in_dim=pe_in_dim,
                            pe_model=pe_model, pe_layers=pe_layers, pe_norm=pe_norm)
    raise ValueError(f"unknown node_encoder: {kind!r}")


def _apply_node_encoder(node_enc: nn.Module, batch: _Batch) -> torch.Tensor:
    if isinstance(node_enc, _FuseEncoder):
        return node_enc(batch.atom_features, batch.pos_enc, batch.atom_mask)
    return node_enc(batch.atom_features)


class _PairTypeAwareRBFSoftplus(nn.Module):
    def __init__(self, num_kernel: int, num_pair_types: int):
        super().__init__()
        self.num_kernel = int(num_kernel)
        self.scale = nn.Embedding(num_pair_types, 1)
        self.shift = nn.Embedding(num_pair_types, 1)
        self.means = nn.Parameter(torch.empty(num_kernel))
        self.log_stds = nn.Parameter(torch.empty(num_kernel))
        nn.init.ones_(self.scale.weight)
        nn.init.zeros_(self.shift.weight)
        nn.init.uniform_(self.means, 0.0, 3.0)
        nn.init.uniform_(self.log_stds, -0.5, 0.5)

    def forward(self, distance: torch.Tensor, pair_type: torch.Tensor) -> torch.Tensor:
        scale = self.scale(pair_type).squeeze(-1)
        shift = self.shift(pair_type).squeeze(-1)
        x = scale * distance + shift
        std = F.softplus(self.log_stds) + 1e-5
        x = x.unsqueeze(-1)
        return torch.exp(-0.5 * ((x - self.means) / std) ** 2) / (
            std * math.sqrt(2.0 * math.pi)
        )


class _RBFToAttentionBias(nn.Module):
    def __init__(self, num_kernel: int, num_heads: int,
                 hidden_dim: Optional[int], num_layers: int, dropout: float):
        super().__init__()
        hidden = int(hidden_dim or num_kernel)
        layers: list[nn.Module] = [nn.LayerNorm(num_kernel)]
        if num_layers == 1:
            layers.append(nn.Linear(num_kernel, num_heads))
        else:
            layers.extend([nn.Linear(num_kernel, hidden), nn.GELU(),
                           nn.Dropout(dropout)])
            for _ in range(num_layers - 2):
                layers.extend([nn.Linear(hidden, hidden), nn.GELU(),
                               nn.Dropout(dropout)])
            layers.append(nn.Linear(hidden, num_heads))
        self.net = nn.Sequential(*layers)

    def forward(self, rbf: torch.Tensor) -> torch.Tensor:
        return self.net(rbf).permute(0, 3, 1, 2).contiguous()


class _Attention(nn.Module):
    def __init__(self, d_model: int, num_heads: int, attn_dropout: float):
        super().__init__()
        self.num_heads = int(num_heads)
        self.head_dim = d_model // num_heads
        self.scale = self.head_dim ** -0.5
        self.qkv = nn.Linear(d_model, 3 * d_model)
        self.out = nn.Linear(d_model, d_model)
        self.attn_dropout = float(attn_dropout)

    def forward(self, x: torch.Tensor, key_padding_mask: torch.Tensor,
                attn_bias: torch.Tensor) -> torch.Tensor:
        B, S, C = x.shape
        qkv = self.qkv(x).view(B, S, 3, self.num_heads, self.head_dim)
        q, k, v = qkv.unbind(dim=2)
        q = q.transpose(1, 2) * self.scale
        k = k.transpose(1, 2)
        v = v.transpose(1, 2)
        logits = torch.matmul(q, k.transpose(-1, -2))
        logits = logits + attn_bias.to(dtype=logits.dtype)
        logits = logits.masked_fill(
            key_padding_mask.view(B, 1, 1, S),
            torch.finfo(logits.dtype).min,
        )
        probs = F.softmax(logits.float(), dim=-1).to(x.dtype)
        probs = F.dropout(probs, p=self.attn_dropout, training=self.training)
        out = torch.matmul(probs, v).transpose(1, 2).contiguous().view(B, S, C)
        return self.out(out)


class _Layer(nn.Module):
    def __init__(self, d_model: int, num_heads: int, ffn_dim: int,
                 dropout: float, attn_dropout: float,
                 activation_dropout: float, activation: str, layer_norm: str):
        super().__init__()
        self.layer_norm = layer_norm
        self.act = F.gelu if activation == "gelu" else F.relu
        self.ln1 = nn.LayerNorm(d_model)
        self.attn = _Attention(d_model, num_heads, attn_dropout)
        self.drop1 = nn.Dropout(dropout)
        self.ln2 = nn.LayerNorm(d_model)
        self.fc1 = nn.Linear(d_model, ffn_dim)
        self.fc2 = nn.Linear(ffn_dim, d_model)
        self.activation_drop = nn.Dropout(activation_dropout)
        self.drop2 = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor, seq_mask: torch.Tensor,
                attn_bias: torch.Tensor) -> torch.Tensor:
        key_padding_mask = ~seq_mask
        if self.layer_norm == "pre":
            h = self.attn(self.ln1(x), key_padding_mask, attn_bias)
            x = x + self.drop1(h)
            x = x * seq_mask.unsqueeze(-1).to(x.dtype)
            h = self.fc2(self.activation_drop(self.act(self.fc1(self.ln2(x)))))
            x = x + self.drop2(h)
            x = x * seq_mask.unsqueeze(-1).to(x.dtype)
            return x
        h = self.attn(x, key_padding_mask, attn_bias)
        x = self.ln1(x + self.drop1(h))
        x = x * seq_mask.unsqueeze(-1).to(x.dtype)
        h = self.fc2(self.activation_drop(self.act(self.fc1(x))))
        x = self.ln2(x + self.drop2(h))
        x = x * seq_mask.unsqueeze(-1).to(x.dtype)
        return x


@dataclass
class Graphormer3DConfig:
    d_model: int = 256
    num_heads: int = 8
    num_layers: int = 6
    ffn_dim: int = 1024
    dropout: float = 0.1
    attn_dropout: float = 0.1
    activation_dropout: float = 0.0
    activation: str = "gelu"
    layer_norm: str = "post"
    num_kernel: int = 128
    gbf_hidden: Optional[int] = None
    gbf_num_layers: int = 2
    use_graph_token: bool = True
    use_centrality: bool = True
    readout: str = "token"
    require_3d: bool = True
    node_encoder: str = "atomic"
    pe_dim: int = 0
    pe_in_dim: int = 0
    pe_model: str = "mlp"
    pe_layers: int = 2
    pe_norm: str = "none"


class Graphormer3D(nn.Module):
    """Graph-level Graphormer3D encoder. Inference-only reproduction."""

    def __init__(self, cfg: Graphormer3DConfig):
        super().__init__()
        self.cfg = cfg
        self.embed_dim = cfg.d_model

        self.node_encoder = _build_node_encoder(
            cfg.node_encoder, dim_emb=cfg.d_model,
            pe_dim=cfg.pe_dim, pe_in_dim=cfg.pe_in_dim,
            pe_model=cfg.pe_model, pe_layers=cfg.pe_layers, pe_norm=cfg.pe_norm,
        )

        self.num_atom_types = ATOM_FEATURE_DIMS[0] + 2
        self.num_pair_types = self.num_atom_types * self.num_atom_types
        self.rbf = _PairTypeAwareRBFSoftplus(cfg.num_kernel, self.num_pair_types)
        self.rbf_to_bias = _RBFToAttentionBias(
            num_kernel=cfg.num_kernel, num_heads=cfg.num_heads,
            hidden_dim=cfg.gbf_hidden, num_layers=cfg.gbf_num_layers,
            dropout=cfg.activation_dropout,
        )

        if cfg.use_centrality:
            centrality_hidden = int(cfg.gbf_hidden or cfg.num_kernel)
            self.centrality_proj = nn.Sequential(
                nn.LayerNorm(cfg.num_kernel),
                nn.Linear(cfg.num_kernel, centrality_hidden),
                nn.GELU(),
                nn.Dropout(cfg.activation_dropout),
                nn.Linear(centrality_hidden, cfg.d_model),
            )
        else:
            self.centrality_proj = None

        if cfg.use_graph_token:
            self.graph_token = nn.Parameter(torch.empty(1, 1, cfg.d_model))
            self.graph_token_bias = nn.Parameter(torch.zeros(1, cfg.num_heads, 1, 1))
            nn.init.normal_(self.graph_token, mean=0.0, std=0.02)
        else:
            self.graph_token = None
            self.graph_token_bias = None

        self.layers = nn.ModuleList([
            _Layer(cfg.d_model, cfg.num_heads, cfg.ffn_dim, cfg.dropout,
                   cfg.attn_dropout, cfg.activation_dropout, cfg.activation,
                   cfg.layer_norm)
            for _ in range(cfg.num_layers)
        ])
        self.final_ln = (nn.LayerNorm(cfg.d_model)
                         if cfg.layer_norm == "pre" else nn.Identity())

    def _atom_pair_inputs(self, batch: _Batch
                          ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        coords = batch.coords.float()
        dist = torch.cdist(coords, coords, p=2)
        atom_type = batch.atom_features[..., 0]
        pair_type = atom_type.unsqueeze(2) * self.num_atom_types + atom_type.unsqueeze(1)
        atom_pair_mask = batch.atom_mask.unsqueeze(1) & batch.atom_mask.unsqueeze(2)
        atom_pair_mask = atom_pair_mask & batch.has_3d.view(-1, 1, 1)
        return dist, pair_type, atom_pair_mask

    def _make_attention_bias(self, rbf: torch.Tensor,
                              atom_pair_mask: torch.Tensor,
                              atom_mask: torch.Tensor
                              ) -> Tuple[torch.Tensor, torch.Tensor]:
        atom_bias = self.rbf_to_bias(rbf)
        atom_bias = atom_bias * atom_pair_mask.unsqueeze(1).to(atom_bias.dtype)
        if not self.cfg.use_graph_token:
            return atom_bias, atom_mask
        B, H, N, _ = atom_bias.shape
        full_bias = atom_bias.new_zeros(B, H, N + 1, N + 1)
        full_bias[:, :, 1:, 1:] = atom_bias
        full_bias[:, :, 0:1, :] = self.graph_token_bias.to(atom_bias.dtype)
        full_bias[:, :, :, 0:1] = self.graph_token_bias.to(atom_bias.dtype)
        seq_mask = torch.cat([
            torch.ones(B, 1, dtype=torch.bool, device=atom_mask.device),
            atom_mask,
        ], dim=1)
        qk_mask = seq_mask.unsqueeze(1) & seq_mask.unsqueeze(2)
        full_bias = full_bias * qk_mask.unsqueeze(1).to(full_bias.dtype)
        return full_bias, seq_mask

    @staticmethod
    def _masked_readout(x: torch.Tensor, mask: torch.Tensor, mode: str) -> torch.Tensor:
        weights = mask.unsqueeze(-1).to(x.dtype)
        summed = (x * weights).sum(dim=1)
        if mode == "sum":
            return summed
        return summed / weights.sum(dim=1).clamp(min=1.0)

    def forward(self, batch: _Batch) -> torch.Tensor:
        x = _apply_node_encoder(self.node_encoder, batch)
        x = x * batch.atom_mask.unsqueeze(-1).to(x.dtype)

        dist, pair_type, atom_pair_mask = self._atom_pair_inputs(batch)
        rbf = self.rbf(dist, pair_type)
        rbf = rbf * atom_pair_mask.unsqueeze(-1).to(rbf.dtype)

        if self.centrality_proj is not None:
            centrality = self.centrality_proj(rbf.sum(dim=2)).to(x.dtype)
            x = x + centrality
            x = x * batch.atom_mask.unsqueeze(-1).to(x.dtype)

        attn_bias, seq_mask = self._make_attention_bias(
            rbf, atom_pair_mask, batch.atom_mask,
        )

        if self.cfg.use_graph_token:
            token = self.graph_token.expand(x.size(0), -1, -1).to(dtype=x.dtype)
            x = torch.cat([token, x], dim=1)

        for layer in self.layers:
            x = layer(x, seq_mask=seq_mask, attn_bias=attn_bias)
        x = self.final_ln(x)
        x = x * seq_mask.unsqueeze(-1).to(x.dtype)

        if self.cfg.use_graph_token:
            nodes = x[:, 1:]
            if self.cfg.readout == "token":
                graph = x[:, 0]
            else:
                graph = self._masked_readout(nodes, batch.atom_mask, self.cfg.readout)
        else:
            graph = self._masked_readout(x, batch.atom_mask, self.cfg.readout)
        return graph


# =============================================================================
# Checkpoint loading + embedding entrypoint
# =============================================================================

_ENCODER_CFG_FIELDS = {f for f in Graphormer3DConfig.__dataclass_fields__}


def _extract_encoder_state_and_config(
    ckpt: dict,
) -> Tuple[dict, Graphormer3DConfig]:
    """Accept both encoder-only and full-regressor checkpoint formats."""
    # Encoder-only: pretrain_graphomer3d_pt.py writes encoder_best.pt as
    #   {"encoder": <Graphormer3D.state_dict>, "encoder_cfg": {...}}
    if "encoder" in ckpt and "encoder_cfg" in ckpt:
        enc_cfg_dict = dict(ckpt["encoder_cfg"])
        enc_cfg = Graphormer3DConfig(**{
            k: v for k, v in enc_cfg_dict.items() if k in _ENCODER_CFG_FIELDS
        })
        return dict(ckpt["encoder"]), enc_cfg

    # Full regressor: src/train.py writes best.pt as
    #   {"model": <Graphormer3DRegressor.state_dict>, "config": cfg}
    # Strip "encoder." prefix from state-dict keys; build the encoder config
    # from cfg["model"].
    if "model" in ckpt and "config" in ckpt:
        cfg = ckpt["config"]
        model_cfg = dict(cfg.get("model", {}))
        enc_cfg = Graphormer3DConfig(**{
            k: v for k, v in model_cfg.items() if k in _ENCODER_CFG_FIELDS
        })
        full_sd = ckpt["model"]
        encoder_sd = {}
        for k, v in full_sd.items():
            if k.startswith("encoder."):
                encoder_sd[k[len("encoder."):]] = v
        if not encoder_sd:
            raise ValueError(
                "Checkpoint has 'model' but no keys with 'encoder.' prefix; "
                "this loader expects a Graphormer3DRegressor state_dict."
            )
        return encoder_sd, enc_cfg

    raise ValueError(
        "Unrecognized checkpoint format. Expected either "
        "{'encoder', 'encoder_cfg', ...} (pretrain encoder_best.pt) or "
        "{'model', 'config', ...} (train.py best.pt)."
    )


def _make_single_batch(
    atom_feats_raw: np.ndarray,
    coords: np.ndarray,
    pe_in_dim: int,
    device: torch.device,
) -> _Batch:
    """Pack one molecule into a B=1 batch matching collator conventions
    (atom feature ids +1-shifted so PAD=0)."""
    n = atom_feats_raw.shape[0]
    atom_features = torch.from_numpy(atom_feats_raw + 1).long().unsqueeze(0)
    atom_mask = torch.ones(1, n, dtype=torch.bool)
    coords_t = torch.from_numpy(coords).float().unsqueeze(0)
    has_3d = torch.ones(1, dtype=torch.bool)
    pos_enc = torch.zeros(1, n, max(pe_in_dim, 0), dtype=torch.float32)
    return _Batch(
        atom_features=atom_features.to(device),
        atom_mask=atom_mask.to(device),
        coords=coords_t.to(device),
        has_3d=has_3d.to(device),
        pos_enc=pos_enc.to(device),
    )


def embed_smiles(
    path_to_checkpoint: str,
    smiles_string: str,
    device: Optional[str] = None,
    n_confs: int = 10,
    mmff: bool = True,
    max_iters: int = 200,
) -> np.ndarray:
    """Run a Graphormer3D checkpoint on a single SMILES, return the graph
    embedding as a 1-D numpy array of shape (d_model,).

    Args:
        path_to_checkpoint: path to a Graphormer3D checkpoint (encoder_best.pt
            or best.pt).
        smiles_string: input SMILES.
        device: torch device string ("cpu", "cuda", "cuda:0", ...). Defaults
            to CUDA when available, otherwise CPU.
        n_confs, mmff, max_iters: 3D conformer search parameters; defaults
            match the training-time featurizer.
    """
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    dev = torch.device(device)

    ckpt = torch.load(path_to_checkpoint, map_location="cpu", weights_only=False)
    encoder_sd, enc_cfg = _extract_encoder_state_and_config(ckpt)

    model = Graphormer3D(enc_cfg)
    missing, unexpected = model.load_state_dict(encoder_sd, strict=False)
    if missing:
        raise RuntimeError(
            f"Missing keys when loading encoder state_dict: {missing[:5]}"
            f"{'...' if len(missing) > 5 else ''}"
        )
    # Unexpected keys are tolerated (e.g. non-persistent buffers); warn instead.
    if unexpected:
        print(f"[embed_smiles] ignoring {len(unexpected)} unexpected keys "
              f"(first few: {unexpected[:3]})", file=sys.stderr)
    model.to(dev).eval()

    atom_feats_raw, coords = featurize_smiles(
        smiles_string, n_confs=n_confs, mmff=mmff, max_iters=max_iters,
    )
    batch = _make_single_batch(atom_feats_raw, coords,
                               pe_in_dim=enc_cfg.pe_in_dim, device=dev)

    with torch.no_grad():
        graph = model(batch)
    return graph.squeeze(0).detach().cpu().numpy()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: embed_graphormer3d.py <checkpoint.pt> <SMILES>",
              file=sys.stderr)
        sys.exit(2)
    emb = embed_smiles(sys.argv[1], sys.argv[2])
    print(f"embedding shape: {emb.shape}")
    print(emb)
