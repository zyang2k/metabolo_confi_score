"""EvidenceSetTransformer: per-spectrum candidate-set scorer.

Three output heads, all returned in one forward pass so we can compare them
against each other and against the GBM baseline:

  (A) equiv_logits  (B, K)         per-candidate P(TP), permutation equivariant
  (B) set_logit     (B,)           single P(any candidate is TP), invariant
  (C) z_set, Z_top  (B, d), (B, d) embeddings for the hybrid GBM head

Categorical features are embedded and concatenated with numeric features before
the encoder, mirroring how the GBM uses one-hot expansions.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, NamedTuple

import torch
import torch.nn as nn

from .modules import ISAB, PMA


class STOutput(NamedTuple):
    equiv_logits: torch.Tensor   # (B, K)  per-candidate logit; padding = -1e9
    set_logit: torch.Tensor      # (B,)    set-level logit
    z_set: torch.Tensor          # (B, d_model)  pooled set embedding
    z_top: torch.Tensor          # (B, d_model)  encoder output at rank-0 (top-by-esim)


@dataclass
class STConfig:
    d_num: int                       # number of numeric features
    cat_cardinalities: List[int]     # vocab size per categorical (including pad/missing)
    cat_emb_dim: int = 8
    d_model: int = 64
    n_heads: int = 4
    n_inducing: int = 8
    n_isab: int = 2
    dropout: float = 0.1


class EvidenceSetTransformer(nn.Module):
    def __init__(self, cfg: STConfig):
        super().__init__()
        self.cfg = cfg
        self.cat_emb = nn.ModuleList([
            nn.Embedding(card, cfg.cat_emb_dim, padding_idx=0)
            for card in cfg.cat_cardinalities
        ])
        d_in = cfg.d_num + cfg.cat_emb_dim * len(cfg.cat_cardinalities)
        self.input_proj = nn.Sequential(
            nn.Linear(d_in, cfg.d_model),
            nn.LayerNorm(cfg.d_model),
        )
        self.encoder = nn.ModuleList([
            ISAB(cfg.d_model, cfg.d_model, cfg.n_heads, cfg.n_inducing, cfg.dropout)
            for _ in range(cfg.n_isab)
        ])
        self.equiv_head = nn.Sequential(
            nn.Linear(cfg.d_model, cfg.d_model),
            nn.GELU(),
            nn.Dropout(cfg.dropout),
            nn.Linear(cfg.d_model, 1),
        )
        self.pma = PMA(cfg.d_model, cfg.n_heads, k=1, dropout=cfg.dropout)
        self.set_head = nn.Sequential(
            nn.Linear(cfg.d_model, cfg.d_model),
            nn.GELU(),
            nn.Dropout(cfg.dropout),
            nn.Linear(cfg.d_model, 1),
        )

    def forward(self,
                X: torch.Tensor,                 # (B, K, d_num)
                C: torch.Tensor,                 # (B, K, n_cat)
                M: torch.Tensor                  # (B, K) bool
                ) -> STOutput:
        emb = [self.cat_emb[i](C[..., i]) for i in range(C.shape[-1])]
        H = torch.cat([X] + emb, dim=-1)
        H = self.input_proj(H)
        # Zero out padded rows so nothing leaks through skip connections
        H = H * M[..., None].to(H.dtype)
        for blk in self.encoder:
            H = blk(H, mask=M)

        equiv_logits = self.equiv_head(H).squeeze(-1)              # (B, K)
        equiv_logits = equiv_logits.masked_fill(~M, -1e9)

        z_set = self.pma(H, mask=M).squeeze(1)                      # (B, d_model)
        set_logit = self.set_head(z_set).squeeze(-1)                # (B,)

        # Top-by-esim is rank 0 (dataset is sorted desc on entropy_similarity).
        z_top = H[:, 0, :]                                          # (B, d_model)

        return STOutput(equiv_logits=equiv_logits,
                        set_logit=set_logit,
                        z_set=z_set,
                        z_top=z_top)


def masked_bce_per_hit(logits: torch.Tensor,
                       y_hit: torch.Tensor,
                       mask: torch.Tensor,
                       pos_weight: float = 1.0) -> torch.Tensor:
    """BCE on per-candidate logits; ignores padding (y_hit == -1).

    Args:
        logits: (B, K) raw logits (use the un-masked logits, not -1e9 padded).
        y_hit:  (B, K) {0, 1, -1}
        mask:   (B, K) bool, True for valid candidates
        pos_weight: scalar for class balancing on the positive class
    """
    valid = mask & (y_hit >= 0)
    if valid.sum() == 0:
        return logits.new_zeros(())
    y = y_hit.clamp(min=0).to(logits.dtype)
    pw = torch.tensor(pos_weight, device=logits.device, dtype=logits.dtype)
    loss = nn.functional.binary_cross_entropy_with_logits(
        logits, y, pos_weight=pw, reduction='none'
    )
    return (loss * valid.to(loss.dtype)).sum() / valid.to(loss.dtype).sum()


__all__ = ['EvidenceSetTransformer', 'STConfig', 'STOutput', 'masked_bce_per_hit']
