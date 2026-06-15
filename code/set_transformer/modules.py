"""Set Transformer building blocks (Lee et al., ICML 2019), masked.

Implements:
  - MAB(X, Y)  : Multihead Attention Block
  - SAB(X)     : MAB(X, X), full self-attention   O(n^2)
  - ISAB(X)    : induced-set attention via m inducing points   O(mn)
  - PMA(X, k)  : pooling by multihead attention to k seed vectors

Mask convention:
  mask[b, i] is True for valid set elements, False for padding.
  Attention scores at padded keys are set to -1e9 before softmax.

We implement multi-head attention manually rather than relying on
nn.MultiheadAttention so masking and shapes stay transparent for a small
research model.
"""

from __future__ import annotations

import math
from typing import Optional

import torch
import torch.nn as nn
import torch.nn.functional as F


def _split_heads(x: torch.Tensor, n_heads: int) -> torch.Tensor:
    """(B, N, D) -> (B, H, N, D/H)."""
    b, n, d = x.shape
    assert d % n_heads == 0, f'd={d} not divisible by n_heads={n_heads}'
    return x.view(b, n, n_heads, d // n_heads).transpose(1, 2)


def _merge_heads(x: torch.Tensor) -> torch.Tensor:
    """(B, H, N, D/H) -> (B, N, D)."""
    b, h, n, dh = x.shape
    return x.transpose(1, 2).contiguous().view(b, n, h * dh)


class MAB(nn.Module):
    """Multihead Attention Block: cross-attention from X to Y with residual + FFN."""

    def __init__(self, d_q: int, d_kv: int, d_model: int, n_heads: int,
                 dropout: float = 0.0, ln: bool = True):
        super().__init__()
        self.n_heads = n_heads
        self.d_model = d_model
        self.q_proj = nn.Linear(d_q, d_model)
        self.k_proj = nn.Linear(d_kv, d_model)
        self.v_proj = nn.Linear(d_kv, d_model)
        self.o_proj = nn.Linear(d_model, d_model)
        self.ffn = nn.Sequential(
            nn.Linear(d_model, d_model * 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(d_model * 2, d_model),
        )
        self.dropout = nn.Dropout(dropout)
        self.ln1 = nn.LayerNorm(d_model) if ln else nn.Identity()
        self.ln2 = nn.LayerNorm(d_model) if ln else nn.Identity()

    def forward(self, X: torch.Tensor, Y: torch.Tensor,
                mask_q: Optional[torch.Tensor] = None,
                mask_k: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        X: (B, Nq, d_q)   queries
        Y: (B, Nk, d_kv)  keys/values
        mask_q: (B, Nq) bool, True = valid (only used to zero output rows)
        mask_k: (B, Nk) bool, True = valid (used inside softmax)
        """
        Xp = self.q_proj(X)                              # (B, Nq, d_model)
        Q = _split_heads(Xp, self.n_heads)               # (B, H, Nq, d/h)
        K = _split_heads(self.k_proj(Y), self.n_heads)   # (B, H, Nk, d/h)
        V = _split_heads(self.v_proj(Y), self.n_heads)

        scale = 1.0 / math.sqrt(Q.shape[-1])
        scores = torch.matmul(Q, K.transpose(-1, -2)) * scale  # (B, H, Nq, Nk)
        if mask_k is not None:
            scores = scores.masked_fill(~mask_k[:, None, None, :], -1e9)

        attn = F.softmax(scores, dim=-1)
        attn = self.dropout(attn)
        out = torch.matmul(attn, V)                       # (B, H, Nq, d/h)
        out = _merge_heads(out)                           # (B, Nq, d_model)
        out = self.o_proj(out)

        H = self.ln1(Xp + out)                            # residual around attention
        H = self.ln2(H + self.ffn(H))                     # residual around FFN

        if mask_q is not None:
            H = H * mask_q[..., None].to(H.dtype)
        return H


class SAB(nn.Module):
    """Self-Attention Block: MAB(X, X). O(n^2)."""

    def __init__(self, d_model: int, n_heads: int, dropout: float = 0.0):
        super().__init__()
        self.mab = MAB(d_model, d_model, d_model, n_heads, dropout)

    def forward(self, X: torch.Tensor, mask: Optional[torch.Tensor] = None) -> torch.Tensor:
        return self.mab(X, X, mask_q=mask, mask_k=mask)


class ISAB(nn.Module):
    """Induced-Set Attention Block. m learned inducing points; O(mn).

    Two MABs:
      H = MAB(I, X)   inducing points attend to X
      Z = MAB(X, H)   X attends back to inducing points
    """

    def __init__(self, d_in: int, d_model: int, n_heads: int, n_inducing: int,
                 dropout: float = 0.0):
        super().__init__()
        self.I = nn.Parameter(torch.randn(1, n_inducing, d_model) * 0.02)
        self.mab1 = MAB(d_model, d_in, d_model, n_heads, dropout)   # I attends to X
        self.mab2 = MAB(d_in, d_model, d_model, n_heads, dropout)   # X attends to H

    def forward(self, X: torch.Tensor, mask: Optional[torch.Tensor] = None) -> torch.Tensor:
        b = X.shape[0]
        I = self.I.expand(b, -1, -1)
        # Inducing points have no padding; their mask is all True.
        H = self.mab1(I, X, mask_q=None, mask_k=mask)               # (B, m, d_model)
        Z = self.mab2(X, H, mask_q=mask, mask_k=None)               # (B, N, d_model)
        return Z


class PMA(nn.Module):
    """Pooling by Multihead Attention: k seed vectors attend to X."""

    def __init__(self, d_model: int, n_heads: int, k: int = 1, dropout: float = 0.0):
        super().__init__()
        self.S = nn.Parameter(torch.randn(1, k, d_model) * 0.02)
        self.mab = MAB(d_model, d_model, d_model, n_heads, dropout)

    def forward(self, X: torch.Tensor, mask: Optional[torch.Tensor] = None) -> torch.Tensor:
        b = X.shape[0]
        S = self.S.expand(b, -1, -1)
        return self.mab(S, X, mask_q=None, mask_k=mask)             # (B, k, d_model)


__all__ = ['MAB', 'SAB', 'ISAB', 'PMA']
