"""Stage-0 smoke test for MolRex (Graphormer3D) encoder.

Confirms:
  - checkpoint loads + forward pass runs on CPU
  - embedding shape is (512,) and finite
  - per-call wall-clock (projects full-cache cost)
  - structural-sanity ordering: cos(caffeine, theobromine) > cos(caffeine, glucose)
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
# --- sibling-import path shim (code/ root) ---
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))

from embed_graphormer3d import embed_smiles  # noqa: E402

CKPT = "data/encoder_best.pt"

PROBES = {
    "caffeine":    "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "theobromine": "Cn1cnc2[nH]c(=O)n(C)c(=O)c12",  # xanthine-family, sibling of caffeine
    "glucose":     "OCC1OC(O)C(O)C(O)C1O",
    "atp":         "Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O",
    "citrate":     "OC(=O)CC(O)(CC(=O)O)C(=O)O",
    "pyridoxine":  "Cc1ncc(CO)c(CO)c1O",
}


def cos(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-12))


def main() -> None:
    print(f"loading checkpoint: {CKPT}")
    embs: dict[str, np.ndarray] = {}
    times: list[float] = []

    for name, smi in PROBES.items():
        t0 = time.time()
        emb = embed_smiles(CKPT, smi, device="cpu")
        dt = time.time() - t0
        times.append(dt)
        finite = bool(np.isfinite(emb).all())
        print(
            f"  {name:12s}  shape={emb.shape}  "
            f"finite={finite}  norm={np.linalg.norm(emb):8.2f}  "
            f"t={dt:5.2f}s"
        )
        embs[name] = emb

    print(f"\nmean per-call: {np.mean(times):.2f}s   "
          f"median: {np.median(times):.2f}s   "
          f"max: {np.max(times):.2f}s")

    print("\nstructural sanity (cosine):")
    pairs = [
        ("caffeine", "theobromine"),  # expect HIGH (xanthine siblings)
        ("caffeine", "glucose"),      # expect LOW
        ("caffeine", "atp"),          # expect MEDIUM (shares purine ring)
        ("glucose", "citrate"),       # expect LOW-MED (both small polar)
        ("atp", "pyridoxine"),        # expect LOW
    ]
    for a, b in pairs:
        print(f"  cos({a:11s}, {b:11s}) = {cos(embs[a], embs[b]):+.4f}")

    caff_theo = cos(embs["caffeine"], embs["theobromine"])
    caff_glc = cos(embs["caffeine"], embs["glucose"])
    print(
        f"\nordering check: cos(caff,theo)={caff_theo:+.4f} "
        f"{'>' if caff_theo > caff_glc else '<='} "
        f"cos(caff,glc)={caff_glc:+.4f}  "
        f"-> {'PASS' if caff_theo > caff_glc else 'FAIL'}"
    )


if __name__ == "__main__":
    main()
