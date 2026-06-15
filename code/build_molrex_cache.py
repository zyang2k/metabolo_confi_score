"""Stage 2: build IK14 -> MolRex (Graphormer3D) embedding cache.

Multiprocessing pool, each worker loads the encoder checkpoint once and
embeds SMILES from `data/ik14_to_smiles_labeled.csv` in chunks.

Output:
  data/molrex_embeddings.npz   keys: ik14 (object array), embedding (N, 512)
  data/molrex_failed_ik14.csv  IK14s where conformer search or forward failed
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from embed_graphormer3d import (  # noqa: E402
    Graphormer3D,
    _extract_encoder_state_and_config,
    _make_single_batch,
    featurize_smiles,
)

CKPT_PATH = "data/encoder_best.pt"
SMILES_CSV = "data/ik14_to_smiles_labeled.csv"
OUT_NPZ = "data/molrex_embeddings.npz"
FAILED_CSV = "data/molrex_failed_ik14.csv"

_MODEL: Graphormer3D | None = None
_PE_IN_DIM: int = 0
_DEVICE = torch.device("cpu")


def _init_worker(ckpt_path: str) -> None:
    """Load encoder once per worker; reused across all tasks in this process."""
    global _MODEL, _PE_IN_DIM
    torch.set_num_threads(1)  # avoid oversubscription across workers
    ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    sd, cfg = _extract_encoder_state_and_config(ckpt)
    model = Graphormer3D(cfg)
    missing, _ = model.load_state_dict(sd, strict=False)
    if missing:
        raise RuntimeError(f"missing keys when loading encoder: {missing[:5]}")
    model.to(_DEVICE).eval()
    _MODEL = model
    _PE_IN_DIM = cfg.pe_in_dim


def _embed_one(args: tuple[str, str]) -> tuple[str, np.ndarray | None, str | None]:
    ik14, smi = args
    try:
        feats, coords = featurize_smiles(smi)
        batch = _make_single_batch(feats, coords, _PE_IN_DIM, _DEVICE)
        with torch.no_grad():
            graph = _MODEL(batch)  # type: ignore[misc]
        return ik14, graph.squeeze(0).numpy().astype(np.float32), None
    except Exception as e:  # conformer fail, SMILES parse fail, etc.
        return ik14, None, type(e).__name__ + ": " + str(e)[:200]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--limit", type=int, default=None,
                    help="embed only the first N SMILES (smoke / debug).")
    ap.add_argument("--chunksize", type=int, default=20,
                    help="Pool.imap_unordered chunksize.")
    args = ap.parse_args()

    df = pd.read_csv(SMILES_CSV)
    if args.limit is not None:
        df = df.head(args.limit)
    tasks = list(zip(df.hit_ik14, df.smiles))
    print(f"embedding {len(tasks):,} SMILES with {args.workers} workers "
          f"(chunksize={args.chunksize}) ...", flush=True)

    embs = np.zeros((len(tasks), 512), dtype=np.float32)
    ik14s_out: list[str] = []
    failed: list[dict] = []

    t0 = time.time()
    log_every = max(500, len(tasks) // 50)
    with mp.Pool(args.workers, initializer=_init_worker,
                 initargs=(CKPT_PATH,)) as pool:
        for i, (ik14, emb, err) in enumerate(
            pool.imap_unordered(_embed_one, tasks, chunksize=args.chunksize)
        ):
            if err is None and emb is not None:
                embs[len(ik14s_out)] = emb
                ik14s_out.append(ik14)
            else:
                failed.append({"hit_ik14": ik14, "error": err})
            if (i + 1) % log_every == 0 or (i + 1) == len(tasks):
                elapsed = time.time() - t0
                rate = (i + 1) / elapsed
                eta = (len(tasks) - i - 1) / max(rate, 1e-6)
                print(f"  [{i+1:>6,d}/{len(tasks):,}]  "
                      f"ok={len(ik14s_out):>6,d}  "
                      f"failed={len(failed):>4,d}  "
                      f"rate={rate:5.1f}/s  "
                      f"eta={eta/60:5.1f} min", flush=True)

    embs = embs[: len(ik14s_out)]
    ik14_arr = np.array(ik14s_out, dtype=object)
    np.savez_compressed(OUT_NPZ, ik14=ik14_arr, embedding=embs)
    print(f"\nwrote {OUT_NPZ}  shape={embs.shape}  "
          f"size={Path(OUT_NPZ).stat().st_size/1e6:.1f} MB", flush=True)

    pd.DataFrame(failed).to_csv(FAILED_CSV, index=False)
    print(f"wrote {FAILED_CSV}  ({len(failed):,} failures)", flush=True)

    total = time.time() - t0
    print(f"\ntotal: {total/60:.1f} min   "
          f"mean rate={len(tasks)/total:.1f} mol/s   "
          f"success: {len(ik14s_out)/len(tasks):.1%}", flush=True)


if __name__ == "__main__":
    main()
