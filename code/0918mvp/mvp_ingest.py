
# mvp_ingest.py
import pandas as pd
import numpy as np

# -------- Module 01 · Ingest & Normalize --------

def load_spectra(path: str) -> pd.DataFrame:
    df = pd.read_parquet(path) if path.endswith((".parquet", ".pq")) else pd.read_csv(path)
    # Minimal columns contract
    required = ["wiki_id", "rt", "precursor_mz"]
    for c in required:
        if c not in df.columns:
            raise ValueError(f"load_spectra: missing required column: {c}")
    # Fill optional columns
    if "assay" not in df.columns: df["assay"] = "unknown"
    if "polarity" not in df.columns: df["polarity"] = "unknown"
    if "is_manual_annotated" not in df.columns:
        # default False if not present
        df["is_manual_annotated"] = False
    # Coerce types
    df["wiki_id"] = df["wiki_id"].astype(str)
    df["rt"] = pd.to_numeric(df["rt"], errors="coerce")
    df["precursor_mz"] = pd.to_numeric(df["precursor_mz"], errors="coerce")
    return df[["wiki_id","rt","precursor_mz","assay","polarity","is_manual_annotated"]]

def load_hits(path: str) -> pd.DataFrame:
    df = pd.read_parquet(path) if path.endswith((".parquet", ".pq")) else pd.read_csv(path)
    # Try to detect/rename common columns
    rename_map = {}
    if "id" in df.columns and "library_id" not in df.columns:
        rename_map["id"] = "library_id"
    if "precursor_mz" in df.columns and "lib_precursor_mz" not in df.columns:
        rename_map["precursor_mz"] = "lib_precursor_mz"
    if "rt" in df.columns and "lib_rt" not in df.columns:
        rename_map["rt"] = "lib_rt"
    if rename_map:
        df = df.rename(columns=rename_map)
    # Minimal columns
    required = ["wiki_id", "library_id", "name", "adduct", "lib_precursor_mz", "entropy_similarity"]
    for c in required:
        if c not in df.columns:
            raise ValueError(f"load_hits: missing required column: {c}")
    # Optional lib_rt
    if "lib_rt" not in df.columns:
        df["lib_rt"] = np.nan
    # Types
    df["wiki_id"] = df["wiki_id"].astype(str)
    df["library_id"] = df["library_id"].astype(str)
    df["name"] = df["name"].astype(str)
    df["adduct"] = df["adduct"].astype(str)
    df["lib_precursor_mz"] = pd.to_numeric(df["lib_precursor_mz"], errors="coerce")
    df["entropy_similarity"] = pd.to_numeric(df["entropy_similarity"], errors="coerce")
    df["lib_rt"] = pd.to_numeric(df["lib_rt"], errors="coerce")
    # Keep minimal + passthrough db/library_type if present
    keep = ["wiki_id","library_id","name","adduct","lib_precursor_mz","lib_rt","entropy_similarity"]
    for extra in ["db","library_type","rank"]:
        if extra in df.columns: keep.append(extra)
    return df[keep]

# -------- Module 02 · Join & Confusable Sets --------

def join_hits_to_spectra(hits: pd.DataFrame, spectra: pd.DataFrame) -> pd.DataFrame:
    spec_cols = ["wiki_id","rt","precursor_mz","assay","polarity","is_manual_annotated"]
    missing = [c for c in spec_cols if c not in spectra.columns]
    if missing:
        raise ValueError(f"join_hits_to_spectra: spectra missing {missing}")
    joint = hits.merge(spectra[spec_cols], on="wiki_id", how="left", validate="m:1")
    # Compute Δppm
    joint["delta_ppm"] = joint["precursor_mz"] - joint["lib_precursor_mz"]
    # ΔRT and zRT placeholder; filled in build_confusable_sets
    return joint

def build_confusable_sets(joint: pd.DataFrame, ms2_min: float = 0.75, zrt_k: float = 3.0, rt_sigma_sec: float = 2.0) -> pd.DataFrame:
    df = joint.copy()
    # Filter by  entropy similarity
    df = df[df["entropy_similarity"].astype(float) >= ms2_min].copy()
    # ΔRT and zRT if lib_rt is available
    df["delta_rt"] = np.where(df["lib_rt"].notna(), df["rt"] - df["lib_rt"], np.nan)
    df["zrt"] = df["delta_rt"] / rt_sigma_sec
    # Keep hits with |zRT|<=k when zrt is finite; let NaN pass (flag later)
    mask_zrt = df["zrt"].abs() <= zrt_k
    df = pd.concat([df[df["zrt"].isna()], df[mask_zrt]], axis=0)
    # Flags
    df["rt_missing"] = df["lib_rt"].isna()
    return df.reset_index(drop=True)

def collapse_duplicates_by_structure(cset: pd.DataFrame, key_cols=("library_id",)) -> pd.DataFrame:
    # In MVP, use library_id as structure key proxy (replace with InChIKey when available)
    key_cols = list(key_cols)
    group_cols = ["wiki_id"] + key_cols
    agg = {
        "name":"first",
        "adduct": lambda s: ";".join(sorted(pd.Series(s, dtype=str).dropna().unique())),
        "lib_precursor_mz":"median",
        "lib_rt":"median",
        "entropy_similarity":"max",
        "delta_ppm":"median",
        "delta_rt":"median",
        "zrt":"median",
        "rt_missing":"min",
        "is_manual_annotated":"max",
        "precursor_mz":"median",
    }
    extra_cols = [c for c in ["db","library_type","rank"] if c in cset.columns]
    for c in extra_cols: agg[c] = "first"
    agg["member_hits"] = ("library_id", "count")
    df = cset.groupby(group_cols, as_index=False).agg(agg)
    # Fix multiindex output if any
    if isinstance(df.columns, pd.MultiIndex):
        df.columns = ["_".join([str(x) for x in col if x]) for col in df.columns.values]
    if "member_hits" not in df.columns and ("library_id_count" in df.columns):
        df = df.rename(columns={"library_id_count":"member_hits"})
    return df
