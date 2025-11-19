#!/usr/bin/env python3
"""
MAVE Ingestion Module - Phase B: Load and Parse Raw MaveDB Data

This module handles loading raw MaveDB CSV/TSV exports and converting them to
a standardized canonical schema for downstream processing.

Key Responsibilities:
  - Load raw MaveDB CSV/TSV files
  - Auto-detect HGVS variant and score columns
  - Parse HGVS notation (protein: p.Met1Val, nucleotide, etc.)
  - Normalize amino acids to one-letter codes
  - Filter to single amino acid substitutions
  - Map to canonical schema (dataset_id, gene, pos, wt_aa, mut_aa, score, etc.)
  - Output: {gene}_{dataset}_raw.parquet

Extensibility:
  - Add new HGVS formats by extending parse_hgvs_variant()
  - Add new column detection in ingest_dataset() hgvs_candidates/score_candidates
  - Support new file formats by adding elif in read logic

Usage:
  from src.mave.pipeline.ingest import ingest_all_datasets
  results = ingest_all_datasets(Path("config/mave_datasets.yaml"))

MaveDB ingestion: Load raw exports and convert to canonical schema.

This module handles:
1. Loading raw TSV/CSV exports from MaveDB
2. Mapping raw columns to canonical schema fields
3. Normalizing amino acids to one-letter codes
4. Handling multi-condition aggregation
5. Validating data integrity
"""

import logging
import json
from pathlib import Path
from typing import Dict, Optional, Any
import pandas as pd
import numpy as np
from src.mave.utilities.mavedb_loader import ingest_mavedb_csv, find_mavedb_files

logger = logging.getLogger(__name__)

# Amino acid three-letter to one-letter mapping
AA_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "STOP": "*", "TER": "*", "*": "*",
}


def parse_hgvs_variant(hgvs_str: str) -> Optional[Dict[str, Any]]:
    """
    Parse HGVS-style variant string.
    
    Handles formats like:
    - p.Met1775Lys (three-letter aa codes)
    - p.M1775K (one-letter aa codes)
    - M1775K (no 'p.' prefix)
    - 1775M>K (position:wt>mut)
    """
    if not hgvs_str or not isinstance(hgvs_str, str):
        return None
    
    hgvs_str = hgvs_str.strip()
    
    # Remove 'p.' prefix if present
    if hgvs_str.startswith("p."):
        hgvs_str = hgvs_str[2:]
    
    # Format: position:wt>mut or positionWT>MUT or positionWT>MUT
    if ">" in hgvs_str:
        import re
        
        # Try position:wt>mut first
        if ":" in hgvs_str:
            parts = hgvs_str.split(":")
            if len(parts) == 2:
                try:
                    pos = int(parts[0])
                    aa_parts = parts[1].split(">")
                    if len(aa_parts) == 2:
                        wt_aa = aa_parts[0].upper()
                        mut_aa = aa_parts[1].upper()
                        return {"pos": pos, "wt_aa": wt_aa, "mut_aa": mut_aa}
                except ValueError:
                    pass
        
        # Try positionWT>MUT format like "1775M>K"
        match = re.match(r"(\d+)([A-Z*]+)>([A-Z*]+)", hgvs_str)
        if match:
            pos_str, wt_aa, mut_aa = match.groups()
            try:
                pos = int(pos_str)
                return {
                    "pos": pos,
                    "wt_aa": wt_aa.upper(),
                    "mut_aa": mut_aa.upper(),
                }
            except ValueError:
                pass
    
    # Format: three-letter codes like Met1775Lys or M1775K
    import re
    match = re.match(r"([A-Za-z*]+)(\d+)([A-Za-z*]+)", hgvs_str)
    if match:
        wt_str, pos_str, mut_str = match.groups()
        try:
            pos = int(pos_str)
            # Convert if three-letter
            wt_aa = AA_THREE_TO_ONE.get(wt_str.upper(), wt_str.upper())
            mut_aa = AA_THREE_TO_ONE.get(mut_str.upper(), mut_str.upper())
            
            # Only return if we got single characters
            if len(wt_aa) == 1 and len(mut_aa) == 1:
                return {
                    "pos": pos,
                    "wt_aa": wt_aa.upper(),
                    "mut_aa": mut_aa.upper(),
                }
        except (ValueError, KeyError):
            pass
    
    return None


def ingest_dataset(
    dataset_cfg: Dict[str, Any],
    raw_path: Path,
    output_dir: Path,
    gene_config: Optional[Dict[str, Any]] = None,
) -> Optional[pd.DataFrame]:
    """
    Ingest a single MaveDB dataset and write to canonical schema.
    
    Args:
        dataset_cfg: Configuration dict with keys:
            - dataset_id, gene, assay_type, mavedb_accession
            - raw_score_col: Which column contains the raw functional score
            - hgvs_col: Which column contains variant notation
        raw_path: Path to raw MaveDB export (TSV/CSV)
        output_dir: Directory to write canonical parquet
        gene_config: Optional gene-specific config
    
    Returns:
        DataFrame in canonical schema, or None if ingestion failed
    """
    dataset_id = dataset_cfg.get("dataset_id")
    gene = dataset_cfg.get("gene")
    mavedb_accession = dataset_cfg.get("mavedb_accession")
    assay_type = dataset_cfg.get("assay_type")
    
    logger.info(f"Ingesting {dataset_id} from {raw_path.name}")
    
    if not raw_path.exists():
        logger.error(f"Raw file not found: {raw_path}")
        return None
    
    try:
        # Detect delimiter
        if raw_path.suffix == ".csv":
            df_raw = pd.read_csv(raw_path)
        else:  # .tsv or others
            df_raw = pd.read_csv(raw_path, sep="\t")
        
        logger.info(f"Loaded {len(df_raw)} rows from raw file")
    except Exception as e:
        logger.error(f"Failed to read {raw_path}: {e}")
        return None
    
    # Identify columns
    # Common column name patterns - PREFER protein notation
    hgvs_candidates = ["hgvs_pro", "hgvs_nt", "hgvs", "variant", "protein_variant"]
    score_candidates = ["score", "mean_score", "median_score", "functional_score", "effect"]
    
    hgvs_col = None
    for cand in hgvs_candidates:
        if cand in df_raw.columns:
            hgvs_col = cand
            break
    
    score_col = None
    for cand in score_candidates:
        if cand in df_raw.columns:
            score_col = cand
            break
    
    if not hgvs_col or not score_col:
        logger.error(
            f"Could not identify HGVS column ({hgvs_col}) or score column ({score_col}). "
            f"Available columns: {list(df_raw.columns)}"
        )
        return None
    
    logger.info(f"Using HGVS column: {hgvs_col}, Score column: {score_col}")
    
    # Parse variants
    variants = []
    unparseable = 0
    
    for idx, row in df_raw.iterrows():
        hgvs_str = str(row.get(hgvs_col, "")).strip()
        score_val = row.get(score_col)
        
        # Skip missing scores
        if pd.isna(score_val):
            continue
        
        # Parse HGVS
        parsed = parse_hgvs_variant(hgvs_str)
        if not parsed:
            unparseable += 1
            continue
        
        wt_aa = parsed["wt_aa"]
        pos = parsed["pos"]
        mut_aa = parsed["mut_aa"]
        
        # Skip if WT == MUT (synonymous, not useful)
        if wt_aa == mut_aa:
            continue
        
        # Build variant_id
        variant_id = f"{gene}:p.{wt_aa}{pos}{mut_aa}"
        
        # Collect metadata if available
        metadata = {}
        for col in df_raw.columns:
            if col not in [hgvs_col, score_col] and pd.notna(row.get(col)):
                val = row.get(col)
                # Only store simple types
                if isinstance(val, (str, int, float, bool)):
                    metadata[col] = val
        
        variants.append({
            "dataset_id": dataset_id,
            "gene": gene,
            "wt_aa": wt_aa,
            "pos": int(pos),
            "mut_aa": mut_aa,
            "variant_id": variant_id,
            "functional_score_raw": float(score_val),
            "functional_score": float(score_val),  # Will be normalized in next phase
            "assay_type": assay_type,
            "mavedb_accession": mavedb_accession,
            "metadata_json": json.dumps(metadata) if metadata else None,
        })
    
    if not variants:
        logger.error(f"No variants extracted from {dataset_id}")
        return None
    
    df_ingested = pd.DataFrame(variants)
    
    # Stats and validation
    n_raw = len(df_raw)
    n_single_aa_subs = len(df_ingested)
    unique_combos = len(df_ingested[["wt_aa", "pos", "mut_aa"]].drop_duplicates())
    
    logger.info(f"  Raw rows: {n_raw}")
    logger.info(f"  Single AA substitutions: {n_single_aa_subs}")
    logger.info(f"  Unique (WT, pos, MUT) combinations: {unique_combos}")
    logger.info(f"  Unparseable variants: {unparseable}")
    
    # Write output
    output_dir.mkdir(parents=True, exist_ok=True)
    output_parquet = output_dir / f"{dataset_id}_raw.parquet"
    
    try:
        df_ingested.to_parquet(output_parquet, index=False)
        logger.info(f"Wrote {len(df_ingested)} rows to {output_parquet}")
    except Exception as e:
        logger.error(f"Failed to write {output_parquet}: {e}")
        return None
    
    return df_ingested


def ingest_mavedb_dump_auto(
    dump_dir: Path = Path("data_raw/mave/mavedb-dump.20250612164404"),
    limit: int = 5,
) -> str:
    """
    Auto-detect MaveDB files from dump and generate config.
    
    Returns YAML config string with datasets found.
    """
    from .mavedb_loader import auto_generate_config
    
    yaml_config = auto_generate_config(dump_dir, limit=limit)
    logger.info(f"Generated config for {limit} datasets from MaveDB dump")
    
    return yaml_config


def ingest_all_datasets(
    datasets_config_path: Path,
    output_dir: Path = Path("data_processed/mave"),
) -> Dict[str, pd.DataFrame]:
    """
    Ingest all datasets specified in config.
    
    Args:
        datasets_config_path: Path to mave_datasets.yaml
        output_dir: Base output directory
    
    Returns:
        Dict mapping dataset_id to DataFrame
    """
    import yaml
    
    if not datasets_config_path.exists():
        logger.error(f"Config not found: {datasets_config_path}")
        return {}
    
    try:
        with open(datasets_config_path) as f:
            config = yaml.safe_load(f)
    except Exception as e:
        logger.error(f"Failed to load config: {e}")
        return {}
    
    datasets = config.get("datasets", [])
    if not datasets:
        logger.warning("No datasets in config")
        return {}
    
    results = {}
    for dataset_cfg in datasets:
        dataset_id = dataset_cfg.get("dataset_id")
        if not dataset_id:
            logger.warning("Skipping dataset with no dataset_id")
            continue
        
        raw_path = Path(dataset_cfg.get("raw_data_path", ""))
        if not raw_path.is_absolute():
            # If path doesn't start with data_raw, prepend it
            if not str(raw_path).startswith("data_raw"):
                raw_path = Path("data_raw/mave") / raw_path
            else:
                # Otherwise use as-is (might already have data_raw in path)
                raw_path = Path(raw_path)
        
        dataset_output_dir = output_dir / dataset_id
        
        df = ingest_dataset(
            dataset_cfg,
            raw_path,
            dataset_output_dir,
        )
        
        if df is not None:
            results[dataset_id] = df
    
    return results

