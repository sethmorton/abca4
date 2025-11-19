#!/usr/bin/env python3
"""
Sanity checks and guardrails for MAVE data processing.

These checks help catch bad ingestions or misdefined hit sets early.
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def check_ingested_dataset(
    dataset_id: str,
    raw_parquet: Path,
) -> Tuple[bool, List[str]]:
    """
    Quick sanity checks on raw ingested data.
    
    Args:
        dataset_id: Dataset identifier
        raw_parquet: Path to {dataset_id}_raw.parquet
    
    Returns:
        (passed: bool, issues: List[str])
    """
    issues = []
    
    if not raw_parquet.exists():
        return False, [f"File not found: {raw_parquet}"]
    
    try:
        df = pd.read_parquet(raw_parquet)
    except Exception as e:
        return False, [f"Failed to load parquet: {e}"]
    
    # Check required columns
    required_cols = [
        "dataset_id", "gene", "wt_aa", "pos", "mut_aa", "variant_id",
        "functional_score_raw", "assay_type", "mavedb_accession"
    ]
    
    for col in required_cols:
        if col not in df.columns:
            issues.append(f"Missing required column: {col}")
    
    # Check for nulls in critical columns
    for col in ["gene", "wt_aa", "pos", "mut_aa", "variant_id", "functional_score_raw"]:
        if col in df.columns and df[col].isna().any():
            null_count = df[col].isna().sum()
            issues.append(f"Column {col} has {null_count} nulls")
    
    # Check data types
    if "pos" in df.columns and not pd.api.types.is_integer_dtype(df["pos"]):
        issues.append(f"pos column should be integer, got {df['pos'].dtype}")
    
    if "functional_score_raw" in df.columns and not pd.api.types.is_numeric_dtype(df["functional_score_raw"]):
        issues.append(f"functional_score_raw should be numeric, got {df['functional_score_raw'].dtype}")
    
    # Check amino acids are single characters
    if "wt_aa" in df.columns:
        bad_wt = df[df["wt_aa"].str.len() != 1]
        if not bad_wt.empty:
            issues.append(f"wt_aa has {len(bad_wt)} entries with != 1 character")
    
    if "mut_aa" in df.columns:
        bad_mut = df[df["mut_aa"].str.len() != 1]
        if not bad_mut.empty:
            issues.append(f"mut_aa has {len(bad_mut)} entries with != 1 character")
    
    # Check no synonymous (wt == mut)
    if "wt_aa" in df.columns and "mut_aa" in df.columns:
        synonymous = (df["wt_aa"] == df["mut_aa"]).sum()
        if synonymous > 0:
            issues.append(f"Found {synonymous} synonymous variants (wt_aa == mut_aa)")
    
    # Check functional scores are reasonable
    if "functional_score_raw" in df.columns:
        n_infinite = np.isinf(df["functional_score_raw"]).sum()
        if n_infinite > 0:
            issues.append(f"Found {n_infinite} infinite values in functional_score_raw")
        
        if df["functional_score_raw"].empty:
            issues.append("functional_score_raw is empty")
    
    passed = len(issues) == 0
    
    if passed:
        logger.info(f"✅ {dataset_id} ingestion checks passed")
    else:
        logger.error(f"❌ {dataset_id} ingestion checks found {len(issues)} issues:")
        for issue in issues:
            logger.error(f"   - {issue}")
    
    return passed, issues


def check_normalized_dataset(
    dataset_id: str,
    normalized_parquet: Path,
    config: Dict,
) -> Tuple[bool, List[str]]:
    """
    Sanity checks on normalized dataset.
    
    Args:
        dataset_id: Dataset identifier
        normalized_parquet: Path to {dataset_id}_normalized.parquet
        config: Dataset config from mave_datasets.yaml
    
    Returns:
        (passed: bool, issues: List[str])
    """
    issues = []
    
    if not normalized_parquet.exists():
        return False, [f"File not found: {normalized_parquet}"]
    
    try:
        df = pd.read_parquet(normalized_parquet)
    except Exception as e:
        return False, [f"Failed to load parquet: {e}"]
    
    # Check required columns added during normalization
    for col in ["functional_score", "is_hit", "hit_label"]:
        if col not in df.columns:
            issues.append(f"Missing expected column: {col}")
    
    # Check functional_score is normalized properly
    if "functional_score" in df.columns:
        # Depending on method, should be roughly in [-4, 4] for zscore, or [0, 1] for quantile
        min_score = df["functional_score"].min()
        max_score = df["functional_score"].max()
        
        norm_method = config.get("normalization", {}).get("method", "identity")
        
        if norm_method == "quantile":
            if min_score < 0 or max_score > 1:
                issues.append(
                    f"Quantile normalization should be [0, 1], got [{min_score:.4f}, {max_score:.4f}]"
                )
        elif norm_method == "minmax":
            if min_score < 0 or max_score > 1:
                issues.append(
                    f"Minmax normalization should be ~[0, 1], got [{min_score:.4f}, {max_score:.4f}]"
                )
    
    # Check hit fraction matches expectation
    if "is_hit" in df.columns:
        hit_rate = df["is_hit"].mean()
        expected_threshold = config.get("hit_definition", {}).get("percentile_threshold", 0.2)
        
        deviation = abs(hit_rate - expected_threshold)
        if deviation > 0.1:  # 10% tolerance for ties, truncation
            issues.append(
                f"Hit rate {hit_rate:.1%} deviates from expected {expected_threshold:.1%} "
                f"by {deviation:.1%} (may be due to ties)"
            )
    
    # Check hit labels are sensible
    if "hit_label" in df.columns:
        unique_labels = set(df["hit_label"].dropna().unique())
        valid_labels = {"lof", "non_lof", "gof", "non_gof", "unknown"}
        invalid = unique_labels - valid_labels
        if invalid:
            issues.append(f"Unexpected hit labels: {invalid}")
    
    passed = len(issues) == 0
    
    if passed:
        logger.info(f"✅ {dataset_id} normalization checks passed")
    else:
        logger.error(f"❌ {dataset_id} normalization checks found {len(issues)} issues:")
        for issue in issues:
            logger.error(f"   - {issue}")
    
    return passed, issues


def check_full_dataset(
    dataset_id: str,
    full_parquet: Path,
) -> Tuple[bool, List[str]]:
    """
    Sanity checks on full dataset (with features + functional data).
    
    Args:
        dataset_id: Dataset identifier
        full_parquet: Path to {dataset_id}_full.parquet
    
    Returns:
        (passed: bool, issues: List[str])
    """
    issues = []
    
    if not full_parquet.exists():
        return False, [f"File not found: {full_parquet}"]
    
    try:
        df = pd.read_parquet(full_parquet)
    except Exception as e:
        return False, [f"Failed to load parquet: {e}"]
    
    # Check that important columns are present
    important_cols = [
        "variant_id", "gene", "pos", "ref_aa", "alt_aa",
        "functional_score", "is_hit",
        "impact_score",
    ]
    
    for col in important_cols:
        if col not in df.columns:
            issues.append(f"Missing column: {col}")
    
    # Check for feature columns being all zero (likely error)
    feature_cols = ["impact_score", "conservation", "alphamissense_score", "spliceai_score"]
    for col in feature_cols:
        if col in df.columns:
            if (df[col] == 0).all():
                # This might be OK for some features, but is suspicious
                logger.warning(f"Feature column {col} is all zeros; might be a data issue")
    
    # Check row count preservation through pipeline
    # (should be <= original, not lose many rows)
    if len(df) == 0:
        issues.append("Dataset has 0 rows")
    
    # Check for duplicates
    dup_count = df.duplicated(subset=["variant_id"]).sum()
    if dup_count > 0:
        issues.append(f"Found {dup_count} duplicate variant IDs")
    
    # Check merge integrity: all variants should have gene, pos, etc.
    for col in ["gene", "pos", "variant_id"]:
        if col in df.columns and df[col].isna().any():
            null_count = df[col].isna().sum()
            issues.append(f"Column {col} has {null_count} nulls (should be no nulls)")
    
    passed = len(issues) == 0
    
    if passed:
        logger.info(f"✅ {dataset_id} full dataset checks passed")
    else:
        logger.error(f"❌ {dataset_id} full dataset checks found {len(issues)} issues:")
        for issue in issues:
            logger.error(f"   - {issue}")
    
    return passed, issues


def check_all_datasets(
    base_dir: Path = Path("data_processed/mave"),
    datasets_config_path: Path = Path("config/mave_datasets.yaml"),
) -> Dict[str, Dict[str, Tuple[bool, List[str]]]]:
    """
    Run all sanity checks on all datasets.
    
    Returns:
        Dict mapping dataset_id to dict of check results:
        {
            "ingestion": (passed: bool, issues: List[str]),
            "normalization": (passed: bool, issues: List[str]),
            "full": (passed: bool, issues: List[str]),
        }
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
    results = {}
    
    for dataset_cfg in datasets:
        dataset_id = dataset_cfg.get("dataset_id")
        if not dataset_id:
            continue
        
        dataset_dir = base_dir / dataset_id
        logger.info(f"\nChecking {dataset_id}...")
        
        results[dataset_id] = {
            "ingestion": check_ingested_dataset(
                dataset_id,
                dataset_dir / f"{dataset_id}_raw.parquet"
            ),
            "normalization": check_normalized_dataset(
                dataset_id,
                dataset_dir / f"{dataset_id}_normalized.parquet",
                dataset_cfg,
            ),
            "full": check_full_dataset(
                dataset_id,
                dataset_dir / f"{dataset_id}_full.parquet",
            ),
        }
    
    return results


def print_check_report(
    check_results: Dict[str, Dict[str, Tuple[bool, List[str]]]],
) -> None:
    """
    Print a formatted report of all sanity checks.
    
    Args:
        check_results: Results from check_all_datasets
    """
    total_checks = len(check_results) * 3
    passed_checks = sum(
        1 for dataset_results in check_results.values()
        for check_name, (passed, _) in dataset_results.items()
        if passed
    )
    
    logger.info(f"\n{'=' * 60}")
    logger.info(f"SANITY CHECK REPORT")
    logger.info(f"{'=' * 60}")
    logger.info(f"Overall: {passed_checks} / {total_checks} checks passed")
    
    for dataset_id in sorted(check_results.keys()):
        dataset_results = check_results[dataset_id]
        
        statuses = []
        for check_name in ["ingestion", "normalization", "full"]:
            passed, issues = dataset_results.get(check_name, (False, []))
            status = "✅" if passed else "❌"
            statuses.append(f"{check_name}:{status}")
        
        logger.info(f"\n{dataset_id}: {' | '.join(statuses)}")
        
        for check_name, (passed, issues) in dataset_results.items():
            if issues:
                logger.info(f"  {check_name}:")
                for issue in issues:
                    logger.info(f"    - {issue}")

