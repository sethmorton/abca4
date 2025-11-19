#!/usr/bin/env python3
"""
MAVE Normalization Module - Phase C: Standardize Scores and Define Hits

Normalizes raw functional scores to a standardized scale and defines what
constitutes a "hit" (e.g., loss-of-function variant) based on configurable
percentile thresholds.

Key Responsibilities:
  - Apply score normalization methods (z-score, quantile, min-max, identity)
  - Handle score direction (higher/lower = better)
  - Compute hit labels and thresholds
  - Validate normalization (sanity checks)
  - Output: {gene}_{dataset}_normalized.parquet with is_hit column

Extensibility:
  - Add new normalization methods in normalize_functional_score()
  - Customize hit definition in define_hits() by percent ile or absolute threshold
  - Easy to add domain-specific hit logic

Usage:
  from src.mave.pipeline.normalize import normalize_all_datasets
  results = normalize_all_datasets(Path("config/mave_datasets.yaml"))

Score normalization and hit definition for MAVE datasets.

This module:
1. Applies normalization (zscore, quantile, minmax)
2. Orients scores so higher = more functional
3. Defines hit labels based on percentile thresholds
4. Validates hit definitions
"""

import logging
from pathlib import Path
from typing import Dict, Any, Optional
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def normalize_functional_score(
    df: pd.DataFrame,
    method: str = "zscore",
    flip_sign: bool = False,
    clip_low: Optional[float] = None,
    clip_high: Optional[float] = None,
) -> pd.Series:
    """
    Normalize functional scores within a dataset.
    
    Args:
        df: DataFrame with 'functional_score_raw' column
        method: 'zscore', 'quantile', 'minmax', or 'identity'
        flip_sign: If True, multiply by -1 first
        clip_low: Lower percentile for clipping (before normalization)
        clip_high: Upper percentile for clipping (before normalization)
    
    Returns:
        Series of normalized scores
    """
    scores = df["functional_score_raw"].copy()
    
    # Flip sign if requested (inverts interpretation)
    if flip_sign:
        scores = -scores
        logger.debug("Flipped sign of scores")
    
    # Clip by percentile if requested
    if clip_low is not None:
        lower_bound = scores.quantile(clip_low)
        scores = scores.clip(lower=lower_bound)
        logger.debug(f"Clipped lower at percentile {clip_low}: {lower_bound:.4f}")
    
    if clip_high is not None:
        upper_bound = scores.quantile(clip_high)
        scores = scores.clip(upper=upper_bound)
        logger.debug(f"Clipped upper at percentile {clip_high}: {upper_bound:.4f}")
    
    # Apply normalization method
    if method == "identity":
        normalized = scores
    
    elif method == "zscore":
        mean_val = scores.mean()
        std_val = scores.std()
        if std_val == 0:
            logger.warning("Standard deviation is 0; using identity normalization")
            normalized = scores
        else:
            normalized = (scores - mean_val) / std_val
    
    elif method == "minmax":
        min_val = scores.min()
        max_val = scores.max()
        if max_val == min_val:
            logger.warning("Min and max are equal; using identity normalization")
            normalized = scores
        else:
            normalized = (scores - min_val) / (max_val - min_val)
    
    elif method == "quantile":
        # Map to [0, 1] by rank
        normalized = scores.rank(method="average", na_option="keep") / len(scores)
    
    else:
        logger.warning(f"Unknown normalization method: {method}; using identity")
        normalized = scores
    
    return normalized


def define_hits(
    df: pd.DataFrame,
    direction: str = "low_is_loss_of_function",
    percentile_threshold: float = 0.2,
) -> tuple[pd.Series, pd.Series]:
    """
    Define hit labels based on percentile threshold.
    
    Args:
        df: DataFrame with 'functional_score' column
        direction: 'low_is_loss_of_function' or 'high_is_gain_of_function'
        percentile_threshold: Percentile threshold (e.g., 0.2 for bottom 20%)
    
    Returns:
        (is_hit Series, hit_label Series)
    """
    scores = df["functional_score"]
    
    if direction == "low_is_loss_of_function":
        # Hits are in the lower tail
        threshold = scores.quantile(percentile_threshold)
        is_hit = scores <= threshold
        hit_label = is_hit.map({True: "lof", False: "non_lof"})
        logger.debug(
            f"Hit definition: low_is_loss_of_function, "
            f"threshold={threshold:.4f}, hits={is_hit.sum()}"
        )
    
    elif direction == "high_is_gain_of_function":
        # Hits are in the upper tail
        threshold = scores.quantile(1.0 - percentile_threshold)
        is_hit = scores >= threshold
        hit_label = is_hit.map({True: "gof", False: "non_gof"})
        logger.debug(
            f"Hit definition: high_is_gain_of_function, "
            f"threshold={threshold:.4f}, hits={is_hit.sum()}"
        )
    
    else:
        logger.error(f"Unknown direction: {direction}")
        is_hit = pd.Series(False, index=df.index)
        hit_label = pd.Series("unknown", index=df.index)
    
    return is_hit, hit_label


def normalize_dataset(
    dataset_id: str,
    raw_parquet: Path,
    config: Dict[str, Any],
    output_dir: Path,
) -> Optional[pd.DataFrame]:
    """
    Load raw ingested data, apply normalization, and define hits.
    
    Args:
        dataset_id: Dataset identifier
        raw_parquet: Path to dms_raw.parquet
        config: Dataset config from mave_datasets.yaml
        output_dir: Directory for output
    
    Returns:
        DataFrame with normalized scores and hit labels, or None if failed
    """
    if not raw_parquet.exists():
        logger.error(f"Raw parquet not found: {raw_parquet}")
        return None
    
    try:
        df = pd.read_parquet(raw_parquet)
        logger.info(f"Loaded {len(df)} rows from {raw_parquet.name}")
    except Exception as e:
        logger.error(f"Failed to read {raw_parquet}: {e}")
        return None
    
    # Extract normalization and hit config
    normalization_cfg = config.get("normalization", {})
    hit_cfg = config.get("hit_definition", {})
    
    norm_method = normalization_cfg.get("method", "identity")
    flip_sign = normalization_cfg.get("flip_sign", False)
    clip_low = normalization_cfg.get("clip_low")
    clip_high = normalization_cfg.get("clip_high")
    
    hit_direction = hit_cfg.get("direction", "low_is_loss_of_function")
    percentile_threshold = hit_cfg.get("percentile_threshold", 0.2)
    
    # Apply normalization
    logger.info(f"Normalizing with method={norm_method}, flip_sign={flip_sign}")
    df["functional_score"] = normalize_functional_score(
        df,
        method=norm_method,
        flip_sign=flip_sign,
        clip_low=clip_low,
        clip_high=clip_high,
    )
    
    # Define hits
    logger.info(
        f"Defining hits: direction={hit_direction}, "
        f"percentile_threshold={percentile_threshold}"
    )
    is_hit, hit_label = define_hits(
        df,
        direction=hit_direction,
        percentile_threshold=percentile_threshold,
    )
    df["is_hit"] = is_hit
    df["hit_label"] = hit_label
    
    # Stats
    n_total = len(df)
    n_hits = is_hit.sum()
    hit_rate = n_hits / n_total if n_total > 0 else 0.0
    
    logger.info(f"Dataset stats:")
    logger.info(f"  Total variants: {n_total}")
    logger.info(f"  Hits: {n_hits} ({hit_rate:.1%})")
    logger.info(f"  Score range: [{df['functional_score'].min():.4f}, {df['functional_score'].max():.4f}]")
    logger.info(f"  Score mean: {df['functional_score'].mean():.4f}")
    logger.info(f"  Score std: {df['functional_score'].std():.4f}")
    
    # Sanity check
    expected_hit_rate = percentile_threshold
    actual_hit_rate = hit_rate
    deviation = abs(actual_hit_rate - expected_hit_rate)
    
    if deviation > 0.05:  # 5% tolerance
        logger.warning(
            f"Hit rate deviation: expected {expected_hit_rate:.1%}, got {actual_hit_rate:.1%}. "
            f"This may indicate ties or truncation; please review."
        )
    
    # Write output
    output_dir.mkdir(parents=True, exist_ok=True)
    output_parquet = output_dir / f"{dataset_id}_normalized.parquet"
    
    try:
        df.to_parquet(output_parquet, index=False)
        logger.info(f"Wrote {len(df)} rows to {output_parquet}")
    except Exception as e:
        logger.error(f"Failed to write {output_parquet}: {e}")
        return None
    
    return df


def normalize_all_datasets(
    datasets_config_path: Path,
    base_dir: Path = Path("data_processed/mave"),
) -> Dict[str, pd.DataFrame]:
    """
    Normalize all ingested datasets.
    
    Args:
        datasets_config_path: Path to mave_datasets.yaml
        base_dir: Base directory for input/output
    
    Returns:
        Dict mapping dataset_id to normalized DataFrame
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
        
        dataset_dir = base_dir / dataset_id
        raw_parquet = dataset_dir / f"{dataset_id}_raw.parquet"
        
        logger.info(f"\nNormalizing {dataset_id}...")
        df = normalize_dataset(dataset_id, raw_parquet, dataset_cfg, dataset_dir)
        
        if df is not None:
            results[dataset_id] = df
    
    return results


