#!/usr/bin/env python3
"""
Evaluation metrics for MAVE variant selection.

Computes:
- Hit recall and precision
- Domain/cluster coverage
- Mean functional score
"""

import logging
from pathlib import Path
from typing import Dict, Any
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def compute_metrics(
    df_all: pd.DataFrame,
    df_selected: pd.DataFrame,
) -> Dict[str, float]:
    """
    Compute evaluation metrics for a selection strategy.
    
    Args:
        df_all: All variants with 'is_hit' column
        df_selected: Selected variants
    
    Returns:
        Dict with metrics:
        - hit_recall: Fraction of total hits in selection
        - hit_precision: Fraction of selected that are hits
        - n_hits_total: Total number of hits
        - n_hits_selected: Number of hits in selection
        - n_selected: Number selected
        - mean_functional_score: Mean score of selected
    """
    n_selected = len(df_selected)
    n_hits_total = df_all["is_hit"].sum() if "is_hit" in df_all.columns else 0
    
    # Ensure we're comparing the same variant IDs
    selected_ids = set(df_selected.get("variant_id", df_selected.index))
    all_ids = set(df_all.get("variant_id", df_all.index))
    
    # Find which selected variants are actually hits
    if "is_hit" in df_selected.columns:
        n_hits_selected = df_selected["is_hit"].sum()
    else:
        # Re-lookup hit status from df_all
        n_hits_selected = 0
        for selected_id in selected_ids:
            hit_rows = df_all[df_all.get("variant_id") == selected_id]
            if not hit_rows.empty and hit_rows.iloc[0].get("is_hit"):
                n_hits_selected += 1
    
    # Compute metrics
    hit_recall = n_hits_selected / n_hits_total if n_hits_total > 0 else 0.0
    hit_precision = n_hits_selected / n_selected if n_selected > 0 else 0.0
    
    mean_functional_score = 0.0
    if "functional_score" in df_selected.columns:
        mean_functional_score = df_selected["functional_score"].mean()
    
    # Cluster/domain coverage
    cluster_coverage = {}
    if "cluster_id" in df_selected.columns:
        for cluster_id in df_selected["cluster_id"].unique():
            if pd.isna(cluster_id):
                continue
            n_in_cluster_total = (df_all["cluster_id"] == cluster_id).sum()
            n_in_cluster_selected = (df_selected["cluster_id"] == cluster_id).sum()
            if n_in_cluster_total > 0:
                coverage = n_in_cluster_selected / n_in_cluster_total
                cluster_coverage[f"cluster_{int(cluster_id)}_coverage"] = coverage
    
    metrics = {
        "hit_recall": hit_recall,
        "hit_precision": hit_precision,
        "n_hits_total": n_hits_total,
        "n_hits_selected": n_hits_selected,
        "n_selected": n_selected,
        "mean_functional_score": mean_functional_score,
    }
    
    metrics.update(cluster_coverage)
    
    return metrics


def evaluate_dataset(
    dataset_id: str,
    df_full: pd.DataFrame,
    strategy_results: Dict[str, pd.DataFrame],
    k_values: list = None,
) -> pd.DataFrame:
    """
    Evaluate all strategies on a dataset.
    
    Args:
        dataset_id: Dataset identifier
        df_full: Full dataset with variants and functional scores
        strategy_results: Dict mapping strategy name to selected DataFrame
        k_values: K values used (for logging)
    
    Returns:
        DataFrame with rows: one per strategy, columns for each metric
    """
    if k_values is None:
        k_values = [len(df_selected) for df_selected in strategy_results.values()]
    
    rows = []
    for strategy_name, df_selected in strategy_results.items():
        metrics = compute_metrics(df_full, df_selected)
        
        # Add identifiers
        metrics["dataset_id"] = dataset_id
        metrics["strategy"] = strategy_name
        metrics["k"] = len(df_selected)
        
        rows.append(metrics)
        
        # Log
        logger.info(
            f"  {strategy_name}: "
            f"recall={metrics['hit_recall']:.3f}, "
            f"precision={metrics['hit_precision']:.3f}, "
            f"mean_score={metrics['mean_functional_score']:.3f}"
        )
    
    df_metrics = pd.DataFrame(rows)
    
    return df_metrics


def evaluate_all_datasets(
    base_dir: Path = Path("data_processed/mave"),
    k_values: list = None,
) -> pd.DataFrame:
    """
    Evaluate all datasets and strategies.
    
    This is a placeholder that assumes:
    1. Each dataset has {dataset_id}_full.parquet
    2. Selection results are stored (TODO: define format)
    
    Args:
        base_dir: Base directory with datasets
        k_values: K values to evaluate (e.g., [30, 50])
    
    Returns:
        Aggregated metrics DataFrame
    """
    if k_values is None:
        k_values = [30, 50]
    
    logger.info(f"Evaluating all datasets with k={k_values}")
    
    all_metrics = []
    
    # Find all dataset directories
    dataset_dirs = [d for d in base_dir.iterdir() if d.is_dir()]
    
    for dataset_dir in sorted(dataset_dirs):
        dataset_id = dataset_dir.name
        full_parquet = dataset_dir / f"{dataset_id}_full.parquet"
        
        if not full_parquet.exists():
            logger.warning(f"Skipping {dataset_id}: no full parquet")
            continue
        
        try:
            df_full = pd.read_parquet(full_parquet)
            logger.info(f"Loaded {len(df_full)} variants from {dataset_id}")
        except Exception as e:
            logger.error(f"Failed to load {full_parquet}: {e}")
            continue
        
        # TODO: Load selection results for this dataset
        # For now, this is a template for the junior engineer to fill in
        
        # # Placeholder: assume strategy_results loaded from saved files
        # strategy_results = {
        #     "random": df_selected_random,
        #     "top_model_score": df_selected_model,
        #     ...
        # }
        # 
        # df_metrics = evaluate_dataset(dataset_id, df_full, strategy_results)
        # all_metrics.append(df_metrics)
    
    if not all_metrics:
        logger.warning("No evaluation results collected")
        return pd.DataFrame()
    
    df_all_metrics = pd.concat(all_metrics, ignore_index=True)
    
    logger.info(f"Evaluation complete: {len(df_all_metrics)} rows")
    
    return df_all_metrics


def save_metrics(
    metrics_df: pd.DataFrame,
    dataset_id: str,
    output_dir: Path,
) -> Path:
    """
    Save metrics to CSV.
    
    Args:
        metrics_df: Metrics DataFrame
        dataset_id: Dataset identifier
        output_dir: Output directory
    
    Returns:
        Path to saved file
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"mave_{dataset_id}_metrics.csv"
    
    try:
        metrics_df.to_csv(output_path, index=False)
        logger.info(f"Saved metrics to {output_path}")
    except Exception as e:
        logger.error(f"Failed to save metrics: {e}")
        return None
    
    return output_path


def sanity_check_metrics(
    metrics_df: pd.DataFrame,
    dataset_id: str,
) -> list:
    """
    Perform sanity checks on computed metrics.
    
    Returns list of warnings/issues found.
    
    Args:
        metrics_df: Metrics DataFrame
        dataset_id: Dataset identifier
    
    Returns:
        List of issue strings
    """
    issues = []
    
    # Check 1: Oracle recall should be high
    oracle_row = metrics_df[metrics_df["strategy"] == "oracle_functional"]
    if not oracle_row.empty:
        oracle_recall = oracle_row.iloc[0]["hit_recall"]
        if oracle_recall < 0.6:
            issues.append(
                f"Oracle recall is low ({oracle_recall:.3f}); "
                f"possible hit misdefinition?"
            )
    
    # Check 2: Random recall should roughly match hit fraction
    random_row = metrics_df[metrics_df["strategy"] == "random"]
    if not random_row.empty and not oracle_row.empty:
        n_hits_total = oracle_row.iloc[0]["n_hits_total"]
        n_selected = random_row.iloc[0]["n_selected"]
        
        if n_hits_total > 0:
            expected_random_hits = (n_selected / len(metrics_df)) * n_hits_total  # Rough estimate
            actual_random_hits = random_row.iloc[0]["n_hits_selected"]
            
            if actual_random_hits == 0 and n_hits_total > 0:
                issues.append(
                    f"Random strategy found {actual_random_hits} hits "
                    f"but total hits = {n_hits_total}; something is wrong"
                )
    
    # Check 3: Strand should not be worse than random
    strand_row = metrics_df[metrics_df["strategy"] == "strand"]
    if not strand_row.empty and not random_row.empty:
        strand_recall = strand_row.iloc[0]["hit_recall"]
        random_recall = random_row.iloc[0]["hit_recall"]
        
        if strand_recall < random_recall * 0.9:  # 10% tolerance
            issues.append(
                f"Strand recall ({strand_recall:.3f}) is much worse than random "
                f"({random_recall:.3f}); check selection code"
            )
    
    return issues


