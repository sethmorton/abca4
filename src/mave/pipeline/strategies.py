#!/usr/bin/env python3
"""
Selection strategies for MAVE variant evaluation.

Implements:
- Baseline strategies (random, top-model, top-conservation, oracle)
- Strand selection (reuses existing optimizer)
- Shared infrastructure for comparison
"""

import logging
from typing import Dict, List, Optional
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def select_random(
    df: pd.DataFrame,
    k: int,
    random_state: int = 42,
) -> pd.DataFrame:
    """
    Random selection baseline.
    
    Args:
        df: DataFrame with variants
        k: Number of variants to select
        random_state: Random seed for reproducibility
    
    Returns:
        Selected variants with 'selected' and 'rank' columns added
    """
    if k > len(df):
        logger.warning(f"k={k} > n={len(df)}; selecting all")
        k = len(df)
    
    np.random.seed(random_state)
    indices = np.random.choice(len(df), size=k, replace=False)
    
    df_selected = df.iloc[indices].copy()
    df_selected["selected"] = True
    df_selected["rank"] = range(1, len(df_selected) + 1)
    df_selected["strategy"] = "random"
    
    logger.info(f"Random: selected {len(df_selected)} / {len(df)} variants")
    
    return df_selected


def select_top_model_score(
    df: pd.DataFrame,
    k: int,
    score_col: str = "impact_score",
) -> pd.DataFrame:
    """
    Select top variants by model score.
    
    Args:
        df: DataFrame with variants
        k: Number of variants to select
        score_col: Column name for model score
    
    Returns:
        Top-K variants sorted by score
    """
    if score_col not in df.columns:
        # Fallback to other common names
        for alt_col in ["alphamissense_score", "model_score"]:
            if alt_col in df.columns:
                score_col = alt_col
                break
        else:
            raise ValueError(f"Score column {score_col} not found")
    
    if k > len(df):
        logger.warning(f"k={k} > n={len(df)}; selecting all")
        k = len(df)
    
    df_sorted = df.sort_values(score_col, ascending=False)
    df_selected = df_sorted.head(k).copy()
    df_selected["selected"] = True
    df_selected["rank"] = range(1, len(df_selected) + 1)
    df_selected["strategy"] = "top_model_score"
    
    logger.info(f"Top model score: selected {len(df_selected)} / {len(df)} variants")
    
    return df_selected


def select_top_conservation(
    df: pd.DataFrame,
    k: int,
    conservation_col: str = "conservation",
) -> pd.DataFrame:
    """
    Select top variants by conservation score.
    
    Args:
        df: DataFrame with variants
        k: Number of variants to select
        conservation_col: Column name for conservation
    
    Returns:
        Top-K variants by conservation
    """
    if conservation_col not in df.columns:
        for alt_col in ["phyloP100way", "phastCons100way"]:
            if alt_col in df.columns:
                conservation_col = alt_col
                break
        else:
            logger.warning(f"Conservation column not found, using model score fallback")
            return select_top_model_score(df, k)
    
    if k > len(df):
        logger.warning(f"k={k} > n={len(df)}; selecting all")
        k = len(df)
    
    df_sorted = df.sort_values(conservation_col, ascending=False, na_position="last")
    df_selected = df_sorted.head(k).copy()
    df_selected["selected"] = True
    df_selected["rank"] = range(1, len(df_selected) + 1)
    df_selected["strategy"] = "top_conservation"
    
    logger.info(f"Top conservation: selected {len(df_selected)} / {len(df)} variants")
    
    return df_selected


def select_top_functional_score(
    df: pd.DataFrame,
    k: int,
    score_col: str = "functional_score",
) -> pd.DataFrame:
    """
    Oracle: Select variants with LOWEST functional scores (loss-of-function hits).
    
    This is an oracle baseline that uses the ground truth functional scores.
    For loss-of-function assays (low_is_loss_of_function), lower scores = hits.
    It provides an upper bound on achievable performance.
    
    Args:
        df: DataFrame with variants
        k: Number of variants to select
        score_col: Column with true functional score
    
    Returns:
        Top-K variants with lowest functional scores (most loss-of-function)
    """
    if score_col not in df.columns:
        raise ValueError(f"Oracle score column {score_col} not found")
    
    if k > len(df):
        logger.warning(f"k={k} > n={len(df)}; selecting all")
        k = len(df)
    
    # IMPORTANT: For loss-of-function hits, LOWER scores = hits
    # So we sort ASCENDING and take the LOWEST k variants
    df_sorted = df.sort_values(score_col, ascending=True, na_position="last")
    df_selected = df_sorted.head(k).copy()
    df_selected["selected"] = True
    df_selected["rank"] = range(1, len(df_selected) + 1)
    df_selected["strategy"] = "oracle_functional"
    
    logger.info(f"Oracle (functional): selected {len(df_selected)} / {len(df)} variants (lowest scores)")
    
    return df_selected


def select_strand(
    df: pd.DataFrame,
    k: int,
    config: Optional[Dict] = None,
    select_greedy_func=None,
    lambda_penalty: float = 0.6,
) -> pd.DataFrame:
    """
    Strand-style selection with coverage constraints (your brain algorithm).
    
    This is the Strand selection algorithm that combines:
    - Model scores (impact_score from ABCA4)
    - Coverage penalties (cluster/domain targets)
    - Greedy optimization
    
    Args:
        df: DataFrame with variants (must have impact_score)
        k: Number of variants to select
        config: Gene/selection config (optional)
        select_greedy_func: Function from src.reward.optimization.VariantOptimizer.select_greedy
        lambda_penalty: Coverage penalty weight (0.6 default)
    
    Returns:
        Selected variants DataFrame
    """
    if k > len(df):
        logger.warning(f"k={k} > n={len(df)}; selecting all")
        k = len(df)
    
    # Use VariantOptimizer.select_greedy if provided, otherwise fallback to model scores
    if select_greedy_func is None:
        logger.warning("No select_greedy_func; using top model scores instead")
        return select_top_model_score(df, k, score_col="impact_score")
    
    # Extract config parameters
    if config:
        lambda_penalty = config.get("selection", {}).get("lambda_penalty", lambda_penalty)
    
    # Prepare cluster targets if available
    cluster_targets = {}
    if "cluster_id" in df.columns and "cluster_target" in df.columns:
        cluster_targets = (
            df[["cluster_id", "cluster_target"]]
            .dropna()
            .drop_duplicates("cluster_id")
            .set_index("cluster_id")["cluster_target"]
            .to_dict()
        )
    
    logger.info(f"Strand: lambda_penalty={lambda_penalty}, clusters={len(cluster_targets)}")
    
    # Call the Strand selector (greedy with coverage constraints)
    df_selected = select_greedy_func(
        df,
        k=k,
        lambda_penalty=lambda_penalty,
        cluster_targets=cluster_targets if cluster_targets else None,
    )
    
    df_selected["strategy"] = "strand"
    
    logger.info(f"Strand: selected {len(df_selected)} / {len(df)} variants")
    
    return df_selected


def run_all_strategies(
    df: pd.DataFrame,
    k: int,
    config: Dict,
    select_greedy_func=None,
) -> Dict[str, pd.DataFrame]:
    """
    Run all selection strategies for comparison.
    
    Args:
        df: DataFrame with variants
        k: Number to select
        config: Gene configuration
        select_greedy_func: Strand selector function
    
    Returns:
        Dict mapping strategy name to selected variants
    """
    results = {}
    
    strategies = [
        ("random", lambda: select_random(df, k)),
        ("top_model_score", lambda: select_top_model_score(df, k)),
        ("top_conservation", lambda: select_top_conservation(df, k)),
        ("oracle_functional", lambda: select_top_functional_score(df, k)),
    ]
    
    # Add Strand if function provided
    if select_greedy_func is not None:
        strategies.append(
            ("strand", lambda: select_strand(df, k, config, select_greedy_func))
        )
    
    logger.info(f"Running {len(strategies)} selection strategies for k={k}")
    
    for strategy_name, strategy_func in strategies:
        try:
            df_selected = strategy_func()
            results[strategy_name] = df_selected
            logger.info(f"  {strategy_name}: OK")
        except Exception as e:
            logger.error(f"  {strategy_name}: FAILED - {e}")
    
    return results


