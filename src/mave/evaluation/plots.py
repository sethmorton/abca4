#!/usr/bin/env python3
"""
Visualization helpers for MAVE evaluation results.

Generates plots for inclusion in HTML reports or standalone viewing.
"""

import logging
from pathlib import Path
from typing import Dict, Optional, List
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def ensure_plotly():
    """Ensure plotly is available; log error if not."""
    try:
        import plotly.graph_objects as go
        import plotly.express as px
        return go, px
    except ImportError:
        logger.error("plotly not installed; plotting disabled")
        return None, None


def plot_hit_recall_by_strategy(
    metrics_df: pd.DataFrame,
    dataset_id: str,
    k: int,
    output_dir: Optional[Path] = None,
) -> Optional[str]:
    """
    Bar chart of hit recall by strategy.
    
    Args:
        metrics_df: Metrics from evaluate_dataset
        dataset_id: Dataset identifier
        k: K value for filtering
        output_dir: If provided, save HTML to this directory
    
    Returns:
        HTML string for embedding, or None if plotting failed
    """
    go, px = ensure_plotly()
    if go is None:
        return None
    
    # Filter for this K
    df = metrics_df[metrics_df["k"] == k].copy() if "k" in metrics_df.columns else metrics_df.copy()
    
    if df.empty:
        logger.warning(f"No metrics for k={k}")
        return None
    
    # Create bar chart
    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=df["strategy"],
        y=df["hit_recall"],
        marker_color="rgba(0, 100, 200, 0.8)",
        text=df["hit_recall"].apply(lambda x: f"{x:.2%}"),
        textposition="auto",
    ))
    
    fig.update_layout(
        title=f"Hit Recall by Strategy - {dataset_id} (K={k})",
        xaxis_title="Strategy",
        yaxis_title="Hit Recall",
        yaxis=dict(range=[0, 1]),
        height=400,
        showlegend=False,
    )
    
    # Save if output_dir provided
    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{dataset_id}_recall_k{k}.html"
        try:
            fig.write_html(str(output_path))
            logger.info(f"Saved plot to {output_path}")
        except Exception as e:
            logger.error(f"Failed to save plot: {e}")
    
    # Return HTML div for embedding
    return fig.to_html(include_plotlyjs="cdn", div_id=f"recall_{dataset_id}_{k}")


def plot_precision_recall_scatter(
    metrics_df: pd.DataFrame,
    dataset_id: str,
    output_dir: Optional[Path] = None,
) -> Optional[str]:
    """
    Scatter plot: precision vs recall for all strategies.
    
    Args:
        metrics_df: Metrics DataFrame
        dataset_id: Dataset identifier
        output_dir: If provided, save HTML
    
    Returns:
        HTML string for embedding, or None
    """
    go, px = ensure_plotly()
    if go is None:
        return None
    
    df = metrics_df.copy()
    if df.empty:
        return None
    
    fig = go.Figure()
    
    # Color by strategy
    strategies = df["strategy"].unique()
    colors = px.colors.qualitative.Plotly
    color_map = {s: colors[i % len(colors)] for i, s in enumerate(strategies)}
    
    for strategy in sorted(strategies):
        df_strat = df[df["strategy"] == strategy]
        fig.add_trace(go.Scatter(
            x=df_strat["hit_recall"],
            y=df_strat["hit_precision"],
            mode="markers",
            name=strategy,
            marker=dict(size=10, color=color_map[strategy]),
        ))
    
    fig.update_layout(
        title=f"Precision vs Recall - {dataset_id}",
        xaxis_title="Hit Recall",
        yaxis_title="Hit Precision",
        xaxis=dict(range=[0, 1]),
        yaxis=dict(range=[0, 1]),
        height=500,
        hovermode="closest",
    )
    
    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{dataset_id}_precision_recall.html"
        try:
            fig.write_html(str(output_path))
            logger.info(f"Saved plot to {output_path}")
        except Exception as e:
            logger.error(f"Failed to save plot: {e}")
    
    return fig.to_html(include_plotlyjs="cdn", div_id=f"pr_{dataset_id}")


def plot_functional_score_distribution(
    df_full: pd.DataFrame,
    dataset_id: str,
    output_dir: Optional[Path] = None,
) -> Optional[str]:
    """
    Histogram of functional scores with hit region highlighted.
    
    Args:
        df_full: Full dataset with functional_score and is_hit columns
        dataset_id: Dataset identifier
        output_dir: If provided, save HTML
    
    Returns:
        HTML string, or None
    """
    go, px = ensure_plotly()
    if go is None:
        return None
    
    if "functional_score" not in df_full.columns:
        logger.warning("No functional_score column")
        return None
    
    df = df_full.copy()
    
    fig = go.Figure()
    
    # Separate hits and non-hits
    if "is_hit" in df.columns:
        df_hits = df[df["is_hit"] == True]
        df_non_hits = df[df["is_hit"] == False]
    else:
        df_hits = pd.DataFrame()
        df_non_hits = df
    
    fig.add_trace(go.Histogram(
        x=df_non_hits["functional_score"],
        name="Non-hits",
        marker_color="rgba(100, 100, 100, 0.7)",
        opacity=0.75,
    ))
    
    if not df_hits.empty:
        fig.add_trace(go.Histogram(
            x=df_hits["functional_score"],
            name="Hits",
            marker_color="rgba(255, 0, 0, 0.8)",
            opacity=0.75,
        ))
    
    fig.update_layout(
        title=f"Functional Score Distribution - {dataset_id}",
        xaxis_title="Functional Score",
        yaxis_title="Count",
        barmode="overlay",
        height=400,
    )
    
    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{dataset_id}_score_dist.html"
        try:
            fig.write_html(str(output_path))
            logger.info(f"Saved plot to {output_path}")
        except Exception as e:
            logger.error(f"Failed to save plot: {e}")
    
    return fig.to_html(include_plotlyjs="cdn", div_id=f"dist_{dataset_id}")


def plot_model_vs_functional_scores(
    df_full: pd.DataFrame,
    dataset_id: str,
    df_selected: Optional[pd.DataFrame] = None,
    output_dir: Optional[Path] = None,
) -> Optional[str]:
    """
    Scatter: model score vs functional score, highlighting selected variants.
    
    Args:
        df_full: Full dataset
        dataset_id: Dataset identifier
        df_selected: Optional selected variants (highlighted in red)
        output_dir: If provided, save HTML
    
    Returns:
        HTML string, or None
    """
    go, px = ensure_plotly()
    if go is None:
        return None
    
    if "impact_score" not in df_full.columns or "functional_score" not in df_full.columns:
        logger.warning("Missing impact_score or functional_score")
        return None
    
    df = df_full.copy()
    
    fig = go.Figure()
    
    # Background: all variants
    fig.add_trace(go.Scatter(
        x=df["impact_score"],
        y=df["functional_score"],
        mode="markers",
        name="All variants",
        marker=dict(size=5, color="rgba(100, 100, 100, 0.3)"),
        text=df.get("variant_id", ""),
        hovertemplate="<b>%{text}</b><br>Model: %{x:.3f}<br>Functional: %{y:.3f}<extra></extra>",
    ))
    
    # Overlay: selected variants
    if df_selected is not None and not df_selected.empty:
        fig.add_trace(go.Scatter(
            x=df_selected["impact_score"],
            y=df_selected["functional_score"],
            mode="markers",
            name="Selected",
            marker=dict(size=8, color="red", symbol="star"),
            text=df_selected.get("variant_id", ""),
            hovertemplate="<b>%{text}</b><br>Model: %{x:.3f}<br>Functional: %{y:.3f}<extra></extra>",
        ))
    
    fig.update_layout(
        title=f"Model Score vs Functional Score - {dataset_id}",
        xaxis_title="Model Score (Impact)",
        yaxis_title="Functional Score",
        height=500,
        hovermode="closest",
    )
    
    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{dataset_id}_model_vs_functional.html"
        try:
            fig.write_html(str(output_path))
            logger.info(f"Saved plot to {output_path}")
        except Exception as e:
            logger.error(f"Failed to save plot: {e}")
    
    return fig.to_html(include_plotlyjs="cdn", div_id=f"scatter_{dataset_id}")


def create_summary_plots(
    base_dir: Path = Path("data_processed/mave"),
    output_dir: Path = Path("fig/mave"),
) -> Dict[str, str]:
    """
    Create all plots for all datasets.
    
    Returns dict mapping dataset_id to dict of plots.
    """
    logger.info(f"Generating plots for all datasets...")
    
    all_plots = {}
    
    # TODO: Implement
    # This is a template for the junior engineer to:
    # 1. Iterate over all dataset directories
    # 2. Load full parquets and metrics
    # 3. Call plot functions
    # 4. Aggregate and return results
    
    return all_plots


