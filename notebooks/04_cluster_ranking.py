#!/usr/bin/env python3
"""
ABCA4 Campaign – Cluster Ranking & Biological Scoring

Analyze variant clusters and rank them by biological relevance score.

Run interactively:  marimo edit notebooks/04_cluster_ranking.py
Run as dashboard:   marimo run notebooks/04_cluster_ranking.py
Run as script:      python notebooks/04_cluster_ranking.py
"""

import marimo

__generated_with = "0.17.8"
app = marimo.App()


@app.cell
def _():
    """Import core libraries."""
    import marimo as mo
    import pandas as pd
    import numpy as np
    from pathlib import Path
    return Path, mo, np, pd


@app.cell
def _(Path):
    """Define data paths."""
    NOTEBOOKS_DIR = Path(__file__).resolve().parent
    CAMPAIGN_ROOT = NOTEBOOKS_DIR.parent
    DATASETS_DIR = CAMPAIGN_ROOT / ".datasets"

    parquet_path = DATASETS_DIR / "variants_scored1.parquet"

    return (parquet_path,)


@app.cell
def _(mo):
    """
    ## Cluster Ranking & Biological Scoring

    Load variant data and rank clusters by biological relevance.
    """
    mo.md(__doc__)
    return


@app.cell
def _(parquet_path, pd):
    """Load variant data."""
    df = pd.read_parquet(parquet_path)
    return (df,)


@app.cell
def _(df, mo):
    """Display data overview."""
    mo.md(f"""
    ### Data Overview

    - **Total variants:** {len(df):,}
    - **Total clusters:** {df['cluster'].nunique()}
    - **Columns:** {len(df.columns)}
    """)

    mo.ui.table(df.head(10))
    return


@app.cell
def _(df):
    """Get unique clusters."""
    unique_clusters = df['cluster'].unique()
    return (unique_clusters,)


@app.cell
def _(df, mo, pd, unique_clusters):
    """Display all available clusters."""
    mo.md("### Available Clusters")

    _clusters_df = pd.DataFrame({
        "Cluster": unique_clusters,
        "Count": [len(df[df['cluster'] == c]) for c in unique_clusters]
    }).sort_values("Count", ascending=False)

    mo.ui.table(_clusters_df)
    return


@app.cell
def _(df):
    """Identify score columns."""
    # Find all columns that might contain scores
    score_columns = [col for col in df.columns if 'score' in col.lower() or 'pred' in col.lower() or 'prob' in col.lower()]

    # Filter to actual numeric score columns
    all_score_cols = [
        'alphamissense_score',
        'missense_combined_score', 
        'spliceai_max_score',
        'phylop_score',
        'phastcons_score',
        'conservation_score',
        'tss_window_score',
        'regulatory_score',
        'model_score'
    ]

    # Only include columns that actually exist in the dataframe
    available_score_cols = [col for col in all_score_cols if col in df.columns]

    return (available_score_cols,)


@app.cell
def _(available_score_cols, df, mo, pd):
    """Display score column information."""
    mo.md("### Available Score Columns")

    _score_info = []
    for col in available_score_cols[:5]:  # Show first 5
        _sample_vals = df[col].dropna().iloc[:5].tolist()
        _score_info.append({
            "Column": col,
            "Sample Values": str(_sample_vals)[:100]  # Truncate for display
        })

    if _score_info:
        mo.ui.table(pd.DataFrame(_score_info))
    return


@app.cell
def _(mo):
    """Score column selection."""
    mo.md("### Select Score Columns for Ranking")

    score_column_options = mo.ui.multiselect(
        options=[
            'alphamissense_score',
            'missense_combined_score', 
            'spliceai_max_score',
            'phylop_score',
            'phastcons_score',
            'conservation_score',
            'tss_window_score',
            'regulatory_score',
            'model_score'
        ],
        value=[
            'alphamissense_score',
            'missense_combined_score', 
            'spliceai_max_score',
            'phylop_score',
            'phastcons_score',
            'conservation_score',
            'tss_window_score',
            'regulatory_score',
            'model_score'
        ],
        label="Score Columns"
    )

    return (score_column_options,)


@app.cell
def _():
    """Define scoring function."""
    def calculate_cluster_score(cluster_df, score_cols):
        """
        Calculate a biologically relevant score for a cluster.
        Score = sum of impact scores - penalty for uncovered important domains
        """
        # Sum of all available impact scores
        impact_sum = 0
        for score_col in score_cols:
            if score_col in cluster_df.columns:
                # Handle NaN values
                impact_sum += cluster_df[score_col].fillna(0).sum()

        # Define important domains that should be covered
        important_domains = ['TMD', 'NBD1', 'NBD2', 'CTD']

        # Check if this cluster covers any important domain
        cluster_name = cluster_df['cluster'].iloc[0]
        covers_important_domain = any(domain in cluster_name for domain in important_domains)

        # Penalty if no important domain is covered
        penalty = 0 if covers_important_domain else 50

        return impact_sum - penalty

    return (calculate_cluster_score,)


@app.cell
def _(calculate_cluster_score, df, score_column_options):
    """Calculate scores for all clusters."""
    cluster_scores = {}

    for cluster_name in df['cluster'].unique():
        cluster_df = df[df['cluster'] == cluster_name]
        c_score = calculate_cluster_score(cluster_df, score_column_options.value)
        cluster_scores[cluster_name] = c_score

    # Rank clusters by score
    ranked_clusters = sorted(cluster_scores.items(), key=lambda x: x[1], reverse=True)

    return (ranked_clusters,)


@app.cell
def _(df, mo, pd, ranked_clusters):
    """Display ranked clusters."""
    mo.md("### Clusters Ranked by Biological Relevance Score")

    _ranking_data = []
    for i, (cluster, score) in enumerate(ranked_clusters, 1):
        _count = len(df[df['cluster'] == cluster])
        _ranking_data.append({
            "Rank": i,
            "Cluster": cluster,
            "Score": f"{score:.2f}",
            "Variant Count": _count
        })

    _ranking_df = pd.DataFrame(_ranking_data)
    mo.ui.table(_ranking_df)
    return


@app.cell
def _():
    """Import plotly for visualizations."""
    import plotly.graph_objects as go
    import plotly.express as px
    return (go,)


@app.cell
def _(df, go, mo, ranked_clusters):
    """Visualize cluster rankings."""
    _ranks = list(range(1, len(ranked_clusters) + 1))
    _scores = [score for _, score in ranked_clusters]
    _clusters = [cluster for cluster, _ in ranked_clusters]
    _counts = [len(df[df['cluster'] == c]) for c in _clusters]

    _fig = go.Figure()

    _fig.add_trace(go.Bar(
        x=_ranks,
        y=_scores,
        text=[f"{c[:30]}..." if len(c) > 30 else c for c in _clusters],
        textposition='outside',
        marker=dict(
            color=_scores,
            colorscale='RdYlGn',
            showscale=True,
            colorbar=dict(title="Score")
        ),
        name="Cluster Score"
    ))

    _fig.update_layout(
        title="Cluster Ranking by Biological Relevance Score",
        xaxis_title="Rank",
        yaxis_title="Score",
        showlegend=False,
        template="plotly_white",
        height=600
    )

    mo.ui.plotly(_fig)
    return


@app.cell
def _(df, go, mo, ranked_clusters):
    """Visualize cluster scores vs variant counts."""
    _scores = [score for _, score in ranked_clusters]
    _clusters = [cluster for cluster, _ in ranked_clusters]
    _counts = [len(df[df['cluster'] == c]) for c in _clusters]

    _fig = go.Figure()

    _fig.add_trace(go.Scatter(
        x=_counts,
        y=_scores,
        mode='markers+text',
        text=[c[:20] + "..." if len(c) > 20 else c for c in _clusters],
        textposition='top center',
        marker=dict(
            size=10,
            color=_scores,
            colorscale='RdYlGn',
            showscale=True,
            colorbar=dict(title="Score")
        ),
        name="Clusters"
    ))

    _fig.update_layout(
        title="Cluster Score vs Variant Count",
        xaxis_title="Variant Count",
        yaxis_title="Biological Relevance Score",
        template="plotly_white",
        height=500
    )

    mo.ui.plotly(_fig)
    return


@app.cell
def _(mo, ranked_clusters):
    """Display top clusters summary."""
    _top_n = 10
    _top_clusters = ranked_clusters[:_top_n]

    mo.md(f"""
    ### Top {_top_n} Clusters

    | Rank | Cluster | Score |
    |------|---------|-------|
    """ + "\n".join([
    f"| {i+1} | {cluster} | {score:.2f} |"
    for i, (cluster, score) in enumerate(_top_clusters)
    ]))
    return


@app.cell
def _(mo, np, ranked_clusters):
    """Summary statistics."""
    _top_10_scores = [score for _, score in ranked_clusters[:10]]
    _all_scores = [score for _, score in ranked_clusters]

    mo.md(f"""
    ### Summary Statistics

    - **Top 10 Average Score:** {np.mean(_top_10_scores):.2f}
    - **Overall Average Score:** {np.mean(_all_scores):.2f}
    - **Highest Score:** {max(_all_scores):.2f}
    - **Lowest Score:** {min(_all_scores):.2f}
    - **Score Range:** {max(_all_scores) - min(_all_scores):.2f}
    """)
    return


if __name__ == "__main__":
    app.run()
