#!/usr/bin/env python3
"""
Generate Plotly figures for the ABCA4 report.

This script reads the scored variant matrix, the selected panel, and CRO
work package definitions, then produces a small set of figures that tell
clear, quantitative stories:

- How the selected panel sits in the impact score distribution
- How impact scores relate to population frequency
- How SpliceAI scores distribute, highlighting splice-heavy variants
- How the optimizer trades a little score for much better cluster coverage
- How the panel is distributed across assay modules

Outputs (written to ``data_processed/reports/figures``):
- figure1a_impact.html / .png
- figure1b_af_vs_impact.html / .png
- figure1c_spliceai.html / .png
- figure2_optimization.html / .png
- figure3_mechanisms.html / .png
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
from scipy.stats import gaussian_kde

from ..config import (
    REPORTS_DIR,
    get_scored_variants_path,
    validate_dataframe,
    logger,
)


CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]
FIGURES_DIR = REPORTS_DIR / "figures"
CRO_DIR = CAMPAIGN_ROOT / "data_processed" / "cro"


@dataclass
class SelectionSummary:
    total_variants: int
    k: int
    clusters_covered: int
    total_clusters: int


def _load_base_data() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load scored variants and selected panel with basic validation."""
    scored_path = get_scored_variants_path()
    all_variants = pd.read_parquet(scored_path)
    validate_dataframe(all_variants, "scored_variants", required_columns=["variant_id", "impact_score"])

    selected_path = REPORTS_DIR / "variants_selected.csv"
    if not selected_path.exists():
        raise FileNotFoundError(
            f"Selected variants not found at {selected_path}. "
            "Run the optimization step before plotting."
        )

    selected = pd.read_csv(selected_path)
    validate_dataframe(
        selected,
        "selected_variants",
        required_columns=["variant_id", "impact_score", "cluster_id", "rank"],
    )

    logger.info("Loaded %s scored variants and %s selected variants", len(all_variants), len(selected))
    return all_variants, selected


def _prefix_unique_counts(values: List[str]) -> List[int]:
    """Return cumulative count of unique values as we walk the list."""
    counts: List[int] = []
    seen: set[str] = set()
    for item in values:
        seen.add(str(item))
        counts.append(len(seen))
    return counts


def _save_figure(fig: go.Figure, basename: str) -> Dict[str, str]:
    """Save Plotly figure to HTML and PNG, returning relative paths."""
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    html_path = FIGURES_DIR / f"{basename}.html"
    png_path = FIGURES_DIR / f"{basename}.png"

    logger.info("Writing figure %s", basename)
    # Interactive HTML
    pio.write_html(fig, file=str(html_path), include_plotlyjs="cdn", full_html=True)

    # Static PNG (best effort; kaleido may not be installed)
    try:
        # High-res export to avoid blurry embeds in the HTML report
        pio.write_image(fig, str(png_path), width=1800, height=1200, scale=3)
    except Exception as exc:  # pragma: no cover - visualization environment dependent
        logger.warning("Could not write PNG for %s: %s", basename, exc)

    # Manifest uses paths relative to reports directory
    rel_html = f"figures/{basename}.html"
    rel_png = f"figures/{basename}.png"
    return {"html": rel_html, "png": rel_png}


def _make_impact_distribution(all_variants: pd.DataFrame, selected: pd.DataFrame) -> go.Figure:
    """Tell the impact story: KDE overlay + ECDF + tail hist."""
    impact_all = all_variants["impact_score"].fillna(0.0).clip(lower=0.0)
    impact_sel = selected["impact_score"].fillna(0.0).clip(lower=0.0)

    fig = make_subplots(
        rows=3,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.08,
        row_heights=[0.4, 0.22, 0.38],
        subplot_titles=(
            "Distribution: all vs panel (KDE overlay with medians)",
            "ECDF: how quickly each group accumulates mass",
            "Tail (impact ≥ 0.10): where decisions happen",
        ),
    )

    # KDE helper
    def _kde(series: pd.Series, xmin: float, xmax: float, gridsize: int = 400):
        xs = np.linspace(xmin, xmax, gridsize)
        kde = gaussian_kde(series.values, bw_method="scott")
        ys = kde(xs)
        return xs, ys

    xmax = max(impact_all.max(), 0.6)
    xs_all, ys_all = _kde(impact_all, 0, xmax)
    xs_sel, ys_sel = _kde(impact_sel, 0, xmax)

    fig.add_trace(
        go.Scatter(
            x=xs_all,
            y=ys_all,
            name="All variants",
            mode="lines",
            line=dict(color="#9EA5B5", width=3),
            fill="tozeroy",
            fillcolor="rgba(158,165,181,0.3)",
            hovertemplate="All<br>Impact=%{x:.3f}<br>Density=%{y:.3f}<extra></extra>",
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=xs_sel,
            y=ys_sel,
            name="Selected panel",
            mode="lines",
            line=dict(color="#2D5A3D", width=3),
            fill="tozeroy",
            fillcolor="rgba(45,90,61,0.35)",
            hovertemplate="Panel<br>Impact=%{x:.3f}<br>Density=%{y:.3f}<extra></extra>",
        ),
        row=1,
        col=1,
    )

    # ECDF helper
    def _ecdf(arr: pd.Series) -> tuple[np.ndarray, np.ndarray]:
        vals = np.sort(arr.values)
        n = len(vals)
        y = np.arange(1, n + 1) / n
        return vals, y

    x_all, y_all = _ecdf(impact_all)
    x_sel, y_sel = _ecdf(impact_sel)

    fig.add_trace(
        go.Scatter(
            x=x_all,
            y=y_all,
            name="All variants (ECDF)",
            mode="lines",
            line=dict(color="#9EA5B5", width=3, dash="dash"),
            hovertemplate="Impact=%{x:.3f}<br>Fraction ≤ x=%{y:.2%}<extra></extra>",
        ),
        row=2,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=x_sel,
            y=y_sel,
            name="Selected panel (ECDF)",
            mode="lines",
            line=dict(color="#2D5A3D", width=3),
            hovertemplate="Impact=%{x:.3f}<br>Fraction ≤ x=%{y:.2%}<extra></extra>",
        ),
        row=2,
        col=1,
    )

    # Tail histogram (>=0.10)
    tail_thresh = 0.10
    tail_all = impact_all[impact_all >= tail_thresh]
    tail_sel = impact_sel[impact_sel >= tail_thresh]
    fig.add_trace(
        go.Histogram(
            x=tail_all,
            nbinsx=35,
            name="All variants (tail)",
            marker=dict(color="#CBD2E0"),
            opacity=0.45,
            hovertemplate="All tail<br>Impact=%{x:.3f}<br>Count=%{y}<extra></extra>",
        ),
        row=3,
        col=1,
    )
    fig.add_trace(
        go.Histogram(
            x=tail_sel,
            nbinsx=35,
            name="Selected panel (tail)",
            marker=dict(color="#2D5A3D"),
            opacity=0.75,
            hovertemplate="Selected tail<br>Impact=%{x:.3f}<br>Count=%{y}<extra></extra>",
        ),
        row=3,
        col=1,
    )

    # Percentile markers
    p90_all = float(impact_all.quantile(0.9))
    median_sel = float(impact_sel.median())
    for xpos, label, color in [
        (p90_all, "p90 all", "#666"),
        (median_sel, "median selected", "#2D5A3D"),
    ]:
        fig.add_vline(
            x=xpos,
            line_dash="dot",
            line_color=color,
            row=1,
            col=1,
        )
        fig.add_vline(x=xpos, line_dash="dot", line_color=color, row=2, col=1)
        fig.add_vline(x=xpos, line_dash="dot", line_color=color, row=3, col=1)

    # Effect size annotation
    med_all = impact_all.median()
    med_sel = impact_sel.median()
    lift = med_sel / med_all if med_all > 0 else float("inf")
    fig.add_annotation(
        x=med_sel,
        y=max(ys_sel) * 0.75,
        text=f"Median lift: {lift:.1f}× (panel {med_sel:.3f} vs all {med_all:.3f})",
        showarrow=True,
        arrowhead=2,
        ax=80,
        ay=-40,
        font=dict(size=13, color="#2D5A3D"),
    )

    fig.update_layout(
        template="simple_white",
        barmode="overlay",
        title="Impact score: shift of selected panel vs all candidates",
        legend=dict(orientation="h", yanchor="bottom", y=1.12, xanchor="right", x=1.0),
        margin=dict(l=90, r=40, t=120, b=90),
        font=dict(size=16, color="#222"),
        height=1100,
    )
    fig.update_xaxes(showline=True, linewidth=1.5, linecolor="#333", mirror=True, row=1, col=1, title="")
    fig.update_xaxes(showline=True, linewidth=1.5, linecolor="#333", mirror=True, row=2, col=1, tickformat=".2f", title="")
    fig.update_xaxes(
        showline=True,
        linewidth=1.5,
        linecolor="#333",
        mirror=True,
        row=3,
        col=1,
        range=[tail_thresh, max(impact_all.max(), 0.6)],
        title="Impact score",
    )
    fig.update_yaxes(showline=True, linewidth=1.5, linecolor="#333", mirror=True, rangemode="tozero", row=1, col=1, title="Density (KDE)")
    fig.update_yaxes(showline=True, linewidth=1.5, linecolor="#333", mirror=True, rangemode="tozero", row=2, col=1, tickformat=".0%", title="Fraction ≤ x")
    fig.update_yaxes(showline=True, linewidth=1.5, linecolor="#333", mirror=True, rangemode="tozero", row=3, col=1, title="Count")
    return fig


def _make_af_vs_impact(all_variants: pd.DataFrame, selected: pd.DataFrame) -> go.Figure:
    """Density + overlay of impact score vs gnomAD AF, on log x."""
    df_all = all_variants.copy()
    df_sel = selected.copy()

    # Replace missing or zero AF with a small floor to put them on the log axis
    floor = 1e-6
    df_all["gnomad_max_af_plot"] = df_all["gnomad_max_af"].fillna(floor).clip(lower=floor)
    df_sel["gnomad_max_af_plot"] = df_sel["gnomad_max_af"].fillna(floor).clip(lower=floor)

    af_cut = 1e-4
    impact_cut = 0.2

    fig = make_subplots(
        rows=2,
        cols=2,
        shared_xaxes=False,
        shared_yaxes=False,
        specs=[[{"type": "xy", "colspan": 2}, None], [{"type": "xy"}, {"type": "xy"}]],
        column_widths=[0.7, 0.3],
        row_heights=[0.7, 0.3],
        horizontal_spacing=0.07,
        vertical_spacing=0.09,
        subplot_titles=(
            "Impact score vs gnomAD AF (density + panel overlay)",
            "AF marginal",
            "Impact marginal",
        ),
    )

    # Background density as hexbin-like heatmap
    hist = go.Histogram2d(
        x=df_all["gnomad_max_af_plot"],
        y=df_all["impact_score"],
        nbinsx=55,
        nbinsy=55,
        colorscale="Greys",
        showscale=False,
        opacity=0.55,
        name="All variants density",
        hovertemplate="AF=%{x:.1e}<br>Impact=%{y:.3f}<br>Count=%{z}<extra></extra>",
    )
    fig.add_trace(hist, row=1, col=1)

    # Selected overlay
    fig.add_trace(
        go.Scatter(
            x=df_sel["gnomad_max_af_plot"],
            y=df_sel["impact_score"],
            mode="markers",
            name="Selected panel",
            marker=dict(size=11, color="#2D5A3D", line=dict(color="white", width=1.5)),
            hovertemplate="Selected<br>AF=%{x:.1e}<br>Impact=%{y:.3f}<extra></extra>",
        ),
        row=1,
        col=1,
    )

    # Highlight the ultra-rare, higher-impact region
    y_max = float(max(df_all["impact_score"].max(), df_sel["impact_score"].max(), 0.5))
    fig.add_shape(
        type="rect",
        x0=floor,
        x1=af_cut,
        y0=impact_cut,
        y1=y_max,
        fillcolor="rgba(45, 90, 61, 0.08)",
        line=dict(width=0),
        layer="below",
    )
    fig.add_annotation(
        x=floor * 5,
        y=y_max,
        text="Ultra-rare, higher-impact region",
        showarrow=False,
        xanchor="left",
        yanchor="top",
        font=dict(size=11, color="#2D5A3D"),
    )

    # Marginals
    fig.add_trace(
        go.Histogram(
            x=df_all["gnomad_max_af_plot"],
            nbinsx=50,
            name="AF all",
            marker=dict(color="#CBD2E0"),
            opacity=0.5,
            showlegend=False,
        ),
        row=2,
        col=1,
    )
    fig.add_trace(
        go.Histogram(
            x=df_sel["gnomad_max_af_plot"],
            nbinsx=30,
            name="AF selected",
            marker=dict(color="#2D5A3D"),
            opacity=0.7,
            showlegend=False,
        ),
        row=2,
        col=1,
    )

    fig.add_trace(
        go.Histogram(
            y=df_all["impact_score"],
            nbinsy=40,
            name="Impact all",
            marker=dict(color="#CBD2E0"),
            opacity=0.5,
            showlegend=False,
            orientation="h",
        ),
        row=2,
        col=2,
    )
    fig.add_trace(
        go.Histogram(
            y=df_sel["impact_score"],
            nbinsy=25,
            name="Impact selected",
            marker=dict(color="#2D5A3D"),
            opacity=0.7,
            showlegend=False,
            orientation="h",
        ),
        row=2,
        col=2,
    )

    # Box count inside rare/high box
    panel_inside = ((df_sel["gnomad_max_af_plot"] <= af_cut) & (df_sel["impact_score"] >= impact_cut)).sum()
    fig.add_annotation(
        x=floor * 5,
        y=y_max,
        text=f"{panel_inside} panel variants in shaded region",
        showarrow=False,
        xanchor="left",
        yanchor="top",
        font=dict(size=12, color="#2D5A3D"),
    )

    fig.update_layout(
        template="simple_white",
        title="Population frequency vs impact score (density + panel overlay)",
        legend=dict(orientation="h", yanchor="bottom", y=1.08, xanchor="right", x=1.0),
        margin=dict(l=100, r=40, t=110, b=90),
        font=dict(size=16, color="#222"),
        height=950,
    )
    fig.update_xaxes(type="log", showline=True, linewidth=1.5, linecolor="#333", mirror=True, row=1, col=1, title="gnomAD max allele frequency (log)")
    fig.update_yaxes(range=[0, y_max], showline=True, linewidth=1.5, linecolor="#333", mirror=True, row=1, col=1, title="Impact score")
    fig.update_xaxes(showline=True, linewidth=1.2, linecolor="#333", mirror=True, row=2, col=1, title="AF distribution")
    fig.update_yaxes(showline=True, linewidth=1.2, linecolor="#333", mirror=True, row=2, col=1, title="Count")
    fig.update_xaxes(showline=True, linewidth=1.2, linecolor="#333", mirror=True, row=2, col=2, title="Count")
    fig.update_yaxes(showline=True, linewidth=1.2, linecolor="#333", mirror=True, row=2, col=2, title="Impact distribution")
    return fig


def _make_spliceai_distribution(all_variants: pd.DataFrame, selected: pd.DataFrame) -> go.Figure:
    """SpliceAI score distribution, focusing on high scoring splice candidates."""
    if "spliceai_max_score" not in all_variants.columns:
        raise ValueError("spliceai_max_score column missing from scored variants")

    splice_all = all_variants["spliceai_max_score"].fillna(0.0)
    splice_sel = selected["spliceai_max_score"].fillna(0.0)

    fig = go.Figure()
    fig.add_histogram(
        x=splice_all,
        nbinsx=50,
        name="All variants",
        marker=dict(color="#CBD2E0"),
        opacity=0.45,
        hovertemplate="All variants<br>SpliceAI=%{x:.2f}<br>Count=%{y}<extra></extra>",
    )
    fig.add_histogram(
        x=splice_sel,
        nbinsx=50,
        name="Selected panel",
        marker=dict(color="#2D5A3D"),
        opacity=0.75,
        hovertemplate="Selected panel<br>SpliceAI=%{x:.2f}<br>Count=%{y}<extra></extra>",
    )

    fig.add_vline(
        x=0.5,
        line_dash="dot",
        line_color="#666666",
        annotation_text="High splice impact ~0.5+",
        annotation_position="top right",
    )

    fig.update_layout(
        template="simple_white",
        barmode="overlay",
        title="SpliceAI scores: candidates vs selected panel",
        xaxis_title="SpliceAI max score",
        yaxis_title="Number of variants",
        legend=dict(orientation="h", yanchor="bottom", y=1.05, xanchor="right", x=1.0),
        margin=dict(l=70, r=30, t=80, b=70),
        font=dict(size=16, color="#222"),
    )
    fig.update_xaxes(range=[0, 1], showline=True, linewidth=1.5, linecolor="#333", mirror=True)
    fig.update_yaxes(showline=True, linewidth=1.5, linecolor="#333", mirror=True, rangemode="tozero")
    return fig


def _make_optimization_tradeoff(all_variants: pd.DataFrame, selected: pd.DataFrame) -> go.Figure:
    """
    Compare naive top-k selection with the coverage-aware panel.

    Row 1: cumulative impact score by rank
    Row 2: distinct clusters covered as panel grows
    """
    k = len(selected)
    all_sorted = all_variants.sort_values("impact_score", ascending=False).reset_index(drop=True)
    naive = all_sorted.head(k).copy()

    panel = selected.sort_values("rank").reset_index(drop=True)

    naive["cumulative_impact"] = naive["impact_score"].cumsum()
    panel["cumulative_impact"] = panel["impact_score"].cumsum()

    naive_clusters = _prefix_unique_counts(list(naive["cluster_id"]))
    panel_clusters = _prefix_unique_counts(list(panel["cluster_id"]))

    ranks = list(range(1, k + 1))

    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.08,
        subplot_titles=("Cumulative impact score", "Distinct clusters covered"),
    )

    fig.add_trace(
        go.Scatter(
            x=ranks,
            y=naive["cumulative_impact"],
            name="Naive top-k (score only)",
            line=dict(color="#9EA5B5", dash="dash", width=3),
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=ranks,
            y=panel["cumulative_impact"],
            name="Coverage-aware panel",
            line=dict(color="#2D5A3D", width=3),
        ),
        row=1,
        col=1,
    )

    fig.add_trace(
        go.Scatter(
            x=ranks,
            y=naive_clusters,
            name="Naive clusters",
            line=dict(color="#9EA5B5", dash="dot", width=3),
        ),
        row=2,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=ranks,
            y=panel_clusters,
            name="Panel clusters",
            line=dict(color="#2D5A3D", width=3),
        ),
        row=2,
        col=1,
    )

    fig.update_xaxes(title_text="Variant rank in panel", row=2, col=1)
    fig.update_yaxes(title_text="Cumulative impact score", row=1, col=1)
    fig.update_yaxes(title_text="Distinct clusters", row=2, col=1)

    fig.update_layout(
        template="simple_white",
        title="Optimization tradeoff: score vs cluster coverage",
        legend=dict(orientation="h", yanchor="bottom", y=-0.06, xanchor="right", x=1.0),
        margin=dict(l=90, r=30, t=80, b=90),
        font=dict(size=16, color="#222"),
    )
    fig.update_xaxes(showline=True, linewidth=1.5, linecolor="#333", mirror=True)
    fig.update_yaxes(showline=True, linewidth=1.5, linecolor="#333", mirror=True, rangemode="tozero")

    # Headline annotation
    score_retained = panel["impact_score"].sum() / naive["impact_score"].sum() if naive["impact_score"].sum() else 0
    coverage_gain = panel_clusters[-1] - naive_clusters[-1]
    fig.add_annotation(
        x=ranks[-1] * 0.65,
        y=naive["cumulative_impact"].iloc[-1] * 0.9,
        text=f"Keeps {score_retained:.0%} of score; +{coverage_gain} clusters vs naive",
        showarrow=True,
        arrowhead=2,
        ax=40,
        ay=-40,
        font=dict(size=13, color="#2D5A3D"),
    )
    return fig


def _load_cro_work_packages() -> List[Dict[str, object]]:
    """Load CRO work packages from JSONL, if available."""
    wp_path = CRO_DIR / "work_packages.jsonl"
    if not wp_path.exists():
        logger.warning("CRO work package file not found at %s; skipping Figure 3", wp_path)
        return []

    work_packages: List[Dict[str, object]] = []
    with wp_path.open() as handle:
        for line in handle:
            if line.strip():
                work_packages.append(json.loads(line))
    return work_packages


def _make_mechanism_coverage(work_packages: List[Dict[str, object]]) -> go.Figure:
    """Heatmap of cluster × assay module plus totals strip."""
    if not work_packages:
        raise ValueError("No work packages provided for mechanism coverage figure")

    # Build module totals and per-variant mapping
    module_map: dict[str, list[str]] = {}
    for wp in work_packages:
        module = str(wp.get("assay_module", "unknown"))
        variant_ids = wp.get("variant_ids", [])
        module_map.setdefault(module, []).extend(variant_ids)

    # Load selected variants with cluster_id
    selected_path = REPORTS_DIR / "variants_selected.csv"
    sel_df = pd.read_csv(selected_path)

    def _to_cro_format(vid: str) -> str:
        parts = vid.split("_")
        if len(parts) == 4:
            return f"{parts[0]}:{parts[1]}:{parts[2]}/{parts[3]}"
        return vid

    sel_df["variant_id_cro"] = sel_df["variant_id"].apply(_to_cro_format)

    records = []
    for module, var_ids in module_map.items():
        for vid in var_ids:
            cluster = sel_df.loc[sel_df["variant_id_cro"] == vid, "cluster_id"]
            if not cluster.empty:
                records.append({"module": module, "cluster_id": cluster.iloc[0]})

    if not records:
        # Fallback to module totals bar if join failed
        modules = list(module_map.keys())
        counts = [len(v) for v in module_map.values()]
        fig = go.Figure()
        fig.add_bar(
            x=modules,
            y=counts,
            marker=dict(color="#2D5A3D"),
            text=counts,
            textposition="outside",
        )
        fig.update_layout(
            template="simple_white",
            title="Assay module coverage for selected variants",
            xaxis_title="Assay module",
            yaxis_title="Variants per module",
            margin=dict(l=70, r=20, t=70, b=60),
        )
        return fig

    cov_df = pd.DataFrame(records)
    # Limit to top clusters by count to keep chart legible
    top_clusters = cov_df["cluster_id"].value_counts().head(12).index.tolist()
    cov_df = cov_df[cov_df["cluster_id"].isin(top_clusters)]

    pivot = (
        cov_df.assign(count=1)
        .pivot_table(index="cluster_id", columns="module", values="count", aggfunc="sum", fill_value=0)
    )
    pivot = pivot.sort_index()

    modules = pivot.columns.tolist()
    clusters = pivot.index.tolist()
    z = pivot.values

    # Totals strip (module totals)
    module_totals = [len(module_map.get(m, [])) for m in modules]

    fig = make_subplots(
        rows=2,
        cols=1,
        row_heights=[0.82, 0.18],
        vertical_spacing=0.08,
        shared_xaxes=True,
        specs=[[{"type": "heatmap"}], [{"type": "bar"}]],
    )

    fig.add_trace(
        go.Heatmap(
            z=z,
            x=modules,
            y=clusters,
            colorscale="Greens",
            colorbar=dict(title="Count"),
            hovertemplate="Module=%{x}<br>Cluster=%{y}<br>Variants=%{z}<extra></extra>",
        ),
        row=1,
        col=1,
    )

    fig.add_trace(
        go.Bar(
            x=modules,
            y=module_totals,
            marker=dict(color="#2D5A3D"),
            text=module_totals,
            textposition="outside",
            name="Module totals",
            hovertemplate="Module=%{x}<br>Total variants=%{y}<extra></extra>",
        ),
        row=2,
        col=1,
    )

    fig.update_layout(
        template="simple_white",
        title="Mechanistic coverage: clusters × assay modules (top 12 clusters)",
        margin=dict(l=120, r=40, t=80, b=80),
        font=dict(size=16, color="#222"),
        height=900,
        showlegend=False,
    )
    fig.update_yaxes(title="Cluster", automargin=True, row=1, col=1)
    fig.update_xaxes(title="", row=1, col=1)
    fig.update_xaxes(title="Assay module totals", row=2, col=1)
    fig.update_yaxes(title="Variants", row=2, col=1)
    return fig


def generate_figures_and_manifest() -> Dict[str, object]:
    """
    Generate all figures and write an updated manifest.

    Returns manifest dict for downstream HTML rendering.
    """
    all_variants, selected = _load_base_data()

    selection_summary = SelectionSummary(
        total_variants=int(len(all_variants)),
        k=int(len(selected)),
        clusters_covered=int(selected["cluster_id"].nunique()),
        total_clusters=int(all_variants["cluster_id"].nunique()),
    )

    figures: Dict[str, Dict[str, str]] = {}

    fig1a = _make_impact_distribution(all_variants, selected)
    figures["figure1a_impact"] = _save_figure(fig1a, "figure1a_impact")

    fig1b = _make_af_vs_impact(all_variants, selected)
    figures["figure1b_af_vs_impact"] = _save_figure(fig1b, "figure1b_af_vs_impact")

    fig1c = _make_spliceai_distribution(all_variants, selected)
    figures["figure1c_spliceai"] = _save_figure(fig1c, "figure1c_spliceai")

    fig2 = _make_optimization_tradeoff(all_variants, selected)
    figures["figure2_optimization"] = _save_figure(fig2, "figure2_optimization")

    work_packages = _load_cro_work_packages()
    if work_packages:
        fig3 = _make_mechanism_coverage(work_packages)
        figures["figure3_mechanisms"] = _save_figure(fig3, "figure3_mechanisms")

    manifest: Dict[str, object] = {
        "metadata": {
            "date": date.today().isoformat(),
            "campaign_name": "ABCA4 Variant Triage",
        },
        "total_variants": selection_summary.total_variants,
        "k": selection_summary.k,
        "clusters_covered": selection_summary.clusters_covered,
        "total_clusters": selection_summary.total_clusters,
        "work_packages": [
            {
                "wp_id": wp.get("wp_id"),
                "assay_module": wp.get("assay_module"),
                "variant_count": len(wp.get("variant_ids", [])),
            }
            for wp in work_packages
        ],
        "figures": [paths["png"] for paths in figures.values()],
    }

    manifest_path = REPORTS_DIR / "abca4_manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)

    logger.info("Wrote manifest to %s", manifest_path)
    return manifest


def main() -> None:
    """CLI entry point."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(name)s | %(levelname)s | %(message)s",
    )
    generate_figures_and_manifest()


if __name__ == "__main__":
    main()
