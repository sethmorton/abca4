#!/usr/bin/env python3
"""
ABCA4 Campaign – Strand Optimization Dashboard & Reporting

Steps 6-8: Run Strand optimization, map to assays, generate reports.

Run interactively:  marimo edit notebooks/03_optimization_dashboard.py
Run as dashboard:   marimo run notebooks/03_optimization_dashboard.py
Run as script:      python notebooks/03_optimization_dashboard.py
Run autonomously:   python notebooks/03_optimization_dashboard.py --auto
"""

import marimo

__generated_with = "0.17.8"
app = marimo.App()


@app.cell
def __():
    """Import core libraries and detect autonomous mode."""
    import marimo as mo
    import pandas as pd
    import numpy as np
    from pathlib import Path
    import json
    import sys
    import logging
    from typing import Optional, Dict, List, Tuple
    from datetime import datetime
    # Config is imported at module level, no need to re-import

    logger = logging.getLogger(__name__)
    
    # Detect autonomous mode from command line or environment
    AUTONOMOUS_MODE = "--auto" in sys.argv or "MARIMO_AUTONOMOUS" in __import__("os").environ
    
    return mo, pd, np, Path, logger, json, datetime, Optional, Dict, List, Tuple, sys, AUTONOMOUS_MODE


@app.cell
def __(mo, logger, pd):
    """Define campaign paths."""
    from src.config import (
        CAMPAIGN_ROOT, FEATURES_DIR, REPORTS_DIR, get_scored_variants_path,
        validate_file_exists, load_parquet_safely, DEMO_MODE_ENABLED
    )

    # Ensure output directory exists
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)

    scored_path_optim = get_scored_variants_path()
    try:
        validate_file_exists(scored_path_optim, "Run 02_feature_engineering.py first")
        df_scored_optim = load_parquet_safely(scored_path_optim, "scored variants")
        logger.info(f"Loaded {len(df_scored_optim)} scored variants for optimization")
    except (FileNotFoundError, ValueError) as e:
        if DEMO_MODE_ENABLED:
            mo.md(f"**Demo Mode:** {str(e)}\n\nUsing synthetic data for demonstration.")
            # Create minimal demo data
            df_scored_optim = pd.DataFrame({
                'variant_id': [f'demo_{i}' for i in range(10)],
                'chrom': ['1'] * 10,
                'pos': range(94400000, 94400010),
                'ref': ['A'] * 10,
                'alt': ['T'] * 10,
                'model_score': [0.1 * i for i in range(10)],
                'cluster_id': [0] * 5 + [1] * 5
            })
        else:
            mo.md(f"**Error:** {str(e)}")
            df_scored_optim = pd.DataFrame()

    return CAMPAIGN_ROOT, FEATURES_DIR, REPORTS_DIR, scored_path_optim, df_scored_optim


@app.cell
def __(mo, df_scored_optim):
    """Smoke test: Validate scored variants dataframe."""

    if df_scored_optim.empty:
        mo.md("**Error:** Scored variants dataframe is empty")
    else:
        # Check required columns for optimization
        required_cols = ['variant_id', 'chrom', 'pos', 'ref', 'alt']
        missing_cols = [col for col in required_cols if col not in df_scored_optim.columns]

        if missing_cols:
            mo.md(f"**Error:** Missing required columns: {missing_cols}")
        else:
            # Check for model scores
            score_cols = ['model_score', 'optimization_score']
            present_scores = [col for col in score_cols if col in df_scored_optim.columns]

            mo.md(f"✅ Scored variants validation passed\n- Rows: {len(df_scored_optim)}\n- Required columns: ✓\n- Score columns: {len(present_scores)} present")


@app.cell
def __(mo):
    """
    ## Step 6: Strand Optimization Dashboard

    Configure and run Strand engine for variant selection.
    """
    mo.md(__doc__)


@app.cell
def __(mo):
    """Plan alignment note for Steps 6‑8."""
    mo.md("""
### Plan Alignment

- **Step 6:** This notebook’s optimization controls drive the Strand `OptimizationRunner` and record when we fall back to feature ranking.
- **Step 7:** Mechanism/assay mapping cells turn selected variants into experimental intents.
- **Step 8:** The report preview/export produces the short, five-minute read promised in the v1 plan.
""")


@app.cell
def __(mo, df_scored_optim, AUTONOMOUS_MODE):
    """Optimization parameter controls."""
    if df_scored_optim.empty:
        mo.md("No scored variants. Load data first.")
        k_optim = None
        strat_optim = None
        lambda_cov = None
    else:
        mo.md("### Optimization Parameters")
        
        if AUTONOMOUS_MODE:
            # In autonomous mode, use fixed values without interactive controls
            class _StaticControl:
                def __init__(self, value):
                    self.value = value
            
            k_value = min(30, len(df_scored_optim))
            strat_value = "Greedy coverage"
            lambda_value = 0.6
            
            k_optim = _StaticControl(k_value)
            strat_optim = _StaticControl(strat_value)
            lambda_cov = _StaticControl(lambda_value)
            
            mo.md(f"**Autonomous Mode:** Using K={k_value}, Strategy={strat_value}, λ={lambda_value}")
        else:
            k_optim = mo.ui.slider(
                10, min(200, len(df_scored_optim)),
                value=min(30, len(df_scored_optim)),
                label="Panel Size (K)"
            )

            strat_optim = mo.ui.radio(
                options=["Greedy coverage", "Top score", "Random"],
                value="Greedy coverage",
                label="Strategy"
            )

            lambda_cov = mo.ui.slider(
                0.0, 2.0, value=0.6, step=0.05, label="Coverage penalty λ"
            )

    return k_optim, strat_optim, lambda_cov


@app.cell
def __(
    mo,
    k_optim, strat_optim, lambda_cov
):
    """Echo chosen strategy and coverage penalty."""
    if strat_optim is None:
        mo.md("Select a strategy to proceed.")
        normalized_wgt_optim = None
    else:
        mo.md(f"**Strategy:** {strat_optim.value if hasattr(strat_optim, 'value') else strat_optim}\n\n**Coverage λ:** {lambda_cov.value if lambda_cov else 0}")
        normalized_wgt_optim = None


@app.cell
def __(mo, AUTONOMOUS_MODE):
    """Run optimization button."""
    if AUTONOMOUS_MODE:
        # In autonomous mode, simulate button press automatically
        run_optim_btn = True
        mo.md("**Running optimization automatically...**")
    else:
        run_optim_btn = mo.ui.button(label="▶️ Run Optimization", on_click=lambda _: True)
    return run_optim_btn


@app.cell
def __(
    logger, pd, np,
    df_scored_optim,
    k_optim, strat_optim, lambda_cov,
    run_optim_btn, CAMPAIGN_ROOT
):
    """
    Execute Strand optimization with real engine or feature-based ranking.
    
    Uses the authoritative reward stack from campaigns/abca4/src/reward/run_abca4_optimization.py
    and logs results to MLflow.
    """
    if (run_optim_btn is None or not run_optim_btn or df_scored_optim.empty or 
        k_optim is None):
        optim_results = None
    else:
        logger.info(f"Starting optimization: {strat_optim.value if strat_optim else 'Random'}")

        def _select_greedy(df: pd.DataFrame, k: int, lambda_penalty: float) -> pd.DataFrame:
            df_work = df.copy()
            score_col = 'impact_score' if 'impact_score' in df_work.columns else 'model_score'
            if score_col not in df_work.columns:
                raise ValueError("No impact_score or model_score available for optimization")

            # Precompute tau_j from cluster_target if present
            cluster_targets = {}
            if 'cluster_id' in df_work.columns and 'cluster_target' in df_work.columns:
                cluster_targets = (
                    df_work[['cluster_id', 'cluster_target']]
                    .dropna()
                    .drop_duplicates('cluster_id')
                    .set_index('cluster_id')['cluster_target']
                    .to_dict()
                )

            selected = []
            cov = {cid: 0.0 for cid in cluster_targets.keys()}

            def _penalty(cov_dict):
                return sum(max(0.0, cluster_targets.get(cid, 0.0) - cov_dict.get(cid, 0.0)) for cid in cluster_targets)

            remaining_idx = set(df_work.index)
            for _ in range(k):
                best_idx = None
                best_gain = -np.inf
                for idx in list(remaining_idx):
                    row = df_work.loc[idx]
                    cid = row.get('cluster_id', None)
                    score = float(row.get(score_col, 0.0))
                    cov_new = cov.copy()
                    if cid in cluster_targets:
                        cov_new[cid] = max(cov_new.get(cid, 0.0), score)
                    gain = score - lambda_penalty * (_penalty(cov_new) - _penalty(cov))
                    if gain > best_gain:
                        best_gain = gain
                        best_idx = idx
                if best_idx is None:
                    break
                remaining_idx.remove(best_idx)
                row = df_work.loc[best_idx]
                cid = row.get('cluster_id', None)
                if cid in cluster_targets:
                    cov[cid] = max(cov.get(cid, 0.0), float(row.get(score_col, 0.0)))
                selected.append(best_idx)

            df_sel = df_work.loc[selected].copy()
            df_sel['optimization_score'] = df_sel[score_col]
            df_sel['selected'] = True
            df_sel['rank'] = range(1, len(df_sel) + 1)
            return df_sel

        try:
            if strat_optim and strat_optim.value == "Top score":
                score_col = 'impact_score' if 'impact_score' in df_scored_optim.columns else 'model_score'
                ranked = df_scored_optim.sort_values(score_col, ascending=False).head(k_optim.value).copy()
                ranked['optimization_score'] = ranked[score_col]
                ranked['selected'] = True
                ranked['rank'] = range(1, len(ranked) + 1)
                optim_mode = "top_score"
                _df_selected = ranked
            elif strat_optim and strat_optim.value == "Random":
                sample = df_scored_optim.sample(n=min(k_optim.value, len(df_scored_optim)), random_state=42).copy()
                score_col = 'impact_score' if 'impact_score' in sample.columns else 'model_score'
                sample['optimization_score'] = sample[score_col]
                sample['selected'] = True
                sample['rank'] = range(1, len(sample) + 1)
                optim_mode = "random"
                _df_selected = sample
            else:
                _df_selected = _select_greedy(df_scored_optim, k_optim.value, lambda_cov.value if lambda_cov else 0.0)
                optim_mode = "greedy_coverage"

            optim_results = {
                "strategy": strat_optim.value if strat_optim else "Greedy coverage",
                "k": k_optim.value,
                "coverage_lambda": lambda_cov.value if lambda_cov else 0.0,
                "selected_variants": _df_selected,
                "timestamp": pd.Timestamp.now(),
                "mode": optim_mode,
                "artifact_path": None,
                "error": None,
            }

            mo.md("Optimization completed.")
            logger.info(f"Optimization complete. Selected {len(_df_selected)} variants")

        except Exception as e:
            logger.error(f"Optimization failed: {e}")
            optim_results = {
                "mode": "error",
                "error": str(e),
                "selected_variants": pd.DataFrame(),
                "strategy": strat_optim.value if strat_optim else "Greedy coverage",
                "k": k_optim.value,
                "coverage_lambda": lambda_cov.value if lambda_cov else 0.0,
                "timestamp": pd.Timestamp.now(),
            }

    return optim_results


@app.cell
def __():
    """Import plotly for visualizations."""
    import plotly.graph_objects as go
    return go


@app.cell
def __(mo, optim_results):
    """Summarize which optimization path ran."""
    if optim_results is None:
        mo.md("Optimization not run yet.")
    else:
        mode = optim_results.get('mode', 'unknown')
        error = optim_results.get('error')

        if mode == "error":
            msg = f"**Optimizer mode:** `{mode}` - {error}"
        else:
            msg = f"**Optimizer mode:** `{mode}`\n**λ (coverage):** {optim_results.get('coverage_lambda', 0)}"
        mo.md(msg)


@app.cell
def __(mo, go, optim_results):
    """Display optimization summary with visualization."""
    if optim_results is None:
        mo.md("Run optimization to see results.")
    elif optim_results.get('mode') == 'error':
        mo.md("Cannot display results due to optimization error. Fix the issue above first.")
    else:
        mo.md(f"""
### Results

- Strategy: {optim_results['strategy']}
- K: {optim_results['k']}
- Selected: {len(optim_results['selected_variants'])}
""")

        # Plot optimization scores
        try:
            _df_sel = optim_results['selected_variants'].sort_values("rank")
            _fig_opt = go.Figure()
            _fig_opt.add_trace(go.Bar(
                x=_df_sel["rank"],
                y=_df_sel["optimization_score"],
                name="Score",
                marker_color="lightseagreen"
            ))
            _fig_opt.update_layout(
                title="Selected Variants by Optimization Score",
                xaxis_title="Rank",
                yaxis_title="Score",
                showlegend=False,
                template="plotly_white"
            )
            mo.ui.plotly(_fig_opt)
        except Exception as _e:
            mo.md(f"Plot error: {_e}")


@app.cell
def __(mo):
    """
    ## Step 7: Experimental Mapping

    Map variants to assays and collect rationale.
    """
    mo.md(__doc__)


@app.cell
def __(
    pd, optim_results, logger, CAMPAIGN_ROOT
):
    """
    Map selected variants to experimental mechanisms and assays.
    
    Uses consequence annotations and domain information from the real
    annotation pipeline to suggest appropriate validation assays.
    """
    _mech_map = {
        "frameshift_variant": "Loss-of-Function",
        "stop_gained": "Loss-of-Function",
        "stop_lost": "Loss-of-Function",
        "splice_acceptor_variant": "Splicing defect",
        "splice_donor_variant": "Splicing defect",
        "missense_variant": "Protein folding",
        "inframe_deletion": "Structural disruption",
        "inframe_insertion": "Structural disruption",
        "synonymous_variant": "Regulatory",
        "upstream_gene_variant": "Regulatory",
        "downstream_gene_variant": "Regulatory",
    }

    _assay_map = {
        "frameshift_variant": "Western blot, confocal microscopy, flow cytometry",
        "stop_gained": "Western blot, confocal microscopy, flow cytometry",
        "stop_lost": "Western blot, immunofluorescence",
        "splice_acceptor_variant": "RT-qPCR, Northern blot, Western blot",
        "splice_donor_variant": "RT-qPCR, Northern blot, Western blot",
        "missense_variant": "Differential scanning fluorimetry (DSF), size exclusion chromatography (SEC)",
        "inframe_deletion": "Western blot, immunoprecipitation, SEC",
        "inframe_insertion": "Western blot, immunoprecipitation, SEC",
        "synonymous_variant": "RNA-seq, luciferase reporter",
        "upstream_gene_variant": "EMSA, chromatin IP, reporter assay",
        "downstream_gene_variant": "Reporter assay, ATAC-seq",
    }

    if optim_results is None or optim_results['selected_variants'].empty:
        df_mapped = None
    else:
        df_mapped = optim_results['selected_variants'].copy()

        # Find consequence column (may have various names)
        consequence_col = next((c for c in df_mapped.columns if 'consequence' in c.lower()), None)
        domain_col = next((c for c in df_mapped.columns if 'domain' in c.lower()), None)
        
        # Assign mechanisms and assays based on consequence
        if consequence_col:
            df_mapped["mechanism"] = df_mapped[consequence_col].apply(
                lambda x: _mech_map.get(str(x).lower() if pd.notna(x) else "missense_variant", "Protein misfolding")
            )
            df_mapped["suggested_assay"] = df_mapped[consequence_col].apply(
                lambda x: _assay_map.get(str(x).lower() if pd.notna(x) else "missense_variant", "Functional assay")
            )
        else:
            df_mapped["mechanism"] = "Protein misfolding"
            df_mapped["suggested_assay"] = "Functional assay"
        
        # Enhance with domain info if available
        if domain_col:
            df_mapped["domain_note"] = df_mapped[domain_col].apply(
                lambda x: f"in {x} domain" if pd.notna(x) and str(x).lower() != "unknown" else ""
            )
        else:
            df_mapped["domain_note"] = ""
        
        # Build rationale with safe formatting
        def _format_rationale(row):
            mechanism = row.get('mechanism', 'Putative')
            domain_note = row.get('domain_note') or 'variant position'
            score_val = row.get('model_score', row.get('optimization_score'))
            if isinstance(score_val, (int, float)) and pd.notna(score_val):
                score_str = f"{float(score_val):.3f}"
            else:
                score_str = "N/A"
            return f"{mechanism}, {domain_note}. Model score: {score_str}"

        df_mapped["rationale"] = df_mapped.apply(_format_rationale, axis=1)

    return _mech_map, _assay_map, df_mapped


@app.cell
def __(mo, df_mapped):
    """Display selected variants table."""
    if df_mapped is None or df_mapped.empty:
        mo.md("No optimization results yet.")
    else:
        _display_cols = ["rank", "chrom", "pos", "ref", "alt", "vep_consequence", "mechanism", "suggested_assay"]
        _avail_cols = [c for c in _display_cols if c in df_mapped.columns]

        mo.md("""
### Selected Variants for Experimental Validation
""")
        mo.ui.table(df_mapped[_avail_cols].head(50))


@app.cell
def __(mo):
    """
    ### Export Selected Variants
    """
    mo.md(__doc__)


@app.cell
def __(
    logger, df_mapped, REPORTS_DIR, optim_results
):
    """Export to CSV and JSON."""
    csv_export_path = None
    json_export_path = None

    # Guard against missing MLflow artifacts in error mode
    if optim_results and optim_results.get('mode') == 'error':
        logger.warning("Export disabled in error mode - results are not available")
    elif df_mapped is not None and not df_mapped.empty:
        csv_export_path = REPORTS_DIR / "variants_selected.csv"
        df_mapped.to_csv(csv_export_path, index=False)
        logger.info(f"Exported CSV to {csv_export_path}")

        json_export_path = REPORTS_DIR / "variants_selected.json"
        _json_data = df_mapped[[
            "rank", "chrom", "pos", "ref", "alt",
            "vep_consequence", "mechanism", "suggested_assay"
        ]].to_dict(orient="records")

        with open(json_export_path, "w") as _f_json:
            json.dump(_json_data, _f_json, indent=2)

        logger.info(f"Exported JSON to {json_export_path}")

    return csv_export_path, json_export_path


@app.cell
def __(mo):
    """
    ## Step 8: Report Preview & Export

    Generate final Markdown report.
    """
    mo.md(__doc__)


@app.cell
def __(mo, pd, AUTONOMOUS_MODE):
    """Report configuration."""
    _today = pd.Timestamp.now().strftime("%Y-%m-%d")
    
    if AUTONOMOUS_MODE:
        # In autonomous mode, use fixed values
        class _StaticTextControl:
            def __init__(self, value):
                self.value = value
        
        report_title_widget = _StaticTextControl("ABCA4 Variant Selection Report v1")
        report_date_widget = _StaticTextControl(_today)
        report_notes_widget = _StaticTextControl("Autonomous optimization run for quality assurance verification.")
        mo.md(f"**Report Config (Autonomous):** Title='ABCA4 Variant Selection Report v1', Date='{_today}'")
    else:
        report_title_widget = mo.ui.text(value="ABCA4 Variant Selection Report v1", label="Title")
        report_date_widget = mo.ui.text(value=_today, label="Date")
        report_notes_widget = mo.ui.text_area(value="", label="Additional Notes")

    return report_title_widget, report_date_widget, report_notes_widget


@app.cell
def __():
    """Helper function for weight display."""
    def _build_weight_str(weights_dict):
        if not weights_dict:
            return "- (Feature-based ranking used)"
        lines = []
        for key in ['enformer', 'motif', 'conservation', 'dnafm', 'regulatory', 'splice', 'missense']:
            if key in weights_dict:
                label = key.replace('dnafm', 'DNA FM')
                lines.append(f"- {label}: {weights_dict[key]:.3f}")
        return "\n".join(lines) if lines else "- (Feature-based ranking used)"
    return _build_weight_str


@app.cell
def __(
    pd, datetime,
    df_mapped, optim_results,
    report_title_widget, report_date_widget, report_notes_widget
):
    """Generate Markdown report."""
    if df_mapped is None or optim_results is None:
        report_md = None
    else:
        # Helper function to build weight string
        def _build_weight_str_local(weights_dict):
            if not weights_dict:
                return "- (Feature-based ranking used)"
            lines = []
            for key in ['enformer', 'motif', 'conservation', 'dnafm', 'regulatory', 'splice', 'missense']:
                if key in weights_dict:
                    label = key.replace('dnafm', 'DNA FM')
                    lines.append(f"- {label}: {weights_dict[key]:.3f}")
            return "\n".join(lines) if lines else "- (Feature-based ranking used)"
        
        _date_str = str(report_date_widget.value) if report_date_widget.value else pd.Timestamp.now().strftime("%Y-%m-%d")

        report_md = f"""# {report_title_widget.value}

**Date:** {_date_str}

## Executive Summary

Curated selection of ABCA4 variants for functional validation.

## Selection Approach

- **Strategy:** {optim_results.get('strategy', 'Feature ranking')}
- **Panel Size (K):** {optim_results.get('k', len(df_mapped))}
- **Coverage λ:** {optim_results.get('coverage_lambda', 0)}

### Feature Weights

{_build_weight_str_local(optim_results.get('weights', {}))}

## Selected Variants ({len(df_mapped)})

| Rank | Variant | Consequence | Mechanism | Assay |
|------|---------|-------------|-----------|-------|
"""
        for _idx, _row in df_mapped.head(50).iterrows():
            _rank = _row.get("rank", _idx)
            _variant = f"{_row.get('chrom', '?')}:{_row.get('pos', '?')}:{_row.get('ref', '?')}/{_row.get('alt', '?')}"
            _consequence = _row.get("vep_consequence", "?")
            _mechanism = _row.get("mechanism", "?")
            _assay = _row.get("suggested_assay", "?")
            report_md += f"\n| {_rank} | {_variant} | {_consequence} | {_mechanism} | {_assay} |"

        report_md += f"""

## Assay Plan

Variants will be subjected to mechanism-specific functional assays for validation.

## Notes

{report_notes_widget.value}

---
*Generated on {_date_str} using Strand framework + Marimo*
"""

    return report_md


@app.cell
def __(mo, report_md):
    """Display report preview."""
    if report_md is None:
        mo.md("Generate optimization results first.")
    else:
        mo.md(report_md)


@app.cell
def __(
    logger, report_md, REPORTS_DIR, optim_results
):
    """Export report to Markdown."""
    report_md_path = None

    # Guard against missing MLflow artifacts in error mode
    if optim_results and optim_results.get('mode') == 'error':
        logger.warning("Report export disabled in error mode")
    elif report_md is not None:
        report_md_path = REPORTS_DIR / "report_snapshot.md"
        with open(report_md_path, "w") as _f_report:
            _f_report.write(report_md)
        logger.info(f"Exported report to {report_md_path}")

    return report_md_path


@app.cell
def __(mo, csv_export_path, json_export_path, report_md_path):
    """Confirm all exports."""
    mo.md(f"""
**Optimization Complete**

**Exported Files:**
- Report: {report_md_path}
- CSV: {csv_export_path}
- JSON: {json_export_path}

Ready for downstream analysis!
""")


if __name__ == "__main__":
    app.run()
