#!/usr/bin/env python3
"""
ABCA4 Campaign – Strand Optimization Dashboard & Reporting

Steps 6-8: Run Strand optimization, map to assays, generate reports.

Run interactively:  marimo edit notebooks/03_optimization_dashboard.py
Run as dashboard:   marimo run notebooks/03_optimization_dashboard.py
Run as script:      python notebooks/03_optimization_dashboard.py
"""

import marimo

__generated_with = "0.17.8"
app = marimo.App()


@app.cell
def __():
    """Import core libraries."""
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

    return mo, pd, np, Path, logger, json, datetime, Optional, Dict, List, Tuple, sys


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
def __(mo, df_scored_optim):
    """Optimization parameter controls."""
    if df_scored_optim.empty:
        mo.md("No scored variants. Load data first.")
        k_optim = None
        iters_optim = None
        strat_optim = None
        wgt_enformer = None
        wgt_motif = None
        wgt_cons = None
        wgt_dna = None
    else:
        mo.md("### Optimization Parameters")

        k_optim = mo.ui.slider(
            10, min(200, len(df_scored_optim)),
            value=min(30, len(df_scored_optim)),
            label="Panel Size (K)"
        )

        iters_optim = mo.ui.slider(
            100, 5000,
            value=1000,
            step=100,
            label="Iterations"
        )

        strat_optim = mo.ui.radio(
            options=["CEM", "GA", "Random"],
            value="CEM",
            label="Strategy"
        )

        mo.md("### Reward Weights")
        wgt_enformer = mo.ui.slider(0, 1, value=0.4, step=0.05, label="Enformer Δ")
        wgt_motif = mo.ui.slider(0, 1, value=0.3, step=0.05, label="Motif Δ")
        wgt_cons = mo.ui.slider(0, 1, value=0.2, step=0.05, label="Conservation")
        wgt_dna = mo.ui.slider(0, 1, value=0.1, step=0.05, label="DNA FM Δ")

    return k_optim, iters_optim, strat_optim, wgt_enformer, wgt_motif, wgt_cons, wgt_dna


@app.cell
def __(
    mo,
    k_optim, iters_optim, strat_optim,
    wgt_enformer, wgt_motif, wgt_cons, wgt_dna
):
    """Display normalized weights."""
    if strat_optim is None:
        normalized_wgt_optim = None
    else:
        _total = (wgt_enformer.value + wgt_motif.value + wgt_cons.value + wgt_dna.value)
        if _total == 0:
            _total = 1.0

        normalized_wgt_optim = {
            "enformer": wgt_enformer.value / _total,
            "motif": wgt_motif.value / _total,
            "conservation": wgt_cons.value / _total,
            "dnafm": wgt_dna.value / _total,
        }

        _wgt_df = pd.DataFrame([
            {"Component": k, "Weight": f"{v:.3f}"}
            for k, v in normalized_wgt_optim.items()
        ])
        mo.md("### Normalized Weights")
        mo.ui.table(_wgt_df)


@app.cell
def __(mo):
    """Run optimization button."""
    run_optim_btn = mo.ui.button(label="▶️ Run Optimization", on_click=lambda _: True)
    return run_optim_btn


@app.cell
def __(
    logger, pd, np,
    df_scored_optim,
    k_optim, iters_optim, strat_optim, normalized_wgt_optim,
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
        runner = None
        try:
            # Try to use the Strand engine if available
            try:
                sys.path.insert(0, str(CAMPAIGN_ROOT))
                from src.reward.run_abca4_optimization import OptimizationRunner

                runner = OptimizationRunner()
                logger.info("Using real Strand optimization via OptimizationRunner")
                
                # Build a mock feature matrix from scored variants
                feature_matrix = df_scored_optim.copy()
                sequences = runner.build_sequences(feature_matrix)
                ranked = runner.compute_scores(sequences)
                
                # Log to MLflow
                runner.log_mlflow(ranked, top_k=k_optim.value)
                
                # Select top-k
                _df_selected = ranked.head(k_optim.value).copy()
                _df_selected["selected"] = True
                _df_selected["rank"] = range(1, len(_df_selected) + 1)

                # Require deterministic reward scores - no random fallbacks
                if "reward" not in _df_selected.columns:
                    raise ValueError("Strand optimization runner did not return 'reward' column - scores would be non-deterministic")

                _df_selected["optimization_score"] = _df_selected["reward"]
                optim_mode = "strand"

            except Exception as e:
                optimization_error = f"Strand optimization engine unavailable: {e}"
                logger.error(optimization_error)
                optim_mode = "error"
                _df_selected = pd.DataFrame()  # Empty result

            optim_results = {
                "strategy": strat_optim.value if strat_optim else "Random",
                "k": k_optim.value,
                "iterations": iters_optim.value if iters_optim else 1000,
                "weights": normalized_wgt_optim if normalized_wgt_optim else {},
                "selected_variants": _df_selected,
                "timestamp": pd.Timestamp.now(),
                "mode": optim_mode,
                "artifact_path": (runner.output_dir / "abca4_top_variants.json") if (locals().get('runner') and optim_mode == "strand") else None,
                "error": locals().get('optimization_error') if optim_mode == "error" else None,
            }

            # Show execution path banner
            if optim_mode == "error":
                error_msg = optim_results.get("error", "Unknown optimization error")
                mo.md(f"""
                **Optimization Error: Cannot Proceed**

                {error_msg}

                **Required Actions:**
                1. Check that the Strand optimization engine is properly installed
                2. Verify MLflow configuration
                3. Enable DEMO_MODE_ENABLED in src/config.py for demo data
                4. Fix the underlying issue before proceeding
                """)
            else:
                mo.md("Strand optimization completed successfully.")

            logger.info(f"Optimization complete. Selected {len(_df_selected)} variants")

        except Exception as e:
            logger.error(f"Optimization failed: {e}")
            optim_results = None

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
        artifact = optim_results.get('artifact_path')

        if mode == "error":
            msg = f"**Optimizer mode:** `{mode}` - {error}"
        else:
            msg = f"**Optimizer mode:** `{mode}`"
            if artifact:
                msg += f"\nMLflow artifact source: `{artifact}`"
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
        _display_cols = ["rank", "chrom", "pos", "ref", "alt", "consequence", "mechanism", "suggested_assay"]
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
            "consequence", "mechanism", "suggested_assay"
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
def __(mo, pd):
    """Report configuration."""
    report_title_widget = mo.ui.text(value="ABCA4 Variant Selection Report v1", label="Title")
    _today = pd.Timestamp.now().strftime("%Y-%m-%d")
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
    report_title_widget, report_date_widget, report_notes_widget,
    _build_weight_str
):
    """Generate Markdown report."""
    if df_mapped is None or optim_results is None:
        report_md = None
    else:
        _date_str = str(report_date_widget.value) if report_date_widget.value else pd.Timestamp.now().strftime("%Y-%m-%d")

        report_md = f"""# {report_title_widget.value}

**Date:** {_date_str}

## Executive Summary

Curated selection of ABCA4 variants for functional validation.

## Selection Approach

- **Strategy:** {optim_results.get('strategy', 'Feature ranking')}
- **Panel Size (K):** {optim_results.get('k', len(df_mapped))}
- **Iterations:** {optim_results.get('iterations', 1)}

### Feature Weights

{_build_weight_str(optim_results.get('weights', {}))}

## Selected Variants ({len(df_mapped)})

| Rank | Variant | Consequence | Mechanism | Assay |
|------|---------|-------------|-----------|-------|
"""
        for _idx, _row in df_mapped.head(50).iterrows():
            _rank = _row.get("rank", _idx)
            _variant = f"{_row.get('chrom', '?')}:{_row.get('pos', '?')}:{_row.get('ref', '?')}/{_row.get('alt', '?')}"
            _consequence = _row.get("consequence", "?")
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
