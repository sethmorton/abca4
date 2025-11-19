#!/usr/bin/env python3
"""
ABCA4 Variant Selection: Optimization Dashboard

Streamlined notebook using modular utilities.
Run interactively:  marimo edit notebooks/03_optimization_dashboard.py
Run as dashboard:   marimo run notebooks/03_optimization_dashboard.py
Export to HTML:    marimo export html notebooks/03_optimization_dashboard.py -o notebooks/03_optimization_dashboard.html
"""

import marimo

__generated_with = "0.17.8"
app = marimo.App()


@app.cell
def __():
    """Setup path and import core libraries."""
    import sys
    from pathlib import Path
    # Add project root to path (works for both notebook and exported script)
    try:
        _file_path = Path(__file__).resolve()
    except NameError:
        _file_path = Path.cwd() / "notebooks" / "03_optimization_dashboard.py"
    _project_root = _file_path.parent.parent
    if str(_project_root) not in sys.path:
        sys.path.insert(0, str(_project_root))
    
    import marimo as mo
    import pandas as pd
    from src.config import logger
    from src.reward.optimization import VariantOptimizer
    from src.features.engineering.docs import PipelineDocs
    return mo, pd, logger, VariantOptimizer, PipelineDocs


@app.cell
def __(mo, PipelineDocs):
    """Notebook title and introduction."""
    mo.md(PipelineDocs.get_optimization_intro())


@app.cell
def __(mo, logger, pd):
    """Load scored variants for optimization."""
    from src.config import FEATURES_DIR

    scored_path = FEATURES_DIR / "variants_scored.parquet"
    if not scored_path.exists():
        mo.md(f"❌ Missing variants_scored.parquet at {scored_path}")
        df_scored_optim = None
    else:
        df_scored_optim = pd.read_parquet(scored_path)
        logger.info(f"Loaded {len(df_scored_optim)} scored variants for optimization")

        # Validate required columns
        required_cols = ['model_score', 'cluster_id', 'cluster_target']
        missing_cols = [col for col in required_cols if col not in df_scored_optim.columns]

        if missing_cols:
            mo.md(f"❌ Missing required columns: {missing_cols}")
            df_scored_optim = None
        else:
            mo.md(f"✅ Loaded {len(df_scored_optim)} variants with {df_scored_optim['cluster_id'].nunique()} clusters")

    return df_scored_optim


@app.cell
def __(mo, df_scored_optim):
    """Optimization parameter controls."""
    if df_scored_optim is None or df_scored_optim.empty:
        mo.md("No scored variants available for optimization.")
        k_optim = strat_optim = lambda_cov = None
    else:
        mo.md("### Optimization Parameters")

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
def __(mo, strat_optim, lambda_cov):
    """Display chosen strategy."""
    if strat_optim is None:
        mo.md("Select a strategy to proceed.")
    else:
        mo.md(f"**Strategy:** {strat_optim.value}\n\n**Coverage λ:** {lambda_cov.value if lambda_cov else 0}")


@app.cell
def __(mo, df_scored_optim, k_optim, strat_optim, lambda_cov, logger):
    """Run optimization using modular utilities."""
    if df_scored_optim is None or k_optim is None or strat_optim is None:
        mo.md("Cannot run optimization - missing data or parameters.")
        optim_results = None
    else:
        try:
            # Prepare data for optimization
            df_opt = VariantOptimizer.prepare_optimization_data(df_scored_optim)
            cluster_targets = df_opt[['cluster_id', 'cluster_target']].drop_duplicates().set_index('cluster_id')['cluster_target'].to_dict()

            # Run optimization based on strategy
            if strat_optim.value == "Greedy coverage":
                selected = VariantOptimizer.select_greedy(
                    df_opt, k_optim.value, lambda_cov.value, cluster_targets
                )
                optim_mode = "greedy_coverage"
            elif strat_optim.value == "Top score":
                # For comparison, always compute top-K
                selected = VariantOptimizer.select_greedy(df_opt, k_optim.value, 0.0, {})  # No coverage penalty
                optim_mode = "top_score"
            elif strat_optim.value == "Random":
                _selected = df_opt.sample(n=k_optim.value, random_state=42).copy()
                _selected['optimization_score'] = _selected['impact_score'] if 'impact_score' in _selected.columns else _selected['model_score']
                _selected['selected'] = True
                _selected['rank'] = range(1, len(_selected) + 1)
                optim_mode = "random"
            else:
                mo.md(f"❌ Unsupported strategy: {strat_optim.value}")
                optim_results = None

            # Always compute top-K for comparison
            top_k_by_score = VariantOptimizer.select_greedy(df_opt, k_optim.value, 0.0, {})  # No coverage penalty

            # Handle different selection strategies
            if strat_optim.value == "Random":
                final_selected = _selected
            else:
                final_selected = selected

            # Compute comparison metrics
            _comparison = VariantOptimizer.compute_comparison_metrics(final_selected, top_k_by_score)

            optim_results = {
                "strategy": strat_optim.value,
                "k": k_optim.value,
                "coverage_lambda": lambda_cov.value if lambda_cov else 0.0,
                "selected_variants": final_selected,
                "top_k_by_score": top_k_by_score,
                "comparison": _comparison,
                "mode": optim_mode
            }

            logger.info(f"Optimization complete: {optim_mode}, selected {len(final_selected)} variants")

        except Exception as e:
            logger.error(f"Optimization failed: {e}")
            mo.md(f"❌ Optimization failed: {e}")
            optim_results = None

    return optim_results


@app.cell
def __(mo, optim_results):
    """Display optimization results and comparison."""
    if optim_results is None:
        mo.md("No optimization results to display.")
    else:
        # Use modular template for results display
        mo.md(PipelineDocs.get_optimization_results_template(optim_results))


@app.cell
def __(mo, optim_results):
    """Display optimization score distribution."""
    if optim_results is None or 'selected_variants' not in optim_results:
        mo.md("No results to plot.")
    else:
        try:
            import plotly.graph_objects as go

            _selected = optim_results['selected_variants']
            fig = go.Figure()

            fig.add_trace(go.Bar(
                x=_selected["rank"],
                y=_selected["optimization_score"],
                name="Selected Variants",
                marker_color="lightseagreen"
            ))

            fig.update_layout(
                title="Selected Variants by Optimization Score",
                xaxis_title="Rank",
                yaxis_title="Score",
                showlegend=False,
                template="plotly_white"
            )

            mo.ui.plotly(fig)

        except Exception as e:
            mo.md(f"❌ Plot error: {e}")


@app.cell
def __(mo, optim_results):
    """Export selected variants."""
    if optim_results is None or 'selected_variants' not in optim_results:
        mo.md("No results to export.")
    else:
        try:
            from src.config import REPORTS_DIR

            _selected = optim_results['selected_variants']

            # Export CSV
            csv_path = REPORTS_DIR / "variants_selected.csv"
            _selected.to_csv(csv_path, index=False)

            # Export JSON
            json_path = REPORTS_DIR / "variants_selected.json"
            _selected.to_json(json_path, orient='records', indent=2)

            mo.md(f"✅ **Exported {len(_selected)} variants to:**\n- {csv_path}\n- {json_path}")

        except Exception as e:
            mo.md(f"❌ Export failed: {e}")


@app.cell
def __(mo, optim_results):
    """Summary and next steps."""
    if optim_results is None:
        mo.md("Optimization incomplete.")
    else:
        _comparison = optim_results.get('comparison', {})
        impact_diff = _comparison.get('impact_diff_pct', 0)

        status = "✅ **SUCCESS**" if abs(impact_diff) <= 5 else "⚠️ **REVIEW NEEDED**"

        mo.md(f"""
## Optimization Summary

{status}

**Panel Size:** {optim_results['k']} variants
**Strategy:** {optim_results['strategy']}
**Coverage λ:** {optim_results['coverage_lambda']}

**Impact vs Baseline:** {impact_diff:+.1f}%

### Next Steps

1. **Review selected variants** in exported CSV/JSON files
2. **Generate assay drafts** using `invoke reporting.drafts` or `invoke run-pipeline`
3. **Proceed to CRO planning** for experimental design (includes LLM-generated assay drafts)
4. **Adjust parameters** if optimization results need improvement

**Assay Drafts Output:** `data_processed/reports/assay_drafts/`
""")


# Run the app
if __name__ == "__main__":
    app.run()
