#!/usr/bin/env python3
"""
ABCA4 Variant Analysis: Feature Engineering

Streamlined notebook using modular utilities.
Run interactively:  marimo edit notebooks/02_feature_engineering.py
Run as dashboard:   marimo run notebooks/02_feature_engineering.py
Export to HTML:    marimo export html notebooks/02_feature_engineering.py -o notebooks/02_feature_engineering.html
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
        _file_path = Path.cwd() / "notebooks" / "02_feature_engineering.py"
    _project_root = _file_path.parent.parent
    if str(_project_root) not in sys.path:
        sys.path.insert(0, str(_project_root))
    
    import marimo as mo
    import pandas as pd
    import numpy as np
    from src.config import logger, CAMPAIGN_ROOT, FEATURES_DIR
    from src.features.feature_engineering import FeatureEngineeringUtils
    from src.features.docs import PipelineDocs
    return mo, pd, np, Path, logger, CAMPAIGN_ROOT, FEATURES_DIR, FeatureEngineeringUtils, PipelineDocs


@app.cell
def __(mo, PipelineDocs):
    """Notebook title and introduction."""
    mo.md(PipelineDocs.get_feature_engineering_intro())


@app.cell
def __(mo, logger, pd):
    """Load and validate annotated variants."""
    from src.config import DEMO_MODE_ENABLED, ANNOTATIONS_DIR

    # Load annotated variants
    annotated_path = ANNOTATIONS_DIR / "abca4_vus_annotated.parquet"
    if not annotated_path.exists():
        mo.md(f"❌ Missing annotated variants at {annotated_path}")
        df_annotated = None
    else:
        df_annotated = pd.read_parquet(annotated_path)
        logger.info(f"Loaded {len(df_annotated)} annotated variants")

        # Smoke test
        if len(df_annotated) == 0:
            mo.md("❌ No variants loaded!")
        else:
            mo.md(f"✅ Loaded {len(df_annotated)} variants successfully")


@app.cell
def __(mo, df_annotated, PipelineDocs):
    """Data quality audit display."""
    if df_annotated is None or df_annotated.empty:
        mo.md("No data to audit.")
    else:
        mo.md(PipelineDocs.get_data_quality_audit())

    # Basic stats
    total_vars = len(df_annotated)
    mo.md(f"""
**Dataset Overview:**
- Total variants: {total_vars:,}
- Columns: {len(df_annotated.columns)}
- Missing protein_change: {(df_annotated['protein_change'].isna().sum() / total_vars * 100):.1f}%
""")


@app.cell
def __(mo, df_annotated, FEATURES_DIR, logger, CAMPAIGN_ROOT, FeatureEngineeringUtils, PipelineDocs):
    """Load and merge all features using modular utilities."""
    if df_annotated is None or df_annotated.empty:
        mo.md("Cannot proceed without annotated variants.")
        df_scored_step3 = None
    else:
        # Define required feature slices
        REQUIRED_FEATURE_SLICES = {
            'missense': {
                'path_relative': "data_processed/features/missense_features.parquet",
                'computer_class': 'src.features.missense.MissenseFeatureComputer',
                'required_columns': ['variant_id', 'alphamissense_score'],
                'description': 'AlphaMissense + ESM features'
            },
            'splice': {
                'path_relative': "data_processed/features/splice_features.parquet",
                'computer_class': 'src.features.splice.SpliceFeatureComputer',
                'required_columns': ['variant_id', 'spliceai_max_score'],
                'description': 'SpliceAI features'
            },
            'conservation': {
                'path_relative': "data_processed/features/conservation_features.parquet",
                'computer_class': 'src.features.conservation.ConservationFeatureComputer',
                'required_columns': ['variant_id', 'phyloP100way'],
                'description': 'phyloP/phastCons conservation scores'
            },
            'regulatory': {
                'path_relative': "data_processed/features/regulatory_features.parquet",
                'computer_class': 'src.features.regulatory.RegulatoryFeatureComputer',
                'required_columns': ['variant_id', 'gnomad_exome_af', 'gnomad_genome_af', 'gnomad_max_af'],
                'description': 'domains + gnomAD regulatory features'
            }
        }

        try:
            # Use modular utility for feature loading and merging
            feature_sources, status_summary = FeatureEngineeringUtils.load_or_compute_required_features(
                REQUIRED_FEATURE_SLICES, CAMPAIGN_ROOT
            )

            # Use modular utility for merging with validation
            df_scored_step3 = FeatureEngineeringUtils.merge_features_with_validation(
                df_annotated, feature_sources, FEATURES_DIR
            )

            # Display status
            mo.md(PipelineDocs.get_feature_status_template(status_summary.get('slices', {})))

            # Display merge validation
            merge_validation = df_scored_step3.attrs.get('merge_validation', [])
            if merge_validation:
                mo.md("**Merge Validation:**")
                mo.md(PipelineDocs.get_merge_validation_template(merge_validation))

            # Display null rates
            null_rates = df_scored_step3.attrs.get('null_rates', {})
            if null_rates:
                mo.md("**Data Quality (Non-Null Rates):**")
                mo.md(PipelineDocs.get_null_rates_template(null_rates))

                min_pct = min(data.get('non_null_pct', 0) for data in null_rates.values())
                if min_pct >= 40:  # Adjusted threshold for AlphaMissense
                    mo.md(f"✅ **All features meet quality thresholds** (minimum: {min_pct:.1f}%)")
                else:
                    mo.md(f"❌ **Quality check failed:** Some features below threshold")

            logger.info(f"Feature engineering complete: {len(df_scored_step3)} variants")

        except Exception as e:
            mo.md(f"❌ Feature engineering failed: {e}")
            logger.error(f"Feature engineering error: {e}")
            df_scored_step3 = None


@app.cell
def __(mo, PipelineDocs):
    """Impact scoring section introduction."""
    mo.md(PipelineDocs.get_impact_scoring_intro())


@app.cell
def __(mo, df_scored_step3, FEATURES_DIR, logger, FeatureEngineeringUtils):
    """Compute impact scores using modular utilities."""
    if df_scored_step3 is None or df_scored_step3.empty:
        mo.md("Cannot compute impact scores without features.")
        df_impact = None
    else:
        try:
            # Use modular utility for impact scoring
            df_impact = FeatureEngineeringUtils.compute_impact_scores(
                df_scored_step3,
                scoring_mode="hand-mix",
                features_dir=FEATURES_DIR
            )

            # Check for errors
            scoring_error = df_impact.attrs.get('scoring_error')
            if scoring_error:
                mo.md(f"⚠️ Scoring warning: {scoring_error}")
            else:
                mo.md("✅ Impact scores computed successfully")

            # Display basic stats
            impact_range = df_impact['impact_score'].agg(['min', 'max'])
            mo.md(f"""
**Impact Score Statistics:**
- Range: {impact_range['min']:.3f} - {impact_range['max']:.3f}
- Mean: {df_impact['impact_score'].mean():.3f}
- Domain flag distribution: {df_impact['domain_flag'].value_counts().to_dict()}
""")

            # Save intermediate results
            features_raw_path = FEATURES_DIR / "variants_features_raw.parquet"
            df_impact.to_parquet(features_raw_path)
            logger.info(f"Saved intermediate features to {features_raw_path}")

        except Exception as e:
            mo.md(f"❌ Impact scoring failed: {e}")
            logger.error(f"Impact scoring error: {e}")
            df_impact = None


@app.cell
def __(mo, df_impact, FeatureEngineeringUtils):
    """QC plots for impact score validation."""
    if df_impact is None or df_impact.empty:
        mo.md("Cannot create QC plots without impact scores.")
    else:
        # Use modular utility for QC plots
        fig_qc, fig_domain = FeatureEngineeringUtils.create_qc_plots(df_impact)

        if fig_qc:
            mo.ui.plotly(fig_qc)
        else:
            mo.md("❌ Could not generate P/LP vs B/LB QC plot")

        if fig_domain:
            mo.ui.plotly(fig_domain)
        else:
            mo.md("❌ Could not generate domain-based QC plot")


@app.cell
def __(mo, PipelineDocs):
    """Clustering section introduction."""
    mo.md(PipelineDocs.get_clustering_intro())


@app.cell
def __(mo, df_impact, logger, FEATURES_DIR, PipelineDocs):
    """Apply clustering using modular utilities."""
    df_clusters = None
    if df_impact is None or df_impact.empty:
        mo.md("Cannot perform clustering without impact scores.")
    else:
        try:
            from src.features.add_clustering import ClusteringProcessor

            # Apply domain boost and clustering
            processor = ClusteringProcessor()

            # Apply domain boost
            df_boosted = processor.apply_domain_boost(df_impact)

            # Apply hotspot overrides
            df_overridden = processor.apply_hotspot_overrides(df_boosted)

            # Assign clusters using domain-based approach
            df_clustered = processor.assign_clusters_domain(df_overridden)

            # Compute biology-aware targets
            df_with_targets = processor.compute_cluster_targets_biology_aware(df_clustered)

            # Compute coverage
            df_clusters = processor.compute_cluster_coverage(df_with_targets)

            # Preserve metadata attributes from original dataframe
            df_clusters.attrs.update(df_impact.attrs)

            # Display summary
            mo.md(PipelineDocs.get_cluster_summary_template(df_clusters))

            # Save clustered variants directly to preserve attributes
            scored_path = FEATURES_DIR / "variants_scored.parquet"
            csv_path = FEATURES_DIR / "variants_scored.csv"

            try:
                df_clusters.to_parquet(scored_path, index=False)
                df_clusters.to_csv(csv_path, index=False)
                logger.info(f"Saved {len(df_clusters)} clustered variants with {len(df_clusters.attrs)} attributes to {scored_path}")

                # Log cluster summary
                if "cluster_id" in df_clusters.columns:
                    cluster_counts = df_clusters["cluster_id"].value_counts()
                    logger.info(f"Cluster distribution:\\n{cluster_counts}")

            except Exception as e:
                logger.error(f"Failed to save clustered variants: {e}")
                mo.md(f"❌ Failed to save results: {e}")

            logger.info(f"Clustering complete: {df_clusters['cluster_id'].nunique()} clusters")

        except Exception as e:
            mo.md(f"❌ Clustering failed: {e}")
            logger.error(f"Clustering error: {e}")
            df_clusters = None


@app.cell
def __(mo, df_clusters, FEATURES_DIR):
    """Save final clustered and scored variants."""
    if df_clusters is None:
        mo.md("Cannot save results without clustered data.")
    else:
        try:
            # Save to variants_scored.parquet
            _scored_path = FEATURES_DIR / "variants_scored.parquet"
            df_clusters.to_parquet(_scored_path)

            # Save CSV version
            _csv_path = FEATURES_DIR / "variants_scored.csv"
            df_clusters.to_csv(_csv_path, index=False)

            mo.md(f"""
✅ **Pipeline Complete!**

**Outputs saved:**
- `variants_scored.parquet/csv`: {len(df_clusters)} clustered variants
- `variants_features_raw.parquet`: Intermediate features with impact scores

**Ready for optimization in next notebook.**
""")

        except Exception as e:
            mo.md(f"❌ Failed to save results: {e}")


# Run the app
if __name__ == "__main__":
    app.run()
