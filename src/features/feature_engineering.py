#!/usr/bin/env python3
"""
Feature engineering utilities for ABCA4 pipeline.

Moved from notebooks/02_feature_engineering.py to keep notebooks thin.
"""

import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import importlib

import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score
import plotly.graph_objects as go

from src.config import logger


class FeatureEngineeringUtils:
    """Utility functions for feature engineering pipeline."""

    @staticmethod
    def load_or_compute_required_features(required_slices: Dict, campaign_root: Path) -> Tuple[List[Tuple[pd.DataFrame, str]], Dict]:
        """
        Load or compute required feature slices with strict validation.

        Raises SystemExit with actionable error message if any required slice fails.
        Returns (feature_sources, status_summary) on success.
        """
        import json
        from pathlib import Path

        feature_sources = []
        status_summary = {"timestamp": pd.Timestamp.now().isoformat(), "slices": {}}

        for slice_name, config in required_slices.items():
            slice_path = campaign_root / config['path_relative']
            description = config['description']
            required_cols = config['required_columns']

            logger.info(f"Loading required feature slice: {slice_name} ({description})")

            try:
                # Check if file exists
                if not slice_path.exists():
                    # Try to compute it
                    logger.info(f"Feature slice {slice_name} missing, attempting computation...")
                    try:
                        module_name, class_name = config['computer_class'].rsplit('.', 1)
                        module = importlib.import_module(module_name)
                        computer_class = getattr(module, class_name)
                        computer = computer_class()

                        if computer.run():
                            logger.info(f"Successfully computed {slice_name} features")
                        else:
                            raise RuntimeError(f"Feature computation failed for {slice_name}")

                    except Exception as compute_error:
                        error_msg = f"""
❌ REQUIRED FEATURE SLICE MISSING: {slice_name}

Description: {description}
Expected file: {slice_path}
Computation error: {compute_error}

ACTION REQUIRED:
Feature slice appears corrupted. Re-run computation for {slice_name}.
"""
                        raise SystemExit(error_msg)

                # Load and validate
                df_slice = pd.read_parquet(slice_path)

                # Check required columns
                missing_cols = [col for col in required_cols if col not in df_slice.columns]
                if missing_cols:
                    error_msg = f"""
❌ REQUIRED COLUMNS MISSING: {slice_name}

Description: {description}
File: {slice_path}
Missing columns: {missing_cols}
Available columns: {list(df_slice.columns)}

ACTION REQUIRED:
Feature slice appears corrupted. Re-run computation for {slice_name}.
"""
                    raise SystemExit(error_msg)

                # Add to feature sources
                feature_sources.append((df_slice, config.get('join_key', 'variant_id')))
                status_summary["slices"][slice_name] = {
                    "status": "loaded",
                    "path": str(slice_path),
                    "rows": len(df_slice),
                    "columns": len(df_slice.columns)
                }

                logger.info(f"✅ Loaded {slice_name}: {len(df_slice)} rows, {len(df_slice.columns)} columns")

            except Exception as e:
                status_summary["slices"][slice_name] = {
                    "status": "error",
                    "error": str(e)
                }
                raise

        return feature_sources, status_summary

    @staticmethod
    def compute_impact_scores(df_scored_step3: pd.DataFrame, scoring_mode: str,
                            alpha_wgt: float = 0.4, splice_wgt: float = 0.3,
                            cons_wgt: float = 0.15, lof_wgt: float = 0.15,
                            features_dir: Optional[Path] = None) -> pd.DataFrame:
        """
        Compute impact scores using hand-mix or logistic regression.
        Moved from notebook to keep it thin.
        """
        df_impact = df_scored_step3.copy()
        # Preserve all attributes from input dataframe (save snapshot – merges below will drop attrs)
        _original_attrs = dict(df_scored_step3.attrs)
        scoring_error = None

        # Merge domain annotations if available
        if features_dir:
            domains_path = features_dir / "variants_with_domains.parquet"
            if domains_path.exists():
                try:
                    df_domains = pd.read_parquet(domains_path)
                    if not df_domains.empty and 'domain' in df_domains.columns and 'variant_id' in df_domains.columns:
                        df_impact = df_impact.merge(
                            df_domains[['variant_id', 'domain']],
                            on='variant_id',
                            how='left',
                            suffixes=('', '_domain')
                        )
                        logger.info(f"Merged domain annotations for impact scoring")
                except Exception as e:
                    logger.warning(f"Could not load domain annotations for impact scoring: {e}")

        if scoring_mode == "hand-mix":
            # Normalize scores to [0, 1]
            def _normalize_score(col):
                if col not in df_impact.columns or df_impact[col].notna().sum() == 0:
                    return np.zeros(len(df_impact))
                s = df_impact[col].fillna(0.0)
                s_min, s_max = s.min(), s.max()
                if s_max > s_min:
                    return (s - s_min) / (s_max - s_min)
                else:
                    return np.zeros(len(s))

            _alpha_norm = _normalize_score("alphamissense_score")
            _splice_norm = _normalize_score("spliceai_max_score")
            _cons_norm = _normalize_score("phylop_score")
            _lof_norm = df_impact["lof_prior"].fillna(0.0)

            _total_wgt = (alpha_wgt + splice_wgt + cons_wgt + lof_wgt)
            if _total_wgt == 0:
                _total_wgt = 1.0

            df_impact["model_score"] = (
                (alpha_wgt * _alpha_norm +
                 splice_wgt * _splice_norm +
                 cons_wgt * _cons_norm +
                 lof_wgt * _lof_norm) / _total_wgt
            )

            logger.info(f"Computed hand-mix impact scores.")

        # Plan-conformant feature vector
        # cons_v: scaled conservation (0-1)
        cons_col = next((c for c in ['phylop_score', 'phyloP100way', 'phyloP100way_z'] if c in df_impact.columns), None)
        if cons_col:
            _cons = df_impact[cons_col].fillna(0.0)
            _cons_range = _cons.max() - _cons.min()
            df_impact['cons_scaled'] = (
                (_cons - _cons.min()) / _cons_range if _cons_range > 0 else np.zeros(len(df_impact))
            ).clip(0, 1)
        else:
            df_impact['cons_scaled'] = 0.0

        # af_penalty: monotonic positive penalty for common variants
        af_col = next((c for c in [
            'gnomad_af_exome', 'gnomad_af_genome', 'gnomad_max_af', 'faf95_max'
        ] if c in df_impact.columns), None)

        if af_col:
            af_raw = df_impact[af_col].fillna(0.0)
            af_clipped = af_raw.clip(lower=1e-6, upper=0.05)
            af_penalty = np.clip(-np.log10(af_clipped) / 6, 0, 1)
            df_impact['af_penalty'] = af_penalty
        else:
            df_impact['af_penalty'] = 0.0

        # domain_flag: 1 if variant sits in a defined domain
        if 'domain' in df_impact.columns:
            df_impact['domain_flag'] = df_impact['domain'].notna().astype(int)
        elif 'in_domain' in df_impact.columns:
            df_impact['domain_flag'] = df_impact['in_domain'].fillna(False).astype(int)
        elif 'domain_label' in df_impact.columns:
            df_impact['domain_flag'] = df_impact['domain_label'].notna().astype(int)
        else:
            df_impact['domain_flag'] = 0

        # splice_prox_flag: consequence mentions splice or close to exon boundary
        def _is_splice(row):
            cons = str(row.get('vep_consequence', '')).lower()
            dist = row.get('intron_distance', np.nan)
            return ('splice' in cons) or (pd.notna(dist) and dist <= 10)

        if 'vep_consequence' in df_impact.columns or 'intron_distance' in df_impact.columns:
            df_impact['splice_prox_flag'] = df_impact.apply(_is_splice, axis=1).astype(int)
        else:
            df_impact['splice_prox_flag'] = 0

        # impact_score: hand-mixed linear score aligned to plan
        df_impact['impact_score'] = (
            0.6 * df_impact.get('model_score', 0) +
            0.2 * df_impact['cons_scaled'] +
            0.1 * df_impact['domain_flag'] +
            0.1 * df_impact['splice_prox_flag'] -
            0.3 * df_impact['af_penalty']
        )
        df_impact['impact_score'] = df_impact['impact_score'].clip(0.0, 1.0)

        if scoring_error:
            df_impact.attrs['scoring_error'] = scoring_error
        else:
            df_impact.attrs['scoring_error'] = None
        # Re-apply original attrs at the very end to restore validation metadata
        try:
            df_impact.attrs.update(_original_attrs)
        except Exception:
            # Best-effort; never fail scoring due to attrs
            pass

        return df_impact

    @staticmethod
    def create_qc_plots(df_impact: pd.DataFrame) -> Tuple[Optional[go.Figure], Optional[go.Figure]]:
        """
        Create QC plots: AUROC for P/LP vs B/LB and density overlays.
        Returns (main_plot, domain_plot) or (None, None) if insufficient data.
        """
        if "impact_score" not in df_impact.columns or df_impact.empty:
            return None, None

        try:
            # Prepare clinical significance labels
            clinsig_col = next((c for c in df_impact.columns if 'clinical_significance' in c.lower()), None)
            if not clinsig_col:
                return None, None

            def _is_pathogenic(sig):
                if not isinstance(sig, str):
                    return False
                sig_lower = str(sig).lower()
                return 'pathogenic' in sig_lower and 'benign' not in sig_lower

            def _is_benign(sig):
                if not isinstance(sig, str):
                    return False
                sig_lower = str(sig).lower()
                return 'benign' in sig_lower and 'pathogenic' not in sig_lower

            df_qc = df_impact.copy()
            df_qc['is_pathogenic'] = df_qc[clinsig_col].apply(_is_pathogenic)
            df_qc['is_benign'] = df_qc[clinsig_col].apply(_is_benign)

            pathogenic_mask = df_qc['is_pathogenic']
            benign_mask = df_qc['is_benign']

            n_pathogenic = pathogenic_mask.sum()
            n_benign = benign_mask.sum()

            if n_pathogenic < 10 or n_benign < 10:
                return None, None

            # Compute AUROC
            y_true = pd.concat([pd.Series([1] * n_pathogenic), pd.Series([0] * n_benign)])
            y_score = pd.concat([df_qc.loc[pathogenic_mask, 'impact_score'], df_qc.loc[benign_mask, 'impact_score']])

            try:
                auroc = roc_auc_score(y_true, y_score)
                logger.info(f"AUROC for P/LP vs B/LB: {auroc:.3f} ({n_pathogenic} P/LP, {n_benign} B/LB variants)")
            except Exception as e:
                logger.warning(f"AUROC calculation failed: {e}")
                return None, None

            # Create density overlay plot
            fig_qc = go.Figure()

            # P/LP distribution
            pathogenic_scores = df_qc.loc[pathogenic_mask, 'impact_score'].dropna()
            if not pathogenic_scores.empty:
                fig_qc.add_trace(go.Histogram(
                    x=pathogenic_scores,
                    nbinsx=30,
                    name="P/LP",
                    marker_color="rgba(239, 85, 59, 0.7)",
                    opacity=0.7
                ))

            # B/LB distribution
            benign_scores = df_qc.loc[benign_mask, 'impact_score'].dropna()
            if not benign_scores.empty:
                fig_qc.add_trace(go.Histogram(
                    x=benign_scores,
                    nbinsx=30,
                    name="B/LB",
                    marker_color="rgba(99, 110, 250, 0.7)",
                    opacity=0.7
                ))

            fig_qc.update_layout(
                title="Impact Score Distribution: P/LP vs B/LB",
                xaxis_title="Impact Score",
                yaxis_title="Frequency",
                hovermode="x unified",
                barmode="overlay",
                template="plotly_white"
            )

            # Domain-based analysis (hotspots vs benign regions)
            fig_domain = None
            if 'domain' in df_qc.columns:
                # Define hotspots (NBD1, NBD2, TMD1, TMD2)
                hotspot_domains = ['NBD1', 'NBD2', 'TMD1', 'TMD2']
                df_qc['is_hotspot'] = df_qc['domain'].isin(hotspot_domains)

                hotspot_pathogenic = df_qc[pathogenic_mask & df_qc['is_hotspot']]
                benign_hotspot = df_qc[benign_mask & df_qc['is_hotspot']]
                benign_nonhotspot = df_qc[benign_mask & ~df_qc['is_hotspot']]

                fig_domain = go.Figure()

                if not hotspot_pathogenic.empty:
                    fig_domain.add_trace(go.Histogram(
                        x=hotspot_pathogenic['impact_score'].dropna(),
                        nbinsx=20,
                        name="P/LP in hotspots",
                        marker_color="rgba(239, 85, 59, 0.8)",
                        opacity=0.8
                    ))

                if not benign_hotspot.empty:
                    fig_domain.add_trace(go.Histogram(
                        x=benign_hotspot['impact_score'].dropna(),
                        nbinsx=20,
                        name="B/LB in hotspots",
                        marker_color="rgba(255, 165, 0, 0.7)",
                        opacity=0.7
                    ))

                if not benign_nonhotspot.empty:
                    fig_domain.add_trace(go.Histogram(
                        x=benign_nonhotspot['impact_score'].dropna(),
                        nbinsx=20,
                        name="B/LB outside hotspots",
                        marker_color="rgba(99, 110, 250, 0.6)",
                        opacity=0.6
                    ))

                fig_domain.update_layout(
                    title="Domain-based QC: Hotspots vs Non-hotspots",
                    xaxis_title="Impact Score",
                    yaxis_title="Frequency",
                    hovermode="x unified",
                    barmode="overlay",
                    template="plotly_white"
                )

            return fig_qc, fig_domain

        except Exception as e:
            logger.error(f"QC analysis error: {e}")
            return None, None

    @staticmethod
    def merge_features_with_validation(df_base: pd.DataFrame, feature_sources: List[Tuple[pd.DataFrame, str]],
                                     features_dir: Path) -> pd.DataFrame:
        """
        Merge feature sources with validation and quality checks.
        Moved from notebook to keep it thin.
        """
        df_result = df_base.copy()

        # Apply conservation column standardization
        for i, (df_features, join_key) in enumerate(feature_sources):
            if any("phylop" in col.lower() or "phastcons" in col.lower() for col in df_features.columns):
                from src.config import standardize_conservation_columns
                feature_sources[i] = (standardize_conservation_columns(df_features), join_key)

        # Join all feature sources (deduplicate first to avoid cartesian product)
        initial_row_count = len(df_result)
        merge_validation = []

        for df_features, join_key in feature_sources:
            if not df_features.empty and join_key in df_features.columns:
                # Check for duplicates before merge
                original_len = len(df_features)
                df_features_dedup = df_features.drop_duplicates(subset=[join_key], keep='first')
                duplicate_count = original_len - len(df_features_dedup)

                # Perform merge
                pre_merge_len = len(df_result)
                df_result = df_result.merge(df_features_dedup, on=join_key, how='left', suffixes=('', '_feature'))
                post_merge_len = len(df_result)

                # Validate merge
                if post_merge_len != initial_row_count:
                    raise RuntimeError(f"❌ Merge validation failed: row count changed from {initial_row_count} to {post_merge_len}")

                merge_validation.append({
                    'feature_source': join_key,
                    'duplicates_removed': int(duplicate_count),
                    'rows_before_merge': int(pre_merge_len),
                    'rows_after_merge': int(post_merge_len)
                })

                logger.info(f"Merged {join_key}: removed {duplicate_count} duplicates, {post_merge_len} total rows")

        # Store merge validation results
        df_result.attrs['merge_validation'] = merge_validation

        # Validate non-null percentages for key columns
        key_columns = {
            'alphamissense_score': 'AlphaMissense',
            'spliceai_max_score': 'SpliceAI',
            'phyloP100way': 'conservation',
            'gnomad_exome_af': 'AF (exome)',
            'gnomad_genome_af': 'AF (genome)',
            'gnomad_max_af': 'AF (max)',
            'domain': 'domain'
        }

        null_rates = {}
        total_rows = len(df_result)

        for col, description in key_columns.items():
            if col in df_result.columns:
                non_null_count = df_result[col].notna().sum()
                null_rate = (total_rows - non_null_count) / total_rows
                null_rates[col] = {
                    'description': description,
                    'non_null_count': int(non_null_count),
                    'null_rate': float(null_rate),
                    'non_null_pct': float((non_null_count / total_rows) * 100)
                }

                # Adjust thresholds for known limitations
                min_pct = 40 if col == 'alphamissense_score' else 95

                if (non_null_count / total_rows) * 100 < min_pct:
                    raise SystemExit(f"❌ FEATURE QUALITY CHECK FAILED: {description} has {(non_null_count/total_rows*100):.1f}% non-null values "
                                   f"({non_null_count}/{total_rows} non-null). Required: ≥{min_pct}% non-null.")

                logger.info(f"✅ {description}: {(non_null_count/total_rows*100):.1f}% non-null ({non_null_count}/{total_rows})")

        df_result.attrs['null_rates'] = null_rates

        # Store column origin metadata for transparency
        df_result.attrs['feature_sources_metadata'] = {
            'missense': {
                'file': str(features_dir / "missense_features.parquet"),
                'columns': ['alphamissense_score', 'alphamissense_pathogenic_category', 'esm_pld_score'],
                'description': 'AlphaMissense and ESM protein language model predictions'
            },
            'splice': {
                'file': str(features_dir / "splice_features.parquet"),
                'columns': ['spliceai_max_score', 'spliceai_ds_ag', 'spliceai_ds_al', 'spliceai_ds_dg', 'spliceai_ds_dl'],
                'description': 'SpliceAI splice site disruption scores'
            },
            'conservation': {
                'file': str(features_dir / "conservation_features.parquet"),
                'columns': ['phyloP100way', 'phastCons100way'],
                'description': 'Phylogenetic conservation (phyloP100way and phastCons100way)'
            },
            'regulatory': {
                'file': str(features_dir / "regulatory_features.parquet"),
                'columns': ['domain', 'domain_name', 'gnomad_af_exome', 'gnomad_ac_exome'],
                'description': 'Protein domains and gnomAD population allele frequencies'
            }
        }

        # Add LoF prior based on consequence
        lof_mapping = {
            "frameshift_variant": 0.95,
            "stop_gained": 0.95,
            "splice_acceptor_variant": 0.95,
            "splice_donor_variant": 0.95,
            "missense_variant": 0.1,
            "synonymous_variant": 0.01,
        }

        if "lof_prior" not in df_result.columns:
            consequence_col = next((c for c in df_result.columns if 'consequence' in c.lower()), 'vep_consequence')
            df_result["lof_prior"] = df_result.get(consequence_col, "missense_variant").apply(
                lambda x: lof_mapping.get(str(x).lower() if pd.notna(x) else "missense_variant", 0.1)
            )

        logger.info(f"Loaded all features. {len(df_result)} variants with {len(df_result.columns)} total columns.")

        return df_result
