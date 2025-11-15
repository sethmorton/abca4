#!/usr/bin/env python3
"""
ABCA4 Campaign – Feature Engineering & Scoring

Run interactively:  marimo edit notebooks/02_feature_engineering.py
Run as dashboard:   marimo run notebooks/02_feature_engineering.py
Run as script:      python notebooks/02_feature_engineering.py
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
    import logging
    from typing import Optional, Dict, List, Tuple
    import sys
    import importlib

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger(__name__)

    return mo, pd, np, Path, logging, logger, Optional, Dict, List, Tuple, sys, importlib


@app.cell
def __(mo, Path, logger, pd):
    """Define campaign paths and load annotated dataset."""
    CAMPAIGN_ROOT = Path(__file__).resolve().parents[0]
    ANNOTATIONS_DIR = CAMPAIGN_ROOT / "data_processed" / "annotations"
    FEATURES_DIR = CAMPAIGN_ROOT / "data_processed" / "features"
    FEATURES_DIR.mkdir(parents=True, exist_ok=True)

    annotated_path = ANNOTATIONS_DIR / "abca4_vus_annotated.parquet"
    if not annotated_path.exists():
        mo.md(f"⚠️ Missing annotated variants at {annotated_path}. Run 01_data_exploration.py first.")
        df_annotated = pd.DataFrame()
    else:
        df_annotated = pd.read_parquet(annotated_path)
        logger.info(f"Loaded {len(df_annotated)} annotated variants")

    return CAMPAIGN_ROOT, ANNOTATIONS_DIR, FEATURES_DIR, annotated_path, df_annotated


@app.cell
def __(mo):
    """
    ## Step 3: Main Model Scoring

    Load AlphaMissense scores, SpliceAI predictions, and construct LoF priors.
    """
    mo.md(__doc__)


@app.cell
def __(
    pd, np, logger,
    df_annotated, CAMPAIGN_ROOT, ANNOTATIONS_DIR, FEATURES_DIR
):
    """
    Load AlphaMissense, SpliceAI, conservation, and regulatory features
    by calling the authoritative feature computation scripts.
    """
    sys.path.insert(0, str(CAMPAIGN_ROOT))
    
    df_scored_step3 = df_annotated.copy()

    feature_logs: List[Dict[str, object]] = []

    # Load features from each specialist script
    feature_sources = []
    
    # 1. Load or compute missense features (AlphaMissense + ESM)
    try:
        from src.features.missense import MissenseFeatureComputer
        missense_path = FEATURES_DIR / "missense_features.parquet"
        status = "cache"
        if missense_path.exists():
            logger.info(f"Loading cached missense features from {missense_path}")
            df_missense = pd.read_parquet(missense_path)
        else:
            logger.info("Computing missense features via AlphaMissense...")
            computer = MissenseFeatureComputer(input_dir=ANNOTATIONS_DIR, output_dir=FEATURES_DIR)
            if computer.run():
                df_missense = pd.read_parquet(missense_path)
                status = "computed"
            else:
                logger.warning("Missense feature computation failed; using fallback")
                df_missense = pd.DataFrame()
                status = "fallback"
        feature_sources.append((df_missense, 'variant_id'))
        feature_logs.append({"feature": "missense", "status": status, "rows": len(df_missense)})
    except Exception as e:
        logger.error(f"Failed to load missense features: {e}")
    
    # 2. Load or compute splice features (SpliceAI)
    try:
        from src.features.splice import SpliceFeatureComputer
        splice_path = FEATURES_DIR / "splice_features.parquet"
        status = "cache"
        if splice_path.exists():
            logger.info(f"Loading cached splice features from {splice_path}")
            df_splice = pd.read_parquet(splice_path)
        else:
            logger.info("Computing splice features via SpliceAI...")
            computer = SpliceFeatureComputer(input_dir=ANNOTATIONS_DIR, output_dir=FEATURES_DIR)
            if computer.run():
                df_splice = pd.read_parquet(splice_path)
                status = "computed"
            else:
                logger.warning("Splice feature computation failed; using fallback")
                df_splice = pd.DataFrame()
                status = "fallback"
        feature_sources.append((df_splice, 'variant_id'))
        feature_logs.append({"feature": "splice", "status": status, "rows": len(df_splice)})
    except Exception as e:
        logger.error(f"Failed to load splice features: {e}")
    
    # 3. Load or compute conservation features (phyloP/phastCons)
    try:
        from src.features.conservation import ConservationFeatureComputer
        cons_path = FEATURES_DIR / "conservation_features.parquet"
        status = "cache"
        if cons_path.exists():
            logger.info(f"Loading cached conservation features from {cons_path}")
            df_cons = pd.read_parquet(cons_path)
        else:
            logger.info("Computing conservation features via UCSC...")
            computer = ConservationFeatureComputer(annotations_dir=ANNOTATIONS_DIR, output_dir=FEATURES_DIR)
            if computer.run():
                df_cons = pd.read_parquet(cons_path)
                status = "computed"
            else:
                logger.warning("Conservation feature computation failed; using fallback")
                df_cons = pd.DataFrame()
                status = "fallback"
        feature_sources.append((df_cons, 'variant_id'))
        feature_logs.append({"feature": "conservation", "status": status, "rows": len(df_cons)})
    except Exception as e:
        logger.error(f"Failed to load conservation features: {e}")
    
    # 4. Load or compute regulatory features (domains + gnomAD)
    try:
        from src.features.regulatory import RegulatoryFeatureComputer
        reg_path = FEATURES_DIR / "regulatory_features.parquet"
        status = "cache"
        if reg_path.exists():
            logger.info(f"Loading cached regulatory features from {reg_path}")
            df_reg = pd.read_parquet(reg_path)
        else:
            logger.info("Computing regulatory features (domains + gnomAD)...")
            computer = RegulatoryFeatureComputer(annotations_dir=ANNOTATIONS_DIR, output_dir=FEATURES_DIR)
            if computer.run():
                df_reg = pd.read_parquet(reg_path)
                status = "computed"
            else:
                logger.warning("Regulatory feature computation failed; using fallback")
                df_reg = pd.DataFrame()
                status = "fallback"
        feature_sources.append((df_reg, 'variant_id'))
        feature_logs.append({"feature": "regulatory", "status": status, "rows": len(df_reg)})
    except Exception as e:
        logger.error(f"Failed to load regulatory features: {e}")
    
    # Join all feature sources
    for df_features, join_key in feature_sources:
        if not df_features.empty and join_key in df_features.columns:
            df_scored_step3 = df_scored_step3.merge(df_features, on=join_key, how='left', suffixes=('', '_feature'))
    
    # Add LoF prior based on consequence
    _consequence_lof = {
        "frameshift_variant": 0.95,
        "stop_gained": 0.95,
        "splice_acceptor_variant": 0.95,
        "splice_donor_variant": 0.95,
        "missense_variant": 0.1,
        "synonymous_variant": 0.01,
    }

    if "lof_prior" not in df_scored_step3.columns:
        consequence_col = next((c for c in df_scored_step3.columns if 'consequence' in c.lower()), 'vep_consequence')
        df_scored_step3["lof_prior"] = df_scored_step3.get(consequence_col, "missense_variant").apply(
            lambda x: _consequence_lof.get(str(x).lower() if pd.notna(x) else "missense_variant", 0.1)
        )

    logger.info(f"Loaded all features. {len(df_scored_step3)} variants with {len(df_scored_step3.columns)} total columns.")
    df_scored_step3.attrs['feature_logs'] = feature_logs

    return df_scored_step3


@app.cell
def __(
    pd, df_scored_step3, logger, FEATURES_DIR
):
    """Save intermediate features."""
    _features_raw_path = FEATURES_DIR / "variants_features_raw.parquet"
    df_scored_step3.to_parquet(_features_raw_path)
    logger.info(f"Wrote raw features to {_features_raw_path}")
    return _features_raw_path


@app.cell
def __(mo, pd, df_scored_step3):
    """Summarize feature-source status (cache vs computed)."""
    logs = df_scored_step3.attrs.get('feature_logs', [])
    if not logs:
        mo.md("No feature status metadata available.")
        table = pd.DataFrame()
    else:
        table = pd.DataFrame(logs)
    mo.md("### Feature Source Status")
    mo.ui.table(table)


@app.cell
def __(mo):
    """
    ## Step 4: Impact Score Construction

    Choose between two modes: hand-mix or logistic regression.
    """
    mo.md(__doc__)

    scoring_mode_widget = mo.ui.radio(
        options=["hand-mix", "logistic"],
        value="hand-mix",
        label="Scoring Mode"
    )

    return scoring_mode_widget


@app.cell
def __(mo, scoring_mode_widget):
    """Guidance on when to switch modes."""
    if scoring_mode_widget.value == "hand-mix":
        mo.md("Hand-mix is ideal for first-pass sanity checks. Dial the weights until pathogenic vs benign distributions behave, then flip to logistic when you're ready for a fixed model.")
    else:
        mo.md("Logistic mode trains a calibrated score using current ClinVar labels. Re-run hand-mix if you want to experiment before retraining.")


@app.cell
def __(mo):
    """Document the v1 decision on impact scoring."""
    mo.md("""
> **v1 decision:** ship the ABCA4 pipeline with the hand-mix weights above. Logistic regression remains available for v1.1 once we collect more curated LP/B labels. Anytime we regenerate `variants_scored.parquet` for v1 we should keep the radio on *hand-mix* so downstream notebooks stay deterministic.
""")


@app.cell
def __(mo, scoring_mode_widget):
    """Create hand-mix weight sliders."""
    if scoring_mode_widget.value == "hand-mix":
        alpha_wgt = mo.ui.slider(0, 1, value=0.4, step=0.05, label="AlphaMissense Weight")
        splice_wgt = mo.ui.slider(0, 1, value=0.3, step=0.05, label="SpliceAI Weight")
        cons_wgt = mo.ui.slider(0, 1, value=0.15, step=0.05, label="Conservation Weight")
        lof_wgt = mo.ui.slider(0, 1, value=0.15, step=0.05, label="LoF Prior Weight")
    else:
        alpha_wgt = None
        splice_wgt = None
        cons_wgt = None
        lof_wgt = None

    return alpha_wgt, splice_wgt, cons_wgt, lof_wgt


@app.cell
def __(scoring_mode_widget, df_scored_step3, logger, pd, np):
    """
    Train logistic regression model on ClinVar labels (LP/P vs B/LB).
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    logistic_model = None
    logistic_scaler = None

    if scoring_mode_widget.value == "logistic":
        try:
            # Extract features and labels
            train_feature_cols = [c for c in df_scored_step3.columns if any(
                x in c.lower() for x in ['alphamissense', 'spliceai', 'phylop', 'conservation', 'lof_prior']
            )]

            if not train_feature_cols:
                raise ValueError("No feature columns available for logistic regression")

            # Define labels: Pathogenic/Likely Pathogenic vs Benign/Likely Benign
            def _is_pathogenic(sig):
                if not isinstance(sig, str):
                    return 0
                sig_lower = sig.lower()
                return 1 if ('pathogenic' in sig_lower and 'benign' not in sig_lower) else 0

            clinsig_col = next((c for c in df_scored_step3.columns if 'clinical_significance' in c.lower()), None)
            if clinsig_col:
                y = df_scored_step3[clinsig_col].apply(_is_pathogenic).values
            else:
                y = np.zeros(len(df_scored_step3))

            # Prepare data
            X_train = df_scored_step3[train_feature_cols].fillna(0.0).values
            logistic_scaler = StandardScaler()
            X_train_scaled = logistic_scaler.fit_transform(X_train)

            # Train logistic regression
            logistic_model = LogisticRegression(random_state=42, max_iter=1000)
            logistic_model.fit(X_train_scaled, y)

            logger.info(f"Trained logistic regression with {len(train_feature_cols)} features")
            logger.info(f"Model accuracy: {logistic_model.score(X_train_scaled, y):.3f}")

            # Display coefficients
            coef_display = pd.DataFrame({
                'Feature': train_feature_cols,
                'Coefficient': logistic_model.coef_[0],
                'Abs_Coeff': np.abs(logistic_model.coef_[0])
            }).sort_values('Abs_Coeff', ascending=False)

            logger.info("\nTop features by logistic regression coefficient:")
            logger.info(coef_display.head(10).to_string())

        except Exception as e:
            logger.error(f"Logistic regression training failed: {e}")
            logistic_model = None
            logistic_scaler = None

    return logistic_model, logistic_scaler


@app.cell
def __(
    pd, np, logger,
    df_scored_step3, scoring_mode_widget,
    alpha_wgt, splice_wgt, cons_wgt, lof_wgt,
    logistic_model, logistic_scaler
):
    """Compute impact scores using hand-mix or logistic regression."""
    
    df_impact = df_scored_step3.copy()

    if scoring_mode_widget.value == "hand-mix" and alpha_wgt is not None:
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

        _total_wgt = (
            alpha_wgt.value + splice_wgt.value +
            cons_wgt.value + lof_wgt.value
        )
        if _total_wgt == 0:
            _total_wgt = 1.0

        df_impact["model_score"] = (
            (alpha_wgt.value * _alpha_norm +
             splice_wgt.value * _splice_norm +
             cons_wgt.value * _cons_norm +
             lof_wgt.value * _lof_norm) / _total_wgt
        )

        logger.info(f"Computed hand-mix impact scores.")

    elif logistic_model is not None and logistic_scaler is not None:
        # Use trained logistic regression model
        try:
            # Extract same features used in training
            predict_feature_cols = [c for c in df_impact.columns if any(
                x in c.lower() for x in ['alphamissense', 'spliceai', 'phylop', 'conservation', 'lof_prior']
            )]

            X_predict = df_impact[predict_feature_cols].fillna(0.0).values
            X_predict_scaled = logistic_scaler.transform(X_predict)
            
            # Get probability of pathogenic
            proba = logistic_model.predict_proba(X_predict_scaled)
            df_impact["model_score"] = proba[:, 1]  # Probability of pathogenic (class 1)
            
            logger.info(f"Computed logistic regression impact scores. Mean: {df_impact['model_score'].mean():.3f}")
        except Exception as e:
            logger.error(f"Logistic regression scoring failed: {e}; using uniform fallback")
            df_impact["model_score"] = np.random.uniform(0, 1, len(df_impact))
    else:
        # Fallback: uniform random
        df_impact["model_score"] = np.random.uniform(0, 1, len(df_impact))
        logger.info("Using uniform random fallback for impact scores")

    return df_impact


@app.cell
def __():
    """Import plotly for visualizations."""
    import plotly.graph_objects as go
    return go


@app.cell
def __(mo, go, df_impact):
    """Visualize impact score distribution with plotly."""
    if "model_score" in df_impact.columns and not df_impact.empty:
        try:
            _fig_score = go.Figure()
            _fig_score.add_trace(go.Histogram(
                x=df_impact["model_score"].dropna(),
                nbinsx=30,
                name="Impact Score",
                marker_color="rgba(99, 110, 250, 0.7)"
            ))
            _fig_score.update_layout(
                title="Impact Score Distribution",
                xaxis_title="Impact Score",
                yaxis_title="Frequency",
                hovermode="x unified",
                showlegend=False,
                template="plotly_white"
            )
            mo.ui.plotly(_fig_score)
        except Exception as _e:
            mo.md(f"Visualization error: {_e}")
    else:
        mo.md("Impact scores not available.")


@app.cell
def __(mo):
    """
    ## Step 5: Clustering & Coverage Targets

    Define clusters and compute coverage thresholds.
    """
    mo.md(__doc__)

    clustering_widget = mo.ui.radio(
        options=["domain", "consequence", "manual"],
        value="domain",
        label="Clustering Strategy"
    )

    return clustering_widget


@app.cell
def __(
    pd, np, logger,
    df_impact, clustering_widget
):
    """Assign cluster membership."""
    df_clusters = df_impact.copy()

    if clustering_widget.value == "domain":
        if "domain" in df_clusters.columns:
            df_clusters["cluster"] = df_clusters["domain"].fillna("unknown")
        else:
            df_clusters["cluster"] = "unknown"
        logger.info(f"Domain-based clustering: {df_clusters['cluster'].nunique()} clusters")

    elif clustering_widget.value == "consequence":
        if "consequence" in df_clusters.columns:
            df_clusters["cluster"] = df_clusters["consequence"].fillna("unknown")
        else:
            df_clusters["cluster"] = "unknown"
        logger.info(f"Consequence-based clustering: {df_clusters['cluster'].nunique()} clusters")

    else:
        df_clusters["cluster"] = df_clusters.get("domain", "unknown").fillna("unknown")

    return df_clusters


@app.cell
def __(mo, df_clusters):
    """Display cluster membership."""
    if "cluster" in df_clusters.columns and not df_clusters.empty:
        _cluster_counts = df_clusters["cluster"].value_counts().to_frame("count")
        mo.md(f"""
### Cluster Membership

**Total clusters:** {df_clusters['cluster'].nunique()}
""")
        mo.ui.table(_cluster_counts.reset_index())


@app.cell
def __(
    pd, logger,
    df_clusters
):
    """Compute cluster coverage targets."""
    _cluster_targets = {}

    for _cluster_name, _group in df_clusters.groupby("cluster"):
        # Count pathogenic variants
        if "clinical_significance" in _group.columns:
            def _is_pathogenic(sig):
                return "pathogenic" in str(sig).lower() and "benign" not in str(sig).lower()
            _n_pathogenic = _group["clinical_significance"].apply(_is_pathogenic).sum()
        else:
            _n_pathogenic = 0

        _n_total = len(_group)
        _max_score = _group.get("model_score", pd.Series([0.0])).max()

        _cluster_targets[_cluster_name] = {
            "n_variants": _n_total,
            "n_pathogenic": _n_pathogenic,
            "max_score": _max_score,
            "tau_j": _max_score * 0.8,
        }

    logger.info(f"Computed coverage targets for {len(_cluster_targets)} clusters")

    return _cluster_targets


@app.cell
def __(mo, pd, df_clusters):
    """Display cluster targets."""
    _cluster_targets_display = {}

    for _cluster_name, _group in df_clusters.groupby("cluster"):
        if "clinical_significance" in _group.columns:
            def _is_path(sig):
                return "pathogenic" in str(sig).lower() and "benign" not in str(sig).lower()
            _n_path = _group["clinical_significance"].apply(_is_path).sum()
        else:
            _n_path = 0

        _n_tot = len(_group)
        _max_sc = _group.get("model_score", pd.Series([0.0])).max()

        _cluster_targets_display[_cluster_name] = {
            "n_variants": _n_tot,
            "n_pathogenic": _n_path,
            "max_score": _max_sc,
            "tau_j": _max_sc * 0.8,
        }

    _targets_df = pd.DataFrame([
        {
            "Cluster": k,
            "Variants": v["n_variants"],
            "Pathogenic": v["n_pathogenic"],
            "Max Score": f"{v['max_score']:.3f}",
            "Target τⱼ": f"{v['tau_j']:.3f}",
        }
        for k, v in _cluster_targets_display.items()
    ])
    mo.md("### Coverage Targets per Cluster")
    mo.ui.table(_targets_df)


@app.cell
def __(
    pd, df_clusters
):
    """Add cluster info and finalize."""
    df_final_scored = df_clusters.copy()

    # Compute cluster targets inline
    _cluster_tgt_dict = {}
    for _cn, _cg in df_clusters.groupby("cluster"):
        _mx = _cg.get("model_score", pd.Series([0.0])).max()
        _cluster_tgt_dict[_cn] = _mx * 0.8

    df_final_scored["cluster_target"] = df_final_scored["cluster"].map(
        lambda c: _cluster_tgt_dict.get(c, 0.5)
    )
    return df_final_scored


@app.cell
def __(
    logger, df_final_scored, FEATURES_DIR
):
    """Save final scored and clustered variants."""
    _final_path = FEATURES_DIR / "variants_scored.parquet"
    df_final_scored.to_parquet(_final_path)
    logger.info(f"Wrote scored & clustered variants to {_final_path}")
    return _final_path


@app.cell
def __(mo, logger, df_final_scored, FEATURES_DIR):
    """Confirm completion and save."""
    _final_path_confirm = FEATURES_DIR / "variants_scored.parquet"
    df_final_scored.to_parquet(_final_path_confirm)
    logger.info(f"Wrote scored & clustered variants")
    
    mo.md(f"""
✅ **Feature Engineering Complete!**

Saved to: `{_final_path_confirm}`

**Next Step:** Open `03_optimization_dashboard.py` for Strand optimization.
""")


if __name__ == "__main__":
    app.run()
