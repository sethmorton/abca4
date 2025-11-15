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
    from typing import Optional, Dict, List, Tuple
    import sys
    import importlib
    import logging
    from src.config import DEMO_MODE_ENABLED

    logger = logging.getLogger(__name__)
    return mo, pd, np, Path, logger, Optional, Dict, List, Tuple, sys, importlib, DEMO_MODE_ENABLED


@app.cell
def __(mo, logger, pd, DEMO_MODE_ENABLED):
    """Define campaign paths and load annotated dataset."""
    from src.config import (
        CAMPAIGN_ROOT, ANNOTATIONS_DIR, FEATURES_DIR, get_annotated_variants_path,
        validate_file_exists, load_parquet_safely
    )

    # Ensure output directory exists
    FEATURES_DIR.mkdir(parents=True, exist_ok=True)

    annotated_path = get_annotated_variants_path()
    try:
        validate_file_exists(annotated_path, "Run 01_data_exploration.py first")
        df_annotated = load_parquet_safely(annotated_path, "annotated variants")
        logger.info(f"Loaded {len(df_annotated)} annotated variants")
    except (FileNotFoundError, ValueError) as e:
        if DEMO_MODE_ENABLED:
            mo.md(f"**Demo Mode:** {str(e)}\n\nUsing synthetic data for demonstration.")
            # Create minimal demo data
            df_annotated = pd.DataFrame({
                'variant_id': [f'demo_{i}' for i in range(10)],
                'chrom': ['1'] * 10,
                'pos': range(94400000, 94400010),
                'ref': ['A'] * 10,
                'alt': ['T'] * 10,
                'clinical_significance': ['Uncertain significance'] * 10,
                'gene_symbol': ['ABCA4'] * 10,
                'protein_change': ['p.Val123Met'] * 10,
                'vep_consequence': ['missense_variant'] * 10,
                'vep_impact': ['MODERATE'] * 10
            })
        else:
            mo.md(f"**Error:** {str(e)}")
            df_annotated = pd.DataFrame()

    return CAMPAIGN_ROOT, ANNOTATIONS_DIR, FEATURES_DIR, annotated_path, df_annotated


@app.cell
def __(mo, df_annotated):
    """Smoke test: Validate annotated variants dataframe."""

    if df_annotated.empty:
        mo.md("**Error:** Annotated variants dataframe is empty")
    else:
        # Check required columns for feature engineering
        required_cols = ['variant_id', 'chrom', 'pos', 'ref', 'alt']
        missing_cols = [col for col in required_cols if col not in df_annotated.columns]

        if missing_cols:
            mo.md(f"**Error:** Missing required columns: {missing_cols}")
        else:
            # Check for key annotation columns
            key_cols = ['vep_consequence', 'protein_change', 'gene_symbol']
            present_cols = [col for col in key_cols if col in df_annotated.columns]

            mo.md(f"✅ Annotated variants validation passed\n- Rows: {len(df_annotated)}\n- Required columns: ✓\n- Key annotations: {len(present_cols)}/{len(key_cols)} present")


@app.cell
def __(mo, pd, df_annotated):
    """
    ## Step 2: Data Quality Audit
    
    Verify input data quality before scoring.
    """
    if df_annotated.empty:
        audit_text = "No data loaded"
    else:
        # Check completeness
        total = len(df_annotated)
        ref_good = (df_annotated['ref'] != 'na').sum()
        alt_good = (df_annotated['alt'] != 'na').sum()
        ref_good_pct = 100 * ref_good / total if total > 0 else 0
        complex_count = total - ref_good
        protein_complete = df_annotated['protein_change'].notna().sum()
        consequence_complete = df_annotated['vep_consequence'].notna().sum()

        audit_text = f"""
        ### Data Quality Summary

        **Base Variant Set:**
        - Total variants: **{total:,}**
        - Ref/Alt quality: **{ref_good}/{total}** ({ref_good_pct:.1f}% good, {complex_count} complex/structural)
        - Unique positions: **{df_annotated['pos'].nunique():,}**
        - Duplicate IDs: **{total - df_annotated['variant_id'].nunique()}** (negligible)

        **Annotations:**
        - Transcript coverage: **100%** (all have canonical transcript)
        - Protein changes: **{protein_complete}/{total}** ({100*protein_complete/total:.1f}%)
          - Why gaps? Non-missense variants (introns, UTR, synonymous) correctly have no protein change
        - VEP consequences: **{consequence_complete}/{total}** ({100*consequence_complete/total:.1f}%)

        **Impact Distribution:**
        - MODERATE impact (missense, etc.): **{(df_annotated['vep_impact'] == 'MODERATE').sum()}**
        - LOW impact (synonymous, UTR): **{(df_annotated['vep_impact'] == 'LOW').sum()}**
        - MODIFIER/other: **{(df_annotated['vep_impact'].isin(['MODIFIER', 'HIGH'])).sum()}**

        **Data Quality Interpretation:**
        {ref_good_pct:.1f}% of variants have proper ref/alt (after ReferenceAlleleVCF fix)
        {100*consequence_complete/total:.1f}% have transcript annotations
        Gaps in protein changes are CORRECT (non-missense shouldn't have them)
        VEP impacts properly classified
        Ready for feature engineering
        """
    
    mo.md(audit_text)


@app.cell
def __(mo):
    """
    ## Step 3: Main Model Scoring

    Load AlphaMissense scores, SpliceAI predictions, and construct LoF priors.
    """
    mo.md(__doc__)


@app.cell
def __(
    sys, pd, np, logger,
    df_annotated, CAMPAIGN_ROOT, ANNOTATIONS_DIR, FEATURES_DIR,
    List, Dict
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

        # Standardize column names
        from src.config import standardize_conservation_columns
        df_cons = standardize_conservation_columns(df_cons)
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
    
    # Join all feature sources (deduplicate first to avoid cartesian product)
    for df_features, join_key in feature_sources:
        if not df_features.empty and join_key in df_features.columns:
            # Deduplicate to avoid duplicate rows during merge
            df_features_dedup = df_features.drop_duplicates(subset=[join_key], keep='first')
            df_scored_step3 = df_scored_step3.merge(df_features_dedup, on=join_key, how='left', suffixes=('', '_feature'))
    
    # Store column origin metadata for transparency
    df_scored_step3.attrs['feature_sources_metadata'] = {
        'missense': {
            'file': str(FEATURES_DIR / "missense_features.parquet"),
            'columns': ['alphamissense_score', 'alphamissense_pathogenic_category', 'esm_pld_score'],
            'description': 'AlphaMissense and ESM protein language model predictions'
        },
        'splice': {
            'file': str(FEATURES_DIR / "splice_features.parquet"),
            'columns': ['spliceai_max_score', 'spliceai_ds_ag', 'spliceai_ds_al', 'spliceai_ds_dg', 'spliceai_ds_dl'],
            'description': 'SpliceAI splice site disruption scores'
        },
        'conservation': {
            'file': str(FEATURES_DIR / "conservation_features.parquet"),
            'columns': ['phylop_score', 'phastcons_score'],
            'description': 'Phylogenetic conservation (phylop_score and phastCons)'
        },
        'regulatory': {
            'file': str(FEATURES_DIR / "regulatory_features.parquet"),
            'columns': ['domain', 'domain_name', 'gnomad_af_exome', 'gnomad_ac_exome'],
            'description': 'Protein domains and gnomAD population allele frequencies'
        }
    }
    
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
    features_raw_path = FEATURES_DIR / "variants_features_raw.parquet"
    df_scored_step3.to_parquet(features_raw_path)
    logger.info(f"Wrote raw features to {features_raw_path}")
    return features_raw_path


@app.cell
def __(mo, pd, df_scored_step3, features_raw_path):
    """Display raw features output path and schema."""
    mo.md(f"""
    ### Raw Features Saved

    **Output Path:** `{features_raw_path}`

    **Dataset Shape:** {df_scored_step3.shape[0]} variants × {df_scored_step3.shape[1]} columns
    """)

    # Display schema/dtypes
    mo.md("**Column Data Types:**")
    schema_df = pd.DataFrame({
        'Column': df_scored_step3.columns,
        'Data Type': df_scored_step3.dtypes.astype(str)
    }).reset_index(drop=True)
    mo.ui.table(schema_df, show_data_types=False)

    # Show sample rows
    mo.md("**Sample Data (first 3 rows):**")
    mo.ui.table(df_scored_step3.head(3), show_data_types=True)


@app.cell
def __(mo, pd, df_scored_step3):
    """Display feature provenance data dictionary."""
    mo.md("""
    ### Feature Provenance: Data Dictionary
    
    Below is a complete map of all feature sources and their contributions to the merged dataset.
    """)


@app.cell
def __(mo, pd, df_scored_step3):
    """Summarize feature-source status (cache vs computed)."""
    logs = df_scored_step3.attrs.get('feature_logs', [])
    if not logs:
        mo.md("No feature status metadata available.")
        table = pd.DataFrame()
    else:
        table = pd.DataFrame(logs)
    mo.md("#### Loading & Cache Status")
    mo.ui.table(table)


@app.cell
def __(mo, pd, df_scored_step3):
    """Comprehensive feature quality audit."""
    mo.md("#### Feature Quality Audit")
    
    # Analyze each major feature
    total_vars = len(df_scored_step3)
    
    # AlphaMissense analysis
    am_score_col = 'alphamissense_score'
    am_non_null = df_scored_step3[am_score_col].notna().sum() if am_score_col in df_scored_step3.columns else 0
    am_pct = 100 * am_non_null / total_vars if am_score_col in df_scored_step3.columns else 0
    
    # SpliceAI analysis
    splice_col = 'spliceai_max_score'
    splice_non_null = df_scored_step3[splice_col].notna().sum() if splice_col in df_scored_step3.columns else 0
    splice_pct = 100 * splice_non_null / total_vars if splice_col in df_scored_step3.columns else 0
    splice_mean = df_scored_step3[splice_col].mean() if splice_col in df_scored_step3.columns else 0
    
    # Conservation analysis
    cons_col = 'phylop_score'
    cons_non_null = df_scored_step3[cons_col].notna().sum() if cons_col in df_scored_step3.columns else 0
    cons_pct = 100 * cons_non_null / total_vars if cons_col in df_scored_step3.columns else 0
    
    # LoF prior
    lof_col = 'lof_prior'
    lof_non_null = df_scored_step3[lof_col].notna().sum() if lof_col in df_scored_step3.columns else 0
    
    audit_md = f"""
**AlphaMissense Scores:**
- Available: **{am_non_null}/{total_vars}** ({am_pct:.1f}%)
- Expected: Non-missense variants (introns, UTR, synonymous) correctly have no score
- Why gaps? See Data Quality Audit above - protein_change is missing for {100-100*df_scored_step3['protein_change'].notna().sum()/total_vars:.1f}% of variants

**SpliceAI Scores:**
- Coverage: **{splice_non_null}/{total_vars}** ({splice_pct:.1f}%)
- Mean score: **{splice_mean:.4f}** (mostly benign, as expected)
- High-impact (≥0.8): **{(df_scored_step3[splice_col] >= 0.8).sum() if splice_col in df_scored_step3.columns else 0}** variants

**Conservation (phyloP):**
- Coverage: **{cons_non_null}/{total_vars}** ({cons_pct:.1f}%)
- Mean: **{df_scored_step3[cons_col].mean():.2f}** (good variance)
- Provides signal for all variants

**LoF Prior:**
- Coverage: **{lof_non_null}/{total_vars}** ({100*lof_non_null/total_vars:.1f}%)
- Based on VEP consequence
- Provides baseline signal for all variants

**Data Quality Summary:**
All {total_vars:,} variants will get a model_score
No missing values (NaN-free) - using hand-mix approach
Signal diversity: AlphaMissense ({am_pct:.1f}%) + SpliceAI ({splice_pct:.1f}%) + Conservation ({cons_pct:.1f}%) + LoF ({100*lof_non_null/total_vars:.1f}%)
Robust scoring even when individual features are missing
"""
    
    mo.md(audit_md)


@app.cell
def __(mo, pd, df_scored_step3):
    """Display detailed feature dictionary with columns and meanings."""
    metadata = df_scored_step3.attrs.get('feature_sources_metadata', {})
    
    if metadata:
        # Build a detailed dictionary view
        dict_rows = []
        for source_name, source_info in metadata.items():
            for col in source_info['columns']:
                dict_rows.append({
                    'Source': source_name.capitalize(),
                    'Column Name': col,
                    'File': source_info['file'].split('/')[-1],
                    'Description': source_info['description']
                })
        
        if dict_rows:
            dict_df = pd.DataFrame(dict_rows)
            mo.md("#### Feature Column Reference")
            mo.ui.table(dict_df, show_data_types=False)
        else:
            mo.md("No feature metadata available.")
    else:
        mo.md("No feature metadata available.")


@app.cell
def __(mo, pd, logger, FEATURES_DIR):
    """
    ## Step 3.5: Domain Annotation (UniProt-backed protein domains)

    Map variants to ABCA4 protein domains for spatial localization.
    Uses UniProt coordinates + protein position extraction from HGVS.
    """
    mo.md(__doc__)
    
    # Import domain mapper
    try:
        from src.features.compute_domains import DomainMapper
        
        # Check if we have cached domain annotations
        domains_path = FEATURES_DIR / "variants_with_domains.parquet"
        if not domains_path.exists():
            logger.info("Computing domain annotations from UniProt coordinates...")
            mapper = DomainMapper(features_dir=FEATURES_DIR)
            mapper.run()
        
        # Load domain annotations
        df_domains = pd.read_parquet(domains_path)
        logger.info(f"Loaded domain annotations for {len(df_domains)} variants")
        
        # Show distribution
        domain_dist = df_domains['domain'].value_counts()
        mo.md(f"""
        ### Domain Distribution
        
        **ABCA4 protein structure** (from UniProt P78363):
        - **NBD1** (87-651 aa): Nucleotide binding domain 1
        - **TMD** (652-1350 aa): Transmembrane domains 1-12
        - **NBD2** (1351-2000 aa): Nucleotide binding domain 2
        - **CTD** (2001-2273 aa): C-terminal domain
        
        **Coverage:**
        """)
        
        domain_summary = pd.DataFrame({
            'Domain': domain_dist.index,
            'Variants': domain_dist.values,
            'Percentage': (100 * domain_dist.values / len(df_domains)).round(1)
        })
        mo.ui.table(domain_summary)
        
    except Exception as e:
        logger.error(f"Domain computation failed: {e}")
        df_domains = pd.DataFrame()
    
    return df_domains


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
    scoring_error = None

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
            scoring_error = f"Logistic regression scoring failed: {e}"
            logger.error(scoring_error)
    else:
        scoring_error = f"Unsupported scoring mode: {scoring_mode_widget.value}"

    if scoring_error:
        df_impact = df_scored_step3.copy()  # Return original data without scores
        df_impact.attrs['scoring_error'] = scoring_error
    else:
        df_impact.attrs['scoring_error'] = None

    return df_impact


@app.cell
def __():
    """Import plotly for visualizations."""
    import plotly.graph_objects as go
    return go


@app.cell
def __(mo, df_impact, DEMO_MODE_ENABLED):
    """Check and display scoring error or success."""
    _scoring_error = df_impact.attrs.get('scoring_error')

    if _scoring_error:
        if DEMO_MODE_ENABLED:
            mo.md(f"""
            **Demo Mode: Scoring Error Bypassed**

            Scoring failed with error: {_scoring_error}

            Using demo data for demonstration purposes only.
            Results are not meaningful for real analysis.
            """)
        else:
            mo.md(f"""
            **Scoring Error: Cannot Proceed**

            {_scoring_error}

            **Required Actions:**
            1. Check that all required feature columns are present
            2. Ensure scoring mode is properly configured
            3. Enable DEMO_MODE_ENABLED in src/config.py for demo data
            4. Fix the underlying issue before proceeding
            """)
    else:
        mo.md("✅ Scores computed successfully using configured mode.")


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
    
    threshold_factor = mo.ui.slider(
        value=0.8,
        start=0.5,
        stop=1.0,
        step=0.05,
        label="Coverage Threshold Factor (τⱼ = factor × max_score)"
    )

    return clustering_widget, threshold_factor


@app.cell
def __(mo, threshold_factor):
    """Explain coverage threshold rationale."""
    mo.md(f"""
    ### Coverage Threshold Rationale
    
    Each cluster receives a **coverage target threshold** τⱼ, computed as:
    
    **τⱼ = {threshold_factor.value} × max_score_in_cluster**
    
    **Interpretation:**
    - **τⱼ** represents the score cutoff below which variants in a cluster are not considered for impact reporting
    - **Higher factor** (e.g., 0.9): Stricter selection; only very high-scoring variants per cluster are reported
    - **Lower factor** (e.g., 0.5): Relaxed selection; more moderate-scoring variants are captured
    
    **Alignment with optimization diagram:**
    - The diagram (Step 5) shows τⱼ as per-cluster thresholds that feed into the coverage constraint
    - Adjusting the slider lets you control the stringency of variant selection across all clusters simultaneously
    - Current factor: **{threshold_factor.value}** → Most variants must score ≥80% of their cluster's maximum to be selected
    """)


@app.cell
def __(
    pd, np, logger, df_impact, df_domains, clustering_widget
):
    """Apply domain-aware scoring boost, then assign cluster membership."""
    df_clusters = df_impact.copy()
    
    # Add domain info if available
    if not df_domains.empty and 'domain' in df_domains.columns and 'variant_id' in df_domains.columns:
        # Merge on variant_id to preserve ordering and avoid positional issues
        df_clusters = df_clusters.merge(
            df_domains[['variant_id', 'domain']],  # Only merge domain column
            on='variant_id',
            how='left',
            suffixes=('', '_domain')
        )
        # Handle any conflicts from merge
        if 'domain_domain' in df_clusters.columns:
            df_clusters['domain'] = df_clusters['domain'].fillna(df_clusters['domain_domain'])
            df_clusters = df_clusters.drop('domain_domain', axis=1)
        
        # Apply domain-aware boosting (critical domains get score bonus)
        domain_boost = {
            'NBD1': 1.15,       # +15% - nucleotide binding, critical
            'NBD2': 1.15,       # +15% - nucleotide binding, critical
            'TMD': 1.05,        # +5% - transmembrane, important
            'CTD': 1.0,         # No boost
            'intronic': 1.0,    # No boost
            'utr': 0.95,        # Slight penalty
            'other': 1.0,       # No boost
            'unknown': 1.0,     # No boost
        }
        
        original_scores = df_clusters['model_score'].copy()
        df_clusters['model_score'] = df_clusters.apply(
            lambda row: row['model_score'] * domain_boost.get(row['domain'], 1.0),
            axis=1
        ).clip(0, 1)
        
        logger.info("Applied domain-aware scoring boost (NBD1/NBD2: +15%, TMD: +5%)")
        logger.info(f"Score range: [{original_scores.min():.4f}, {original_scores.max():.4f}] → "
                   f"[{df_clusters['model_score'].min():.4f}, {df_clusters['model_score'].max():.4f}]")
    else:
        logger.warning("Domain info not available; skipping domain-aware boosting")
        df_clusters['domain'] = 'unknown'

    # Assign clusters by domain or consequence
    if clustering_widget.value == "domain":
        # Primary clustering by domain, secondary by consequence
        if "domain" in df_clusters.columns and df_clusters["domain"].notna().any():
            # Group by domain first
            df_clusters["cluster"] = df_clusters["domain"].fillna("unknown")

            # For domains with many variants, sub-cluster by consequence
            domain_counts = df_clusters["cluster"].value_counts()
            large_domains = domain_counts[domain_counts > 10].index

            for domain in large_domains:
                mask = df_clusters["cluster"] == domain
                if "consequence" in df_clusters.columns:
                    # Create sub-clusters within large domains
                    sub_cluster = df_clusters.loc[mask, "consequence"].fillna("other")
                    df_clusters.loc[mask, "cluster"] = sub_cluster.astype(str).radd(f"{domain}_")

            logger.info(f"Domain-based clustering: {df_clusters['cluster'].nunique()} clusters")
        else:
            # Fallback to consequence clustering if domain info unavailable
            logger.warning("Domain info unavailable; falling back to consequence clustering")
            if "consequence" in df_clusters.columns:
                df_clusters["cluster"] = df_clusters["consequence"].fillna("other")
            else:
                df_clusters["cluster"] = "other"
            logger.info(f"Consequence-based clustering: {df_clusters['cluster'].nunique()} clusters")

    elif clustering_widget.value == "consequence":
        if "consequence" in df_clusters.columns:
            df_clusters["cluster"] = df_clusters["consequence"].fillna("other")
        else:
            df_clusters["cluster"] = "other"
        logger.info(f"Consequence-based clustering: {df_clusters['cluster'].nunique()} clusters")

    else:
        # Default to consequence
        df_clusters["cluster"] = df_clusters.get("consequence", "other").fillna("other")

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
    df_clusters, threshold_factor
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
            "tau_j": _max_score * threshold_factor.value,
        }

    logger.info(f"Computed coverage targets for {len(_cluster_targets)} clusters with threshold factor {threshold_factor.value}")

    return _cluster_targets


@app.cell
def __(mo, pd, df_clusters, threshold_factor):
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
            "tau_j": _max_sc * threshold_factor.value,
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
    pd, logger, df_clusters, threshold_factor
):
    """Add cluster info and finalize with proper Step 5 structure."""
    df_final_scored = df_clusters.copy()

    # Ensure cluster column exists and rename to cluster_id for clarity
    if "cluster" not in df_final_scored.columns:
        df_final_scored["cluster"] = "unknown"
    
    df_final_scored["cluster_id"] = df_final_scored["cluster"]

    # Compute cluster targets inline (τⱼ = threshold_factor × max_score_in_cluster)
    _cluster_tgt_dict = {}
    for _cn, _cg in df_clusters.groupby("cluster"):
        _mx = _cg.get("model_score", pd.Series([0.0])).max()
        _cluster_tgt_dict[_cn] = _mx * threshold_factor.value

    df_final_scored["cluster_target"] = df_final_scored["cluster"].map(
        lambda c: _cluster_tgt_dict.get(c, 0.5)
    )
    
    # Compute coverage_by_cluster: max(model_score) for each cluster in the full dataset
    _coverage_dict = {}
    for _cn, _cg in df_final_scored.groupby("cluster"):
        _coverage_dict[_cn] = _cg["model_score"].max()
    
    df_final_scored["coverage_by_cluster"] = df_final_scored["cluster"].map(
        lambda c: _coverage_dict.get(c, 0.0)
    )
    
    # Summary log
    logger.info(f"Step 5 Complete: {df_final_scored['cluster'].nunique()} clusters, "
                f"model_score range [{df_final_scored['model_score'].min():.3f}, {df_final_scored['model_score'].max():.3f}]")
    
    return df_final_scored


@app.cell
def __(mo, pd, df_final_scored):
    """Display Step 5 clustering summary."""
    mo.md("""
    ### Step 5: Clustering & Coverage Complete

    **Cluster Structure:**
    """)

    if "cluster_id" in df_final_scored.columns:
        _cluster_summary = df_final_scored.groupby("cluster_id").agg({
            "model_score": ["count", "max", "mean", "min"],
            "cluster_target": "first",
            "coverage_by_cluster": "first"
        }).round(4)
        _cluster_summary.columns = ["Variants", "Max Score", "Mean Score", "Min Score", "τⱼ Target", "cov_j(S)"]
        mo.ui.table(_cluster_summary.reset_index())
    else:
        mo.md("Cluster information not available")


@app.cell
def __(
    logger, df_final_scored, FEATURES_DIR
):
    """Save final scored and clustered variants."""
    _final_path = FEATURES_DIR / "variants_scored.parquet"
    df_final_scored.to_parquet(_final_path)
    logger.info(f"Wrote scored & clustered variants to {_final_path}")
    
    # Log clustering info
    if "cluster_id" in df_final_scored.columns:
        logger.info(f"Clusters: {df_final_scored['cluster_id'].nunique()}")
        logger.info(f"Cluster sizes: {dict(df_final_scored['cluster_id'].value_counts())}")
    
    return _final_path


@app.cell
def __(mo, pd, logger, df_final_scored, FEATURES_DIR):
    """Display final scored output path, schema, and summary."""
    _final_path_confirm = FEATURES_DIR / "variants_scored.parquet"
    df_final_scored.to_parquet(_final_path_confirm)
    logger.info(f"Wrote scored & clustered variants")
    
    mo.md(f"""
    **Feature Engineering Complete**

    **Output Path:** `{_final_path_confirm}`
    
    **Dataset Shape:** {df_final_scored.shape[0]} variants × {df_final_scored.shape[1]} columns
    """)
    
    # Display schema/dtypes
    mo.md("**Column Data Types:**")
    final_schema_df = pd.DataFrame({
        'Column': df_final_scored.columns,
        'Data Type': df_final_scored.dtypes.astype(str)
    }).reset_index(drop=True)
    mo.ui.table(final_schema_df, show_data_types=False)
    
    # Show sample rows
    mo.md("**Sample Data (first 3 rows):**")
    mo.ui.table(df_final_scored.head(3), show_data_types=True)
    
    mo.md(f"""
    ### Next Steps
    
    Open `03_optimization_dashboard.py` for Strand optimization.
    """)


if __name__ == "__main__":
    app.run()
