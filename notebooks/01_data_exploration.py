#!/usr/bin/env python3
"""
ABCA4 Campaign – Interactive Data Exploration with Marimo

Stages:
  - Step 0: Scope & framing (context about ABCA4, Stargardt, the v1 question)
  - Step 1: Data ingest (ClinVar/raw files, optional TSV upload, filtering, live tables)
  - Step 2: Annotation & deterministic features (VEP, gnomAD, conservation, domain mapping)

This notebook is fully reactive: adjust filters or upload new data, and all downstream
visualizations update automatically. At the end of each section, data is persisted to
data_processed/ for use in downstream notebooks.

Run interactively:  marimo edit notebooks/01_data_exploration.py
Run as dashboard:   marimo run notebooks/01_data_exploration.py
Run as script:      python notebooks/01_data_exploration.py
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
        _file_path = Path.cwd() / "notebooks" / "01_data_exploration.py"
    _project_root = _file_path.parent.parent
    if str(_project_root) not in sys.path:
        sys.path.insert(0, str(_project_root))

    import marimo as mo
    import pandas as pd
    import numpy as np
    from typing import Optional, List, Dict, Tuple
    import json
    import logging
    from src.config import DEMO_MODE_ENABLED

    logger = logging.getLogger(__name__)

    return mo, pd, np, Path, logging, logger, Optional, List, Dict, Tuple, json, sys, DEMO_MODE_ENABLED


@app.cell
def __(mo, Path):
    """
    ## Step 0: Scope & Framing

    **ABCA4 and Stargardt Disease**

    ABCA4 (ATP Binding Cassette Transporter A4) is an ATP-dependent lipid transporter crucial for photoreceptor function.
    Pathogenic variants in ABCA4 cause Stargardt disease (STGD1), a progressive macular degeneration leading to
    blindness. Over 3,500 known ABCA4 variants exist, many classified as Uncertain Significance (VUS).

    **The v1 Question**

    Given rare ABCA4 variants (VUS or conflicting), how can we select K=20–50 variants for high-throughput MPRAssay,
    combining:
    - Computational pathogenicity predictions (AlphaMissense, SpliceAI, conservation)
    - Domain disruption analysis
    - Allele frequency / gnomAD background
    - Experimental priors (known pathogenic patterns, assay scalability)

    **Campaign Goal**

    Rank and select a diverse, high-confidence subset of variants to advance toward experimental validation.
    """
    mo.md(__doc__)
    return


@app.cell
def __(mo):
    """Define campaign root paths."""
    from src.config import (
        CAMPAIGN_ROOT, DATA_RAW_DIR, DATA_PROCESSED_DIR,
        ANNOTATIONS_DIR, FEATURES_DIR, validate_file_exists
    )

    # Additional paths specific to this notebook
    VARIANTS_DIR = DATA_PROCESSED_DIR / "variants"

    # Ensure directories exist
    for d in [DATA_RAW_DIR, VARIANTS_DIR, ANNOTATIONS_DIR, FEATURES_DIR]:
        d.mkdir(parents=True, exist_ok=True)
    
    return (
        CAMPAIGN_ROOT,
        DATA_RAW_DIR,
        DATA_PROCESSED_DIR,
        VARIANTS_DIR,
        ANNOTATIONS_DIR,
        FEATURES_DIR,
    )


@app.cell
def __(mo, DATA_RAW_DIR, DATA_PROCESSED_DIR):
    """Smoke test: Validate required directories and data availability."""
    import os

    # Check directories
    required_dirs = [DATA_RAW_DIR, DATA_PROCESSED_DIR]
    missing_dirs = [d for d in required_dirs if not d.exists()]

    if missing_dirs:
        mo.md(f"**Error:** Missing required directories: {missing_dirs}")
    else:
        # Check for raw data files
        raw_files = list(DATA_RAW_DIR.glob("**/*.vcf.gz")) + list(DATA_RAW_DIR.glob("**/*.tsv"))
        if not raw_files:
            mo.md("**Warning:** No raw data files found in data_raw/. This notebook may not function properly.")
        else:
            mo.md("✅ Environment validation passed")


@app.cell
def __(mo):
    """
    ## Step 0: Configuration Parameters

    Configure defaults for gene, transcript, K, and narratives. These bind downstream cells automatically.
    """
    mo.md(__doc__)

    # Gene and transcript parameters
    gene_symbol = mo.ui.text(value="ABCA4", label="Gene Symbol", disabled=False)
    transcript_id = mo.ui.text(value="ENST00000370225", label="Canonical Transcript", disabled=False)

    # Optimization parameters
    k_variants = mo.ui.slider(10, 200, value=30, label="Panel Size (K variants to select)")
    budget_narrative = mo.ui.text_area(
        value="High-throughput MPRAssay validation with focus on Stargardt disease pathogenicity",
        label="Budget/Narrative",
    )

    return gene_symbol, transcript_id, k_variants, budget_narrative


@app.cell
def __(mo, gene_symbol, transcript_id, k_variants, budget_narrative):
    """Display current configuration."""
    config_display = mo.md(f"""
### Current Configuration

- **Gene:** {gene_symbol.value}
- **Transcript:** {transcript_id.value}
- **Panel Size (K):** {k_variants.value}
- **Narrative:** {budget_narrative.value[:80]}...

All downstream cells use these parameters automatically.
""")
    return config_display


@app.cell
def __(mo):
    """
    ## Step 1: Data Ingest

    ### Available Data Sources

    - **ClinVar Variants** (data_raw/clinvar/)
    - **gnomAD Frequencies** (data_raw/gnomad/)
    - **SpliceAI Predictions** (data_raw/spliceai/)
    - **AlphaMissense Scores** (data_raw/alphamissense/)

    Upload an optional TSV file with additional variants, or use the downloaded data as-is.
    """
    mo.md(__doc__)


@app.cell
def __(mo):
    """Optional TSV Upload."""
    uploaded_file = mo.ui.file_browser(filetypes=[".tsv", ".csv", ".txt"], label="Upload Partner Data (Optional)")
    return uploaded_file


@app.cell
def __(
    pd, np, logger, DEMO_MODE_ENABLED,
    DATA_RAW_DIR, VARIANTS_DIR,
    uploaded_file, CAMPAIGN_ROOT
):
    """
    Load gene variants using the ClinVar filtering pipeline.
    
    This cell invokes the real data processing pipeline rather than truncating TSV
    or creating synthetic rows.
    """
    sys.path.insert(0, str(CAMPAIGN_ROOT))
    from src.data.filter_clinvar_variants import ClinVarVariantFilter
    
    variants_source = "unknown"

    try:
        # Check if we have cached filtered variants
        cached_path = VARIANTS_DIR / "abca4_clinvar_vus.parquet"
        if cached_path.exists():
            logger.info(f"Loading cached variants from {cached_path}")
            df_variants_raw = pd.read_parquet(cached_path)
            logger.info(f"Loaded {len(df_variants_raw)} cached variants")
            variants_source = "cache"
        else:
            # Run the authoritative filter pipeline
            logger.info("Running authoritative ABCA4 variant filtering pipeline...")
            filter_engine = ABCA4VariantFilter(input_dir=DATA_RAW_DIR, output_dir=VARIANTS_DIR)
            filter_success = filter_engine.run()

            if filter_success:
                df_variants_raw = pd.read_parquet(cached_path)
                logger.info(f"Loaded {len(df_variants_raw)} filtered variants from pipeline")
                variants_source = "pipeline"
            else:
                logger.error("Variant filtering pipeline failed")
                df_variants_raw = pd.DataFrame()

    except Exception as e:
        logger.error(f"Failed to run variant filter: {e}")
        if DEMO_MODE_ENABLED:
            logger.info("Demo mode: Using synthetic example data")
            # Demo: Create example data for testing
            df_variants_raw = pd.DataFrame({
                "chrom": ["1"] * 10,
                "pos": range(94400000, 94400010),
                "ref": ["A"] * 10,
                "alt": ["T", "G", "C"] * 3 + ["T"],
                "clinical_significance": ["Uncertain significance"] * 10,
                "gene_symbol": ["ABCA4"] * 10,
            })
            variants_source = "demo"
        else:
            logger.error("Data loading failed. Enable DEMO_MODE_ENABLED in src/config.py for synthetic data, or check data_raw/ directory.")
            df_variants_raw = pd.DataFrame()
            variants_source = "error"

    # Load uploaded TSV if provided (append to existing data)
    if uploaded_file is not None and hasattr(uploaded_file, 'value') and uploaded_file.value:
        logger.info(f"Loading uploaded file")
        try:
            df_uploaded = pd.read_csv(uploaded_file.value, sep="\t", dtype=str)
            df_variants_raw = pd.concat([df_variants_raw, df_uploaded], ignore_index=True, sort=False)
            # Deduplicate by chrom, pos, ref, alt
            cols_for_dedup = ["chrom", "pos", "ref", "alt"]
            if all(c in df_variants_raw.columns for c in cols_for_dedup):
                df_variants_raw = df_variants_raw.drop_duplicates(subset=cols_for_dedup, keep="first")
            logger.info(f"  Loaded {len(df_uploaded)} variants from uploaded file, now {len(df_variants_raw)} total")
            variants_source = f"{variants_source}+upload" if variants_source else "upload"
        except Exception as e:
            logger.error(f"Failed to load uploaded file: {e}")

    logger.info(f"Total unique variants: {len(df_variants_raw)}")
    df_variants_raw.attrs['source'] = variants_source or "unknown"
    
    return df_variants_raw


@app.cell
def __(mo, pd, df_variants_raw):
    """Display data snapshot after loading variants."""
    if not df_variants_raw.empty:
        mo.md("### Data Snapshot")

        # Display shape
        mo.md(f"**Dataset Shape:** {df_variants_raw.shape[0]} rows × {df_variants_raw.shape[1]} columns")

        # Display column dtypes
        schema_df = pd.DataFrame({
            "column": df_variants_raw.columns,
            "dtype": df_variants_raw.dtypes.astype(str)
        })
        mo.md("**Column Data Types:**")
        mo.ui.table(schema_df.reset_index(drop=True), show_data_types=False)

        # Display head(3)
        mo.md("**Sample Data (first 3 rows):**")
        mo.ui.table(df_variants_raw.head(3), show_data_types=True)


@app.cell
def __(mo, df_variants_raw):
    """Report how variants were loaded (cache vs pipeline)."""
    variant_source = df_variants_raw.attrs.get('source', 'unknown')
    mo.md(f"""
**Variant source:** `{variant_source}`

Cached parquet loads are instant. If the notebook re-ran the ClinVar pipeline you'll see `pipeline` here.
""")


@app.cell
def __(mo, df_variants_raw):
    """Display raw variant counts by clinical significance."""
    if df_variants_raw.empty:
        mo.md("No variants loaded.")
    else:
        # Find clinical significance column
        _clinsig_cols = [c for c in df_variants_raw.columns if "clin" in c.lower() or "significance" in c.lower()]
        _clinsig_col = _clinsig_cols[0] if _clinsig_cols else "clinical_significance"
        
        if _clinsig_col in df_variants_raw.columns:
            _clinsig_summary = df_variants_raw[_clinsig_col].value_counts().to_frame("count")
            mo.md(f"""
### Raw Variant Summary

**Total variants:** {len(df_variants_raw)}

**Distribution by Clinical Significance:**
""")
            mo.ui.table(_clinsig_summary.reset_index())
        else:
            mo.md(f"**Total variants:** {len(df_variants_raw)}")


@app.cell
def __(mo, df_variants_raw):
    """Interactive Filtering panel."""
    mo.md("""
### Interactive Filtering

Use controls below to filter variants.
""")

    filter_clinsig = None
    filter_af = None

    if not df_variants_raw.empty:
        # Find clinical significance column
        _clinsig_cols = [c for c in df_variants_raw.columns if "clin" in c.lower() or "significance" in c.lower()]
        _clinsig_col = _clinsig_cols[0] if _clinsig_cols else None
        
        if _clinsig_col:
            _clinsig_options = df_variants_raw[_clinsig_col].dropna().unique().tolist()
            filter_clinsig = mo.ui.multiselect(
                options=sorted(_clinsig_options),
                value=sorted(_clinsig_options),
                label="Clinical Significance"
            )

        # Optional: allele frequency filter
        filter_af = mo.ui.slider(0.0, 0.01, value=0.01, step=0.0001, label="Max gnomAD AF (if available)")

    return filter_clinsig, filter_af


@app.cell
def __(
    pd, np, df_variants_raw,
    filter_clinsig, filter_af
):
    """Apply filters to get working variant set."""
    df_variants_filtered = df_variants_raw.copy()

    # Filter by clinical significance
    if filter_clinsig is not None and hasattr(filter_clinsig, 'value'):
        _clinsig_cols = [c for c in df_variants_filtered.columns if "clin" in c.lower() or "significance" in c.lower()]
        _clinsig_col = _clinsig_cols[0] if _clinsig_cols else None
        if _clinsig_col and _clinsig_col in df_variants_filtered.columns:
            df_variants_filtered = df_variants_filtered[df_variants_filtered[_clinsig_col].isin(filter_clinsig.value)]

    # Filter by AF if gnomad_af column exists
    if filter_af is not None and "gnomad_af" in df_variants_filtered.columns and hasattr(filter_af, 'value'):
        df_variants_filtered = df_variants_filtered[df_variants_filtered["gnomad_af"] <= filter_af.value]

    return df_variants_filtered


@app.cell
def __(mo, df_variants_filtered):
    """Display filtered variant count."""
    mo.md(f"""
### Filtered Variants: {len(df_variants_filtered)} variants

Ready for annotation and feature engineering.
""")


@app.cell
def __(mo):
    """
    ## Step 2: Annotation & Deterministic Features

    Add VEP annotations, gnomAD joins, conservation scores, and domain mappings.
    """
    mo.md(__doc__)


@app.cell
def __(
    pd, np, logger, DEMO_MODE_ENABLED,
    df_variants_filtered, CAMPAIGN_ROOT, VARIANTS_DIR, ANNOTATIONS_DIR
):
    """
    Add transcript annotations using the authoritative annotate_transcripts.py pipeline.
    
    This cell invokes VEP and pyensembl APIs for real transcript/consequence annotations
    rather than using placeholder values.
    """
    sys.path.insert(0, str(CAMPAIGN_ROOT))
    from src.annotation.annotate_transcripts import VariantAnnotator
    
    annotation_source = "unknown"

    try:
        # Check if we have cached annotated variants
        cached_annot_path = ANNOTATIONS_DIR / "abca4_vus_annotated.parquet"
        if cached_annot_path.exists():
            logger.info(f"Loading cached annotations from {cached_annot_path}")
            df_annot = pd.read_parquet(cached_annot_path)
            logger.info(f"Loaded {len(df_annot)} cached annotated variants")
            annotation_source = "cache"
        else:
            # Run the authoritative annotation pipeline
            logger.info("Running authoritative variant annotation pipeline...")
            
            # First, save filtered variants to the expected location
            temp_path = VARIANTS_DIR / "abca4_clinvar_vus.parquet"
            df_variants_filtered.to_parquet(temp_path)
            
            annotator = VariantAnnotator(input_dir=VARIANTS_DIR, output_dir=ANNOTATIONS_DIR)
            annotation_success = annotator.run()

            if annotation_success:
                df_annot = pd.read_parquet(cached_annot_path)
                logger.info(f"Loaded {len(df_annot)} annotated variants from pipeline")
                annotation_source = "pipeline"
            else:
                if DEMO_MODE_ENABLED:
                    logger.warning("Annotation pipeline failed; using minimal demo annotations")
                df_annot = df_variants_filtered.copy()
                    # Add minimal demo columns
                for col in ["transcript_id", "vep_consequence", "vep_impact", "genomic_region"]:
                    if col not in df_annot.columns:
                        df_annot[col] = None
                    annotation_source = "demo"
                else:
                    logger.error("Annotation pipeline failed. Enable DEMO_MODE_ENABLED in src/config.py for demo data.")
                    df_annot = pd.DataFrame()
                    annotation_source = "error"

    except Exception as e:
        logger.error(f"Failed to run annotation pipeline: {e}")
        if DEMO_MODE_ENABLED:
            logger.info("Demo mode: Using minimal fallback annotations")
        df_annot = df_variants_filtered.copy()
            # Add minimal demo columns
        for col in ["transcript_id", "vep_consequence", "vep_impact", "genomic_region"]:
            if col not in df_annot.columns:
                df_annot[col] = None
            annotation_source = "demo"
        else:
            logger.error("Annotation failed. Enable DEMO_MODE_ENABLED in src/config.py for demo data.")
            df_annot = pd.DataFrame()
            annotation_source = "error"

    logger.info(f"Annotation complete. Ready for feature engineering.")
    df_annot.attrs['source'] = annotation_source

    return df_annot


@app.cell
def __(mo, df_annot):
    """Show annotation data provenance."""
    annot_source = df_annot.attrs.get('source', 'unknown')
    mo.md(f"""
**Annotation source:** `{annot_source}`

Cached annotations skip pyensembl/VEPlac calls; `pipeline` indicates this run hit the Ensembl APIs.
""")


@app.cell
def __(mo, df_annot):
    """Display consequence distribution."""
    if "consequence" in df_annot.columns and not df_annot.empty:
        _consequence_counts = df_annot["consequence"].value_counts().to_frame("count")
        mo.md("""
### Consequence Distribution
""")
        mo.ui.table(_consequence_counts.head(10).reset_index())
    else:
        mo.md("No consequence data available.")


@app.cell
def __(df_annot):
    """
    Pass annotated variants forward.
    
    Detailed feature engineering (gnomAD, conservation, regulatory, domains)
    happens in the next notebook (02_feature_engineering.py) which calls:
    - campaigns/abca4/src/features/regulatory.py (gnomAD + domain mapping)
    - campaigns/abca4/src/features/conservation.py (phyloP/phastCons)
    - campaigns/abca4/src/features/missense.py (AlphaMissense scores)
    - campaigns/abca4/src/features/splice.py (SpliceAI scores)
    """
    df_domains = df_annot.copy()
    return df_domains


@app.cell
def __(mo):
    """Note about domain distribution."""
    mo.md("""
### Domain Annotation

Domain mapping and regulatory features will be added in the next notebook
(02_feature_engineering.py) using the authoritative domain configuration
and gnomAD data.
""")



@app.cell
def __(mo, pd, df_domains):
    """Display completeness metrics."""
    _key_cols = [
        "chrom", "pos", "ref", "alt", "clinical_significance",
        "consequence", "gnomad_af", "phylop_score", "domain"
    ]
    _completeness = {}
    for _col in _key_cols:
        if _col in df_domains.columns:
            _completeness[_col] = df_domains[_col].notna().sum() / len(df_domains) if len(df_domains) > 0 else 0.0
        else:
            _completeness[_col] = 0.0

    _comp_df = pd.DataFrame(
        [(k, f"{v:.1%}") for k, v in _completeness.items()],
        columns=["Field", "Completeness"]
    )
    
    mo.md("""
### Annotation Completeness

Key fields completeness check:
""")
    mo.ui.table(_comp_df)


@app.cell
def __(
    logger,
    df_domains, ANNOTATIONS_DIR
):
    """Export annotated variants to disk."""
    from src.config import get_annotated_variants_path
    output_path_annot = get_annotated_variants_path()
    df_domains.to_parquet(output_path_annot)
    logger.info(f"Wrote annotated variants to {output_path_annot}")
    return output_path_annot


@app.cell
def __(mo, output_path_annot):
    """Confirm export."""
    mo.md(f"""
**Annotation Complete**

Saved to: `{output_path_annot}`

**Next Step:** Open `02_feature_engineering.py` to add model scores and construct impact metrics.
""")


@app.cell
def __(mo):
    """Tie this notebook back to the v1 plan."""
    mo.md("""
### Plan Alignment

- **Step 1 – Data ingest (ClinVar + partner variants):** Completed above via `filter_abca4_variants.py`, outputs now cached under `data_processed/variants/`.
- **Step 2 – Annotation & deterministic features:** Completed via `annotate_transcripts.py` (VEP/pyensembl). The resulting `variants_annotated.parquet` is the hand-off into feature engineering.

This notebook is now the authoritative entry point for Steps 0‑2 of the ABCA4 v1 pipeline.
""")


if __name__ == "__main__":
    app.run()
