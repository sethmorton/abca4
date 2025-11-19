"""
MAVE (Multiplexed Assay of Variant Effect) data integration layer.

This module provides tools to ingest, normalize, and integrate MaveDB datasets
with the variant selection pipeline.

Architecture:
- ingest.py: Load raw MaveDB exports into canonical schema
- normalize.py: Normalize scores, define hits, apply hit labels
- adapter.py: Convert MAVE variants to pipeline input format
- strategies.py: Selection strategies (random, baselines, Strand)
- eval.py: Evaluation metrics (recall, precision, coverage)
- plots.py: Visualization helpers
"""

# Canonical MAVE row schema (used across all datasets after ingestion)
MAVE_CANONICAL_SCHEMA = {
    "dataset_id": "str",           # Short name: BRCA1_SaturationMutagenesis2018
    "gene": "str",                  # Gene symbol: BRCA1
    "wt_aa": "str",                 # Wild-type amino acid (one letter): M
    "pos": "int",                   # Protein position: 1775
    "mut_aa": "str",                # Mutant amino acid (one letter): K
    "variant_id": "str",            # Stable ID: BRCA1:p.Met1775Lys or BRCA1:1775:M>K
    "functional_score_raw": "float", # Raw score from assay
    "functional_score": "float",    # Normalized score (standardized across dataset)
    "assay_type": "str",            # e.g., yeast_complementation, mammalian_growth
    "mavedb_accession": "str",      # Link to source: MAVE000042
    "metadata_json": "str",         # Optional free-form metadata (JSON string)
    # Added during normalization:
    "is_hit": "bool",               # True if variant meets hit definition
    "hit_label": "str",             # e.g., "lof" or "non_lof" for clarity
    # Added during feature pipeline:
    "impact_score": "float",        # Model-derived pathogenicity score
    "conservation": "float",        # PhyloP or similar
    "domain": "str",                # Protein domain annotation
    "cluster_id": "int",            # Clustering group for selection
    # After selection:
    "selected": "bool",             # Whether this variant was selected
    "rank": "int",                  # Rank in selection (1 = top)
}

# Interpretation notes:
# functional_score: 
#   - After normalization, interpretation depends on hit_definition.direction
#   - If low_is_loss_of_function: lower scores indicate loss of function
#   - If high_is_gain_of_function: higher scores indicate gain of function
#   - A variant is a "hit" if it meets the percentile_threshold in the appropriate tail

# Normalization methods:
#   - zscore: (x - mean) / std across dataset
#   - quantile: Map to [0, 1] by rank; clip_low/high apply before mapping
#   - minmax: (x - min) / (max - min)
#   - identity: No normalization, copy raw score

__all__ = [
    "MAVE_CANONICAL_SCHEMA",
]


