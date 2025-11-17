"""
Generate Assay Drafts using LLM

Reads selected variants panel, validates data contract, and generates
assay protocol drafts using Groq LLM with provenance tracking.
"""

import json
import time
from pathlib import Path
from typing import Dict, Any, List

import pandas as pd

from ..config import (
    logger, REQUIRED_VARIANT_COLUMNS, ASSAY_DRAFTS_DIR,
    LLM_MAX_VARIANTS, RUN_ID, validate_dataframe
)
from .llm_client import generate_assay_markdown, create_llm_config, AssayDraftError, get_prompt_hash


def validate_selected_variants_panel(df: pd.DataFrame) -> None:
    """
    Validate selected variants panel meets data contract.

    Args:
        df: Selected variants dataframe

    Raises:
        ValueError: If validation fails
    """
    logger.info(f"Validating selected variants panel with {len(df)} variants")

    # Check for empty panel
    if df.empty:
        raise ValueError("Selected variants panel is empty")

    # Check required columns
    missing_cols = set(REQUIRED_VARIANT_COLUMNS) - set(df.columns)
    if missing_cols:
        raise ValueError(f"Selected variants panel missing required columns: {missing_cols}")

    # Check for null values in required columns
    null_counts = df[REQUIRED_VARIANT_COLUMNS].isnull().sum()
    if null_counts.any():
        null_cols = null_counts[null_counts > 0].index.tolist()
        raise ValueError(f"Selected variants panel has null values in required columns: {null_cols}")

    # Validate impact_score range (should be 0-1 based on optimization)
    if not df['impact_score'].between(0, 1).all():
        raise ValueError("impact_score values must be between 0 and 1")

    # Validate gnomad_max_af is numeric
    if not pd.api.types.is_numeric_dtype(df['gnomad_max_af']):
        raise ValueError("gnomad_max_af must be numeric")

    logger.info("Selected variants panel validation passed")


def prepare_variant_for_prompt(row: pd.Series) -> Dict[str, Any]:
    """
    Convert dataframe row to variant dict for prompt template.

    Args:
        row: Pandas Series representing a variant

    Returns:
        Dictionary with variant information for template rendering
    """
    # Get consequence and cluster_id
    consequence = str(row['vep_consequence'])
    cluster_id = str(row['cluster_id'])

    # Derive mechanism using a two-step approach:
    # 1. Try domain-based mapping from cluster_id
    # 2. Fall back to consequence-based mapping

    mechanism = 'unspecified_functional_impact'

    # Step 1: Domain-based mechanism from cluster_id
    # cluster_id format is typically "{domain}_{consequence_description}"
    cluster_parts = cluster_id.split('_', 1)
    if len(cluster_parts) >= 1:
        domain = cluster_parts[0].lower()
        domain_to_mechanism = {
            'nbd1': 'nucleotide_binding_disturbance',
            'nbd2': 'nucleotide_binding_disturbance',
            'tmd1': 'transmembrane_transport_disruption',
            'tmd2': 'transmembrane_transport_disruption',
            'ctd': 'cytoplasmic_domain_dysfunction',
        }
        if domain in domain_to_mechanism:
            mechanism = domain_to_mechanism[domain]

    # Step 2: Consequence-based fallback/refinement
    consequence_lower = consequence.lower()
    if 'missense' in consequence_lower:
        if mechanism == 'unspecified_functional_impact':
            mechanism = 'protein_function_alteration'
    elif 'splice' in consequence_lower or 'intron' in consequence_lower:
        mechanism = 'splicing_alteration'
    elif 'stop' in consequence_lower or 'frameshift' in consequence_lower:
        mechanism = 'loss_of_function'
    elif 'synonymous' in consequence_lower:
        mechanism = 'minimal_functional_impact'

    # Handle protein_change - may not exist in all datasets
    protein_change = 'N/A (intronic variant)'
    if 'protein_change' in row.index and pd.notna(row['protein_change']):
        protein_change = str(row['protein_change'])

    return {
        'variant_id': str(row['variant_id']),
        'gene': str(row['gene']),
        'protein_change': protein_change,
        'consequence': consequence,
        'domain_cluster_id': cluster_id,
        'mechanism': mechanism,
        'impact_score': float(row['impact_score']),
        'gnomad_max_af': float(row['gnomad_max_af'])
    }


def generate_assay_drafts(selected_variants_path: Path, output_dir: Path) -> Dict[str, Any]:
    """
    Generate assay drafts for selected variants using LLM.

    Args:
        selected_variants_path: Path to selected variants CSV/Parquet
        output_dir: Directory to write outputs

    Returns:
        Summary statistics and metadata

    Raises:
        SystemExit: On validation failure or critical errors
    """
    logger.info(f"Starting assay draft generation from {selected_variants_path}")

    # Validate configuration
    llm_config = create_llm_config()

    # Load and validate selected variants
    if selected_variants_path.suffix == '.parquet':
        df = pd.read_parquet(selected_variants_path)
    else:
        df = pd.read_csv(selected_variants_path)

    try:
        validate_selected_variants_panel(df)
    except ValueError as e:
        logger.error(f"Selected variants panel validation failed: {e}")
        raise SystemExit(1) from e

    # Sort by impact_score descending and limit to MAX_VARIANTS
    df_sorted = df.sort_values('impact_score', ascending=False).head(LLM_MAX_VARIANTS)
    logger.info(f"Processing top {len(df_sorted)} variants by impact score")

    # Create output directories
    protocol_dir = output_dir / "protocol_drafts"
    protocol_dir.mkdir(parents=True, exist_ok=True)

    # Initialize tracking
    successful_drafts = 0
    failed_drafts = 0
    total_tokens = 0
    draft_metadata = []
    start_time = time.time()

    # Process each variant
    for rank, (_, row) in enumerate(df_sorted.iterrows(), 1):
        variant_id = row['variant_id']
        logger.info(f"Generating draft for variant {variant_id} (rank {rank})")

        try:
            # Prepare variant data for prompt
            variant_data = prepare_variant_for_prompt(row)

            # Generate assay markdown
            assay_markdown, prompt_hash, tokens_used = generate_assay_markdown(variant_data, llm_config)

            # Save individual draft
            filename = f"{rank:02d}_{variant_id}.md"
            draft_path = protocol_dir / filename

            with open(draft_path, 'w', encoding='utf-8') as f:
                f.write(f"# Assay Draft for {variant_id}\n\n")
                f.write(f"**Generated:** {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"**Model:** {llm_config.model}\n")
                f.write(f"**Temperature:** {llm_config.temperature}\n")
                f.write(f"**Max Tokens:** {llm_config.max_tokens}\n\n")
                f.write("---\n\n")
                f.write(assay_markdown)
                f.write("\n\n---\n\n*LLM-generated assay protocol - human review required*")

            # Track metadata
            draft_metadata.append({
                'variant_id': variant_id,
                'filename': filename,
                'model': llm_config.model,
                'temperature': llm_config.temperature,
                'max_tokens': llm_config.max_tokens,
                'prompt_hash': prompt_hash,
                'tokens_used': tokens_used,
                'timestamp': time.time(),
                'git_sha': RUN_ID,
                'rank': rank,
                'impact_score': float(row['impact_score']),
                'status': 'success'
            })

            successful_drafts += 1

        except AssayDraftError as e:
            logger.error(f"Failed to generate draft for {variant_id}: {e}")
            draft_metadata.append({
                'variant_id': variant_id,
                'filename': None,
                'model': llm_config.model,
                'temperature': llm_config.temperature,
                'max_tokens': llm_config.max_tokens,
                'prompt_hash': None,
                'timestamp': time.time(),
                'git_sha': RUN_ID,
                'rank': rank,
                'impact_score': float(row['impact_score']),
                'status': 'failed',
                'error': str(e)
            })
            failed_drafts += 1
            continue

    # Generate index file
    index_path = output_dir / "assay_drafts_index.json"
    with open(index_path, 'w', encoding='utf-8') as f:
        json.dump({
            'metadata': {
                'run_id': RUN_ID,
                'timestamp': time.time(),
                'total_variants': len(df_sorted),
                'successful_drafts': successful_drafts,
                'failed_drafts': failed_drafts,
                'llm_config': {
                    'model': llm_config.model,
                    'temperature': llm_config.temperature,
                    'max_tokens': llm_config.max_tokens
                }
            },
            'drafts': draft_metadata
        }, f, indent=2, default=str)

    # Generate CSV index for easier reading
    index_csv_path = output_dir / "assay_drafts_index.csv"
    index_df = pd.DataFrame(draft_metadata)
    index_df.to_csv(index_csv_path, index=False)

    # Summary statistics
    duration = time.time() - start_time
    summary = {
        'total_variants_processed': len(df_sorted),
        'successful_drafts': successful_drafts,
        'failed_drafts': failed_drafts,
        'total_duration_seconds': duration,
        'average_time_per_variant': duration / len(df_sorted) if df_sorted else 0,
        'output_directory': str(output_dir),
        'run_id': RUN_ID
    }

    logger.info(f"Assay draft generation complete: {successful_drafts} success, {failed_drafts} failed")

    # Fail pipeline if any drafts failed (required step)
    if failed_drafts > 0:
        logger.error(f"Pipeline failed: {failed_drafts} assay drafts could not be generated")
        raise SystemExit(1)

    return summary


def main():
    """CLI entry point for assay draft generation."""
    import sys
    from ..config import REPORTS_DIR

    if len(sys.argv) == 1:
        # No args provided, use default path
        selected_variants_path = REPORTS_DIR / "variants_selected.csv"
    elif len(sys.argv) == 2:
        # Path provided as argument
        selected_variants_path = Path(sys.argv[1])
    else:
        print("Usage: python -m src.reporting.generate_assay_drafts [selected_variants_path]")
        print("If no path provided, defaults to data_processed/reports/variants_selected.csv")
        sys.exit(1)

    if not selected_variants_path.exists():
        logger.error(f"Selected variants file not found: {selected_variants_path}")
        sys.exit(1)

    output_dir = ASSAY_DRAFTS_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        summary = generate_assay_drafts(selected_variants_path, output_dir)
        logger.info(f"Assay drafts generated successfully: {summary}")
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        logger.error(f"Unexpected error during assay draft generation: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
