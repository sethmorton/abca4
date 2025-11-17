#!/usr/bin/env python3
"""
Filter ClinVar variants for ABCA4 gene.

Loads ClinVar VCF/TSV data and filters for ABCA4 variants with
Uncertain significance classifications.
"""

import pandas as pd
from pathlib import Path
import sys
from typing import Optional, List
from src.config import logger

CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]

class ABCA4VariantFilter:
    """Filter ClinVar data for ABCA4 variants of interest."""

    # ABCA4 filtering criteria
    GENE_NAME = "ABCA4"
    DEFAULT_CLINSIG_TERMS = (
        "uncertain significance",
        "conflicting interpretations",
        "conflicting interpretations of pathogenicity"
    )

    def __init__(self, input_dir: Optional[Path] = None, output_dir: Optional[Path] = None,
                 clinsig_terms: Optional[List[str]] = None):
        data_root = CAMPAIGN_ROOT / "data_raw"
        processed_root = CAMPAIGN_ROOT / "data_processed" / "variants"
        self.input_dir = input_dir or data_root
        self.output_dir = output_dir or processed_root
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.clinsig_terms = [term.lower() for term in (clinsig_terms or list(self.DEFAULT_CLINSIG_TERMS))]

    @staticmethod
    def _tokenize_clinsig(value: str) -> List[str]:
        tokens: List[str] = []
        if not isinstance(value, str):
            return tokens
        for raw_token in value.replace('|', ',').split(','):
            token = raw_token.strip().lower()
            if token:
                tokens.append(token)
        return tokens

    def _matches_clinsig(self, value: str) -> bool:
        tokens = self._tokenize_clinsig(value)
        return any(token in self.clinsig_terms for token in tokens)

    def _normalize_clinsig_value(self, value: str) -> str:
        tokens = self._tokenize_clinsig(value)
        for token in tokens:
            if token in self.clinsig_terms:
                return token
        return tokens[0] if tokens else value

    def load_clinvar_tsv(self) -> Optional[pd.DataFrame]:
        """Load ClinVar variant summary TSV."""
        tsv_path = self.input_dir / "clinvar" / "variant_summary.txt.gz"

        if not tsv_path.exists():
            logger.error(f"ClinVar TSV not found: {tsv_path}")
            return None

        logger.info(f"Loading ClinVar TSV: {tsv_path}")

        try:
            # Load TSV with appropriate dtypes
            df = pd.read_csv(
                tsv_path,
                sep='\t',
                compression='gzip',
                low_memory=False,
                dtype={
                    'Chromosome': str,
                    'Start': int,
                    'Stop': int,
                    'ReferenceAllele': str,
                    'AlternateAllele': str,
                    'GeneSymbol': str,
                    'ClinicalSignificance': str,
                    'ReviewStatus': str
                }
            )

            logger.info(f"Loaded {len(df)} ClinVar variants")
            return df

        except Exception as e:
            logger.error(f"Failed to load ClinVar TSV: {e}")
            return None

    def filter_abca4_variants(self, clinvar_df: pd.DataFrame) -> pd.DataFrame:
        """Filter for ABCA4 variants meeting our criteria."""
        logger.info("Filtering for ABCA4 variants...")

        # Filter by gene name
        abca4_df = clinvar_df[clinvar_df['GeneSymbol'] == self.GENE_NAME].copy()
        logger.info(f"Found {len(abca4_df)} ABCA4 variants")

        # Normalize and filter by clinical significance tokens
        mask = abca4_df['ClinicalSignificance'].apply(self._matches_clinsig)
        filtered_df = abca4_df[mask].copy()

        logger.info(f"Filtered to {len(filtered_df)} variants with target clinical significance")

        # Prefer VCF alleles (they have actual ref/alt), fallback to regular columns
        filtered_df['ref_to_use'] = filtered_df['ReferenceAlleleVCF'].fillna(filtered_df['ReferenceAllele'])
        filtered_df['alt_to_use'] = filtered_df['AlternateAlleleVCF'].fillna(filtered_df['AlternateAllele'])
        
        # Add variant ID and clean up columns
        filtered_df['variant_id'] = (
            filtered_df['Chromosome'].astype(str) + '_' +
            filtered_df['Start'].astype(str) + '_' +
            filtered_df['ref_to_use'] + '_' +
            filtered_df['alt_to_use']
        )

        # Select and rename columns for our analysis
        columns_to_keep = [
            'variant_id',
            'Chromosome',
            'Start',
            'Stop',
            'ref_to_use',
            'alt_to_use',
            'GeneSymbol',
            'ClinicalSignificance',
            'ReviewStatus',
            'RS# (dbSNP)',
            'PhenotypeList',
            'Origin',
            'Assembly'
        ]

        # Only keep columns that exist
        available_columns = [col for col in columns_to_keep if col in filtered_df.columns]
        filtered_df = filtered_df[available_columns].copy()

        # Rename columns for consistency
        column_rename = {
            'Chromosome': 'chrom',
            'Start': 'pos',
            'Stop': 'end',
            'ref_to_use': 'ref',
            'alt_to_use': 'alt',
            'GeneSymbol': 'gene',
            'ClinicalSignificance': 'clinical_significance',
            'ReviewStatus': 'review_status',
            'RS# (dbSNP)': 'rs_id',
            'PhenotypeList': 'phenotypes',
            'Origin': 'origin',
            'Assembly': 'assembly'
        }

        filtered_df = filtered_df.rename(columns=column_rename)
        filtered_df['clinical_significance'] = filtered_df['clinical_significance'].apply(
            self._normalize_clinsig_value
        )

        return filtered_df

    def add_genomic_context(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Add genomic context information."""
        logger.info("Adding genomic context...")

        # Add sequence context windows for ML models
        # ABCA4 region: chr1:94,400,000-95,200,000
        variants_df['region_start'] = variants_df['pos'] - 500  # 500bp upstream
        variants_df['region_end'] = variants_df['pos'] + 500    # 500bp downstream

        # Add variant type classification
        variants_df['variant_type'] = 'unknown'
        variants_df.loc[
            (variants_df['ref'].str.len() == 1) & (variants_df['alt'].str.len() == 1),
            'variant_type'
        ] = 'SNP'

        variants_df.loc[
            (variants_df['ref'].str.len() != variants_df['alt'].str.len()),
            'variant_type'
        ] = 'indel'

        # Count variants by type
        type_counts = variants_df['variant_type'].value_counts()
        logger.info(f"Variant types: {dict(type_counts)}")

        return variants_df

    def save_filtered_variants(self, variants_df: pd.DataFrame) -> bool:
        """Save filtered variants to disk with schema and stats artifacts."""
        base_name = "abca4_clinvar_vus"
        output_path = self.output_dir / f"{base_name}.parquet"

        try:
            variants_df.to_parquet(output_path, index=False)
            logger.info(f"Saved {len(variants_df)} ABCA4 variants to {output_path}")

            # Also save as CSV for easier inspection
            csv_path = self.output_dir / f"{base_name}.csv"
            variants_df.to_csv(csv_path, index=False)
            logger.info(f"Also saved as CSV: {csv_path}")

            # Save schema and statistics as JSON artifacts
            self.save_schema_artifacts(variants_df, base_name)

            return True

        except Exception as e:
            logger.error(f"Failed to save variants: {e}")
            return False

    def save_schema_artifacts(self, variants_df: pd.DataFrame, base_name: str) -> None:
        """Save schema and statistics artifacts for downstream validation."""
        import json
        from datetime import datetime

        # Schema information
        schema_info = {
            "dataset": base_name,
            "created_at": datetime.now().isoformat(),
            "total_rows": len(variants_df),
            "total_columns": len(variants_df.columns),
            "columns": {}
        }

        # Column information
        for col in variants_df.columns:
            col_info = {
                "dtype": str(variants_df[col].dtype),
                "nullable": variants_df[col].isnull().any(),
                "null_count": int(variants_df[col].isnull().sum()),
                "unique_values": int(variants_df[col].nunique()) if variants_df[col].dtype == 'object' else None
            }

            # Add range info for numeric columns
            if pd.api.types.is_numeric_dtype(variants_df[col]):
                col_info.update({
                    "min": float(variants_df[col].min()) if not variants_df[col].empty else None,
                    "max": float(variants_df[col].max()) if not variants_df[col].empty else None,
                    "mean": float(variants_df[col].mean()) if not variants_df[col].empty else None,
                    "std": float(variants_df[col].std()) if not variants_df[col].empty else None
                })

            schema_info["columns"][col] = col_info

        # Statistics summary
        stats_summary = {
            "dataset": base_name,
            "created_at": datetime.now().isoformat(),
            "summary_stats": self.generate_summary_stats(variants_df)
        }

        # Save artifacts
        schema_path = self.output_dir / f"{base_name}_schema.json"
        stats_path = self.output_dir / f"{base_name}_stats.json"

        with open(schema_path, 'w') as f:
            json.dump(schema_info, f, indent=2, default=str)
        logger.info(f"Saved schema artifact: {schema_path}")

        with open(stats_path, 'w') as f:
            json.dump(stats_summary, f, indent=2, default=str)
        logger.info(f"Saved statistics artifact: {stats_path}")

    def generate_summary_stats(self, variants_df: pd.DataFrame) -> dict:
        """Generate summary statistics."""
        stats = {
            'total_variants': len(variants_df),
            'unique_positions': variants_df['pos'].nunique(),
            'chromosomes': variants_df['chrom'].unique().tolist(),
            'clinical_significance_breakdown': variants_df['clinical_significance'].value_counts().to_dict(),
            'variant_type_breakdown': variants_df['variant_type'].value_counts().to_dict(),
            'review_status_breakdown': variants_df['review_status'].value_counts().to_dict() if 'review_status' in variants_df.columns else {},
            'has_rs_id': variants_df['rs_id'].notna().sum() if 'rs_id' in variants_df.columns else 0,
        }

        return stats

    def run(self) -> bool:
        """Run the complete filtering process."""
        logger.info("Starting ABCA4 variant filtering process...")
        logger.info(f"Gene: {self.GENE_NAME}")
        logger.info(f"Clinical significance filter: {self.clinsig_terms}")

        # Load ClinVar data
        clinvar_df = self.load_clinvar_tsv()
        if clinvar_df is None:
            return False

        # Filter for ABCA4 variants
        abca4_df = self.filter_abca4_variants(clinvar_df)
        if len(abca4_df) == 0:
            logger.warning("No ABCA4 variants found matching criteria")
            return False

        # Add genomic context
        abca4_df = self.add_genomic_context(abca4_df)

        # Generate and log summary stats
        stats = self.generate_summary_stats(abca4_df)
        logger.info("Summary statistics:")
        for key, value in stats.items():
            logger.info(f"  {key}: {value}")

        # Save results
        if not self.save_filtered_variants(abca4_df):
            return False

        logger.info("ABCA4 variant filtering process completed successfully!")
        return True

def main():
    """Main entry point."""
    filter = ABCA4VariantFilter()
    success = filter.run()
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
