#!/usr/bin/env python3
"""
Gene-agnostic ClinVar variant filtering and processing.

This module loads gene-specific ClinVar variant data (extracted via extract_clinvar_variants.py)
and optionally filters to a specific gene, outputting standardized variant tables for 
downstream processing.

Usage:
    # Process ABCA4 variants
    python -m src.data.filter_clinvar_variants --gene ABCA4
    
    # Process TP53 variants from custom path
    python -m src.data.filter_clinvar_variants --gene TP53 --input-dir data_raw
    
    # Process any gene (will use {gene}_clinvar_vus.parquet from input_dir/variants/)
    python -m src.data.filter_clinvar_variants --gene BRCA1
"""

import logging
from pathlib import Path
from typing import Optional
import pandas as pd
from src.data.constants import ABCA4_GENE_SYMBOL

logger = logging.getLogger(__name__)


class ClinVarVariantFilter:
    """Gene-agnostic class for filtering and processing ClinVar variants."""
    
    def __init__(self, input_dir: Path, output_dir: Path, gene_name: str):
        """
        Initialize the ClinVar variant filter.
        
        Args:
            input_dir: Directory containing extracted gene variants (default: data_processed/variants/)
            output_dir: Directory to save processed variants
            gene_name: Gene symbol to process (e.g., 'ABCA4', 'TP53', 'BRCA1')
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.gene_name = gene_name
        self.output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Initialized ClinVarVariantFilter for gene: {self.gene_name}")
    
    def run(self) -> bool:
        """Run the filtering pipeline. Return True on success."""
        try:
            df = self._load_variants()
            if df is None or df.empty:
                logger.warning(f"No variants loaded for {self.gene_name}")
                return False
            
            logger.info(f"Loaded {len(df)} variants for {self.gene_name}")
            
            # Optional: filter by gene if gene_symbol column exists
            df_filtered = self._filter_by_gene(df)
            
            # Save processed variants
            self._save_variants(df_filtered)
            return True
        except Exception as e:
            logger.error(f"Variant filtering failed: {e}")
            return False
    
    def _load_variants(self) -> Optional[pd.DataFrame]:
        """Load variants from pre-extracted gene-specific ClinVar files."""
        # Primary location: data_processed/variants/{gene}_clinvar_vus.parquet
        parquet_file = self.input_dir / f"{self.gene_name.lower()}_clinvar_vus.parquet"
        
        if parquet_file.exists():
            try:
                logger.info(f"Loading from {parquet_file}")
                return pd.read_parquet(parquet_file)
            except Exception as e:
                logger.warning(f"Could not load {parquet_file}: {e}")
        
        # Fallback: try CSV format
        csv_file = self.input_dir / f"{self.gene_name.lower()}_clinvar_vus.csv"
        if csv_file.exists():
            try:
                logger.info(f"Loading from {csv_file}")
                return pd.read_csv(csv_file)
            except Exception as e:
                logger.warning(f"Could not load {csv_file}: {e}")
        
        logger.error(f"No variant files found for {self.gene_name} in {self.input_dir}")
        logger.info(f"Expected file: {parquet_file}")
        logger.info(f"First extract variants using: python extract_clinvar_variants.py --gene {self.gene_name}")
        return None
    
    def _filter_by_gene(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Filter variants to the gene of interest (if applicable).
        
        Args:
            df: Dataframe of variants
            
        Returns:
            Filtered dataframe (or original if no gene_symbol column)
        """
        if "gene_symbol" in df.columns:
            df_filtered = df[df["gene_symbol"] == self.gene_name].copy()
            logger.info(f"Filtered to {len(df_filtered)} variants for {self.gene_name}")
            return df_filtered
        else:
            # If no gene_symbol column, assume all variants are already for the target gene
            logger.info(f"No gene_symbol column found; assuming all {len(df)} variants are for {self.gene_name}")
            return df.copy()
    
    def _save_variants(self, df: pd.DataFrame) -> None:
        """Save filtered variants to output directory."""
        output_file = self.output_dir / f"{self.gene_name.lower()}_clinvar_vus.parquet"
        df.to_parquet(output_file, index=False)
        logger.info(f"✅ Saved {len(df)} variants to {output_file}")


def main():
    """CLI entry point for gene-agnostic ClinVar variant filtering."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Filter and process ClinVar variants for any gene',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process ABCA4 variants (default)
  python -m src.data.filter_clinvar_variants
  
  # Process TP53 variants
  python -m src.data.filter_clinvar_variants --gene TP53
  
  # Process from custom directories
  python -m src.data.filter_clinvar_variants --gene BRCA1 --input-dir data_processed/variants --output-dir data_processed/variants_filtered
    """
    )
    
    parser.add_argument(
        '--gene',
        type=str,
        default=ABCA4_GENE_SYMBOL,
        help=f'Gene symbol to process (default: {ABCA4_GENE_SYMBOL})'
    )
    parser.add_argument(
        '--input-dir',
        type=Path,
        default=Path('data_processed/variants'),
        help='Input directory with extracted gene variant files (default: data_processed/variants/)'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('data_processed/variants'),
        help='Output directory for processed variants (default: data_processed/variants/)'
    )
    
    args = parser.parse_args()
    
    # Create filter for the specified gene
    filter_engine = ClinVarVariantFilter(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        gene_name=args.gene
    )
    
    success = filter_engine.run()
    return 0 if success else 1


if __name__ == "__main__":
    import sys
    sys.exit(main())

