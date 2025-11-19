#!/usr/bin/env python3
"""
Filter and process ClinVar variants for a specific gene.

This module loads ClinVar VCF/TSV data and filters to a gene of interest,
outputting standardized variant tables for downstream processing.
"""

import logging
from pathlib import Path
from typing import Optional
import pandas as pd
import json

logger = logging.getLogger(__name__)


class ClinVarVariantFilter:
    """Base class for filtering ClinVar variants."""
    
    def __init__(self, input_dir: Path, output_dir: Path, gene_name: str = "ABCA4"):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.gene_name = gene_name
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def run(self) -> bool:
        """Run the filtering pipeline. Return True on success."""
        try:
            df = self._load_variants()
            if df is None or df.empty:
                logger.warning(f"No variants loaded for {self.gene_name}")
                return False
            
            df_filtered = self._filter_by_gene(df)
            self._save_variants(df_filtered)
            return True
        except Exception as e:
            logger.error(f"Variant filtering failed: {e}")
            return False
    
    def _load_variants(self) -> Optional[pd.DataFrame]:
        """Load variants from ClinVar data files."""
        clinvar_file = self.input_dir / "clinvar" / "abca4_clinvar_vus.csv"
        
        if clinvar_file.exists():
            try:
                return pd.read_csv(clinvar_file)
            except Exception as e:
                logger.warning(f"Could not load {clinvar_file}: {e}")
        
        # Try parquet format
        parquet_file = self.input_dir / "clinvar" / "abca4_clinvar_vus.parquet"
        if parquet_file.exists():
            try:
                return pd.read_parquet(parquet_file)
            except Exception as e:
                logger.warning(f"Could not load {parquet_file}: {e}")
        
        logger.warning(f"No ClinVar files found in {self.input_dir / 'clinvar'}")
        return None
    
    def _filter_by_gene(self, df: pd.DataFrame) -> pd.DataFrame:
        """Filter variants to the gene of interest."""
        if "gene_symbol" in df.columns:
            df_filtered = df[df["gene_symbol"] == self.gene_name].copy()
        else:
            # If no gene_symbol column, assume all variants are for the gene
            df_filtered = df.copy()
        
        logger.info(f"Filtered to {len(df_filtered)} variants for {self.gene_name}")
        return df_filtered
    
    def _save_variants(self, df: pd.DataFrame) -> None:
        """Save filtered variants to output directory."""
        output_file = self.output_dir / f"{self.gene_name.lower()}_clinvar_vus.parquet"
        df.to_parquet(output_file, index=False)
        logger.info(f"Saved {len(df)} variants to {output_file}")


class ABCA4VariantFilter(ClinVarVariantFilter):
    """ABCA4-specific variant filtering."""
    
    def __init__(self, input_dir: Path, output_dir: Path):
        super().__init__(input_dir, output_dir, gene_name="ABCA4")
    
    def _load_variants(self) -> Optional[pd.DataFrame]:
        """Load ABCA4-specific ClinVar variants."""
        # Try to load from data_processed first (preprocessed)
        processed_file = Path("data_processed/variants/abca4_clinvar_vus.parquet")
        if processed_file.exists():
            try:
                df = pd.read_parquet(processed_file)
                logger.info(f"Loaded {len(df)} variants from {processed_file}")
                return df
            except Exception as e:
                logger.warning(f"Could not load {processed_file}: {e}")
        
        # Fall back to raw data
        return super()._load_variants()


def main():
    """CLI entry point for filtering ClinVar variants."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Filter ClinVar variants by gene")
    parser.add_argument("--gene", default="ABCA4", help="Gene symbol to filter to")
    parser.add_argument("--input-dir", type=Path, default="data_raw",
                       help="Input directory with ClinVar data")
    parser.add_argument("--output-dir", type=Path, default="data_processed/variants",
                       help="Output directory for filtered variants")
    
    args = parser.parse_args()
    
    if args.gene == "ABCA4":
        filter_engine = ABCA4VariantFilter(
            input_dir=args.input_dir,
            output_dir=args.output_dir
        )
    else:
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

