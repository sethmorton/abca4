#!/usr/bin/env python3
"""
Extract gene-specific variants from ClinVar VCF.

This script reads the ClinVar VCF file and extracts variants that are
annotated for a specific gene in the INFO field.

Usage:
    python extract_clinvar_variants.py --gene ABCA4
    python extract_clinvar_variants.py --gene TP53 --limit 500
    python extract_clinvar_variants.py --gene BRCA1 --vcf data_raw/clinvar/clinvar_20251116.vcf.gz
"""

import gzip
import pandas as pd
from pathlib import Path
from typing import List, Dict, Optional
import argparse
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def extract_gene_variants_from_vcf(
    vcf_path: Path,
    output_path: Path,
    gene_symbol: str,
    variant_limit: Optional[int] = None,
    filter_vus_only: bool = True
) -> int:
    """
    Extract gene-specific variants from ClinVar VCF file.
    
    Args:
        vcf_path: Path to clinvar_*.vcf.gz file
        output_path: Where to save the extracted variants as parquet
        gene_symbol: Gene symbol to extract (e.g., 'ABCA4', 'TP53', 'BRCA1')
        variant_limit: Maximum number of variants to extract (None = no limit)
        filter_vus_only: If True, only extract VUS (variants of uncertain significance)
        
    Returns:
        Number of variants extracted
    """
    
    variants = []
    count = 0
    skipped = 0
    
    logger.info(f"Extracting variants for gene: {gene_symbol}")
    logger.info(f"Reading from {vcf_path}...")
    
    if not vcf_path.exists():
        logger.error(f"VCF file not found: {vcf_path}")
        return 0
    
    with gzip.open(vcf_path, 'rt') as f:
        for line_num, line in enumerate(f):
            # Skip header lines
            if line.startswith('##'):
                continue
            
            if line.startswith('#CHROM'):
                continue
            
            # Parse VCF line
            parts = line.strip().split('\t')
            if len(parts) < 8:
                skipped += 1
                continue
            
            chrom, pos, var_id, ref, alt, qual, filt, info = parts[:8]
            
            # Parse INFO field to check for GENEINFO containing gene_symbol
            info_fields = {}
            for item in info.split(';'):
                if '=' in item:
                    key, val = item.split('=', 1)
                    info_fields[key] = val
            
            # Check if this variant is associated with the target gene
            geneinfo = info_fields.get('GENEINFO', '')
            if gene_symbol not in geneinfo:
                skipped += 1
                continue
            
            # Extract relevant fields
            clnsig = info_fields.get('CLNSIG', '')
            clnvc = info_fields.get('CLNVC', '')
            
            # Filter for VUS if requested
            if filter_vus_only:
                is_vus = 'Uncertain_significance' in clnsig or clnvc.lower() == 'single_nucleotide_variant'
                if not is_vus:
                    skipped += 1
                    continue
            
            # Create variant record
            variant = {
                'chrom': chrom,
                'pos': int(pos),
                'ref': ref,
                'alt': alt,
                'variant_id': f"var_{count + 1}",
                'clinvar_id': var_id,
                'clinsig': clnsig,
                'clnvc': clnvc,
            }
            
            variants.append(variant)
            count += 1
            
            # Log progress
            if (line_num + 1) % 100000 == 0:
                logger.info(f"  Processed {line_num + 1:,} lines, found {count} {gene_symbol} variants so far...")
            
            # Stop if we've reached the limit
            if variant_limit and count >= variant_limit:
                logger.info(f"Reached limit of {variant_limit} variants, stopping...")
                break
    
    if not variants:
        logger.warning(f"No {gene_symbol} variants found in VCF (processed {line_num + 1:,} lines, skipped {skipped})")
        return 0
    
    # Create DataFrame and save
    df = pd.DataFrame(variants)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(output_path, index=False)
    
    logger.info(f"✅ Extracted {count} {gene_symbol} variants")
    logger.info(f"   Saved to {output_path}")
    logger.info(f"   Skipped: {skipped}")
    
    return count

def main():
    """Main entry point with CLI argument parsing."""
    parser = argparse.ArgumentParser(
        description='Extract gene-specific variants from ClinVar VCF',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python extract_clinvar_variants.py --gene ABCA4
  python extract_clinvar_variants.py --gene TP53 --limit 500
  python extract_clinvar_variants.py --gene BRCA1 --no-vus-filter
        """
    )
    
    parser.add_argument(
        '--gene',
        type=str,
        required=True,
        help='Gene symbol to extract (e.g., ABCA4, TP53, BRCA1)'
    )
    parser.add_argument(
        '--vcf',
        type=Path,
        default=Path('data_raw/clinvar/clinvar_20251116.vcf.gz'),
        help='Path to ClinVar VCF file (default: data_raw/clinvar/clinvar_20251116.vcf.gz)'
    )
    parser.add_argument(
        '--output',
        type=Path,
        default=None,
        help='Output path (default: data_processed/variants/{gene}_clinvar_vus.parquet)'
    )
    parser.add_argument(
        '--limit',
        type=int,
        default=None,
        help='Maximum number of variants to extract (default: no limit)'
    )
    parser.add_argument(
        '--no-vus-filter',
        action='store_true',
        help='Extract all variants, not just VUS (variants of uncertain significance)'
    )
    
    args = parser.parse_args()
    
    # Set output path if not specified
    if args.output is None:
        args.output = Path('data_processed/variants') / f"{args.gene.lower()}_clinvar_vus.parquet"
    
    # Run extraction
    count = extract_gene_variants_from_vcf(
        vcf_path=args.vcf,
        output_path=args.output,
        gene_symbol=args.gene,
        variant_limit=args.limit,
        filter_vus_only=not args.no_vus_filter
    )
    
    if count > 0:
        # Show sample
        df = pd.read_parquet(args.output)
        logger.info(f"\nSample of extracted variants:")
        logger.info(f"{df.head()}")
        logger.info(f"\nDataFrame shape: {df.shape}")
    
    exit(0 if count > 0 else 1)

if __name__ == '__main__':
    main()

