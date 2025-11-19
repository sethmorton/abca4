#!/usr/bin/env python3
"""
Compute protein domain assignments for variants.

Uses UniProt domain coordinates (provided via config) to map coding variants to their domains.
For non-coding variants, uses consequence as fallback (e.g., "intronic", "utr").

Domain configuration is gene-specific and loaded from config at runtime.
"""

import re
import sys
from pathlib import Path
from typing import Optional, Tuple, Dict

import pandas as pd
from src.config import logger
from src.config import load_gene_config

CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]


class DomainMapper:
    """
    Map variants to protein domains based on gene configuration.
    
    This class is gene-agnostic. Domain boundaries and boost factors are loaded
    from gene-specific config at initialization. If domains are not defined in
    config, the mapper falls back to consequence-based domain assignment.
    """

    def __init__(self, gene_name: str, features_dir: Optional[Path] = None, config: Optional[Dict] = None):
        """
        Initialize domain mapper.
        
        Args:
            gene_name: Gene symbol (e.g., "ABCA4") - REQUIRED, no default
            features_dir: Directory to save outputs
            config: Optional pre-loaded gene config. If not provided, loads from config file.
        
        Raises:
            ValueError: If gene_name is empty or domains not defined in config
        """
        if not gene_name or not isinstance(gene_name, str) or gene_name.strip() == "":
            raise ValueError("gene_name is required and must be a non-empty string")
        
        self.gene_name = gene_name
        self.config = config or load_gene_config(gene_name)
        self.domains = self.config.get("domains", {})
        
        if not self.domains:
            raise ValueError(
                f"No domains defined in config for {gene_name}. "
                f"Domain assignment requires domain boundaries in config/abca4.yaml"
            )
        
        processed_root = CAMPAIGN_ROOT / "data_processed"
        self.features_dir = features_dir or (processed_root / "features")
        self.features_dir.mkdir(parents=True, exist_ok=True)

    def load_annotated_variants(self) -> Optional[pd.DataFrame]:
        """Load annotated variants using config path."""
        annotated_path = CAMPAIGN_ROOT / self.config.get("output_paths", {}).get("annotated_variants", f"data_processed/annotations/{self.gene_name.lower()}_vus_annotated.parquet")
        if not annotated_path.exists():
            logger.error(f"Missing annotated variants at {annotated_path}")
            return None
        
        try:
            df = pd.read_parquet(annotated_path)
            logger.info(f"Loaded {len(df)} annotated variants for {self.gene_name}")
            return df
        except Exception as e:
            logger.error(f"Failed to load annotated variants for {self.gene_name}: {e}")
            return None

    def extract_protein_position(self, protein_change_str: str) -> Optional[int]:
        """
        Extract amino acid position from HGVS protein notation.
        
        Examples:
          "ENSP00000359245.3:p.Ala1028Val" → 1028
          "p.M123V" → 123
          "p.Gly*" → extraction fails → None
        """
        if not isinstance(protein_change_str, str) or pd.isna(protein_change_str):
            return None
        
        # Try to extract position after "p." or just look for digits in context
        match = re.search(r'p\.[A-Z][a-z]*(\d+)', protein_change_str)
        if match:
            try:
                return int(match.group(1))
            except (ValueError, IndexError):
                return None
        
        return None

    def position_to_domain(self, protein_position: int) -> str:
        """Map protein position to domain based on config."""
        for domain_name, (start, end) in self.domains.items():
            if start <= protein_position <= end:
                return domain_name
        
        # Out of bounds (rare, shouldn't happen with valid variants)
        return "other"

    def consequence_to_fallback_domain(self, consequence_str: str) -> str:
        """Map consequence to fallback domain for non-coding variants."""
        if not isinstance(consequence_str, str) or pd.isna(consequence_str):
            return "unknown"
        
        consequence_lower = consequence_str.lower()
        
        # Non-coding variants
        if "intron" in consequence_lower:
            return "intronic"
        if "utr" in consequence_lower:
            return "utr"
        if "splice" in consequence_lower and "region" in consequence_lower:
            return "splice_region"
        if "downstream" in consequence_lower or "upstream" in consequence_lower:
            return "regulatory"
        if "synonymous" in consequence_lower:
            return "synonymous"
        
        return "other"

    def compute_domains(self, df: pd.DataFrame) -> pd.DataFrame:
        """Compute domain assignment for all variants."""
        df_with_domains = df.copy()
        
        domains = []
        for idx, row in df_with_domains.iterrows():
            # Try to extract from protein change first (most reliable)
            protein_change = row.get('protein_change')
            protein_pos = self.extract_protein_position(protein_change)
            
            if protein_pos is not None:
                # Coding variant with extractable position
                domain = self.position_to_domain(protein_pos)
                domains.append(domain)
            else:
                # Non-coding or position not extractable - use consequence
                consequence = row.get('vep_consequence')
                domain = self.consequence_to_fallback_domain(consequence)
                domains.append(domain)
        
        df_with_domains['domain'] = domains
        
        # Log summary
        domain_dist = pd.Series(domains).value_counts()
        logger.info("Domain distribution:")
        for domain, count in domain_dist.items():
            pct = 100 * count / len(domains)
            logger.info(f"  {domain}: {count} ({pct:.1f}%)")
        
        return df_with_domains

    def save_with_domains(self, df: pd.DataFrame) -> bool:
        """Save variants with domain annotations."""
        file_prefix = self.config.get("file_prefix", self.gene_name.lower())
        output_path = self.features_dir / f"{file_prefix}_variants_with_domains.parquet"
        
        try:
            df.to_parquet(output_path, index=False)
            logger.info(f"Saved variants with domains to {output_path}")
            return True
        except Exception as e:
            logger.error(f"Failed to save: {e}")
            return False

    def get_hotspot_mapping(self) -> dict:
        """
        Get mapping of known hotspot positions to override domain assignments.
        Used by notebooks for mechanism-aware clustering.

        Returns dict of (start_pos, end_pos): override_domain
        """
        # Example Stargardt hotspots - to be filled with literature-based positions
        return {
            # (588, 588): 'NBD1_hotspot',  # Example: p.Gly196Arg hotspot
            # (768, 768): 'TMD1_hotspot',  # Example: p.Gly863Ala hotspot
            # Add more based on literature review
        }

    def apply_hotspot_overrides(self, df: pd.DataFrame, hotspot_mapping: dict = None) -> pd.DataFrame:
        """
        Apply hotspot overrides to domain assignments.
        Used by notebooks to override domains for known mechanism hotspots.
        """
        if hotspot_mapping is None:
            hotspot_mapping = self.get_hotspot_mapping()

        df_with_overrides = df.copy()

        if 'pos' in df_with_overrides.columns and hotspot_mapping:
            for (start_pos, end_pos), override_domain in hotspot_mapping.items():
                pos_mask = (df_with_overrides['pos'] >= start_pos) & (df_with_overrides['pos'] <= end_pos)
                df_with_overrides.loc[pos_mask, 'domain'] = override_domain
                logger.info(f"Applied hotspot override: positions {start_pos}-{end_pos} → {override_domain}")

        return df_with_overrides

    def run(self) -> bool:
        """Execute full domain computation pipeline."""
        logger.info("Starting domain mapper...")

        # Load
        df = self.load_annotated_variants()
        if df is None:
            return False

        # Compute domains
        logger.info("Computing domain assignments...")
        df = self.compute_domains(df)

        # Save
        logger.info("Saving results...")
        return self.save_with_domains(df)


def main():
    """Main entry point."""
    mapper = DomainMapper()
    success = mapper.run()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
