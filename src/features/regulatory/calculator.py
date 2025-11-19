#!/usr/bin/env python3
"""
Compute regulatory heuristics for ABCA4 variants.

Derives scores based on domain importance, TSS proximity, and population frequency.
"""

import json
import pandas as pd
import numpy as np
import pysam
from pathlib import Path
import sys
from typing import Optional, Dict, Any
from src.config import logger, ENSEMBL_RELEASE, CANONICAL_TRANSCRIPT
from .constants import (
    PRIORITY_SCORE_WEIGHT, TSS_DISTANCE_WEIGHT, RARITY_WEIGHT,
    REGULATORY_CLIP_THRESHOLD, TSS_WINDOW_SIZE, MAX_TSS_SCORE
)

CAMPAIGN_ROOT = Path(__file__).resolve().parents[3]

class RegulatoryFeatureComputer:
    """Compute regulatory heuristics for variants."""

    def __init__(self, input_dir: Optional[Path] = None, output_dir: Optional[Path] = None):
        processed_root = CAMPAIGN_ROOT / "data_processed"
        self.input_dir = input_dir or (processed_root / "annotations")
        self.raw_dir = CAMPAIGN_ROOT / "data_raw"
        self.output_dir = output_dir or (processed_root / "features")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.domains_config = self._load_domains_config()

    def _load_domains_config(self) -> Dict[str, Any]:
        """Load domain configuration."""
        config_path = self.raw_dir / "abca4_domains.json"
        if not config_path.exists():
            logger.warning(f"Domain config not found at {config_path}, using empty config")
            return {}
        
        try:
            with open(config_path, 'r') as f:
                return json.load(f)
        except Exception as e:
            logger.error(f"Failed to load domain config: {e}")
            return {}

    def load_annotated_variants(self) -> Optional[pd.DataFrame]:
        """Load annotated ABCA4 variants."""
        variants_path = self.input_dir / "abca4_vus_annotated.parquet"

        if not variants_path.exists():
            logger.error(f"Annotated variants not found: {variants_path}")
            return None

        logger.info(f"Loading annotated variants from {variants_path}")
        try:
            df = pd.read_parquet(variants_path)
            return df
        except Exception as e:
            logger.error(f"Failed to load annotated variants: {e}")
            return None

    def _calculate_domain_priority(self, row: pd.Series) -> float:
        """Calculate priority score based on domain location."""
        # This is a heuristic placeholder. In a real scenario, this would check
        # if the variant falls within critical domains defined in self.domains_config.
        # For now, we assign higher priority to exonic variants in key domains.
        
        consequence = str(row.get('vep_consequence', '')).lower()
        
        if 'missense' in consequence:
            return 0.8
        elif 'stop' in consequence or 'frameshift' in consequence:
            return 1.0
        elif 'splice' in consequence:
            return 0.9
        elif 'synonymous' in consequence:
            return 0.1
        elif 'intron' in consequence:
            return 0.2
        
        return 0.5

    def _calculate_tss_score(self, pos: int) -> float:
        """Calculate score based on distance to Transcription Start Site."""
        # ABCA4 TSS (approximate for hg38)
        tss_pos = 94458394 
        distance = abs(pos - tss_pos)
        
        if distance > TSS_WINDOW_SIZE:
            return 0.0
        
        # Linear decay
        return MAX_TSS_SCORE * (1 - (distance / TSS_WINDOW_SIZE))

    def _calculate_rarity_score(self, af: float) -> float:
        """Calculate rarity score (inverse of allele frequency)."""
        if pd.isna(af):
            return 1.0
        
        # Logarithmic scaling for rarity
        if af <= 0:
            return 1.0
            
        return min(1.0, -0.1 * np.log10(af))

    def compute_regulatory_features(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Compute derived regulatory features."""
        logger.info("Computing derived regulatory features...")

        features_df = variants_df.copy()
        
        # 1. Domain Priority
        features_df['domain_priority_score'] = features_df.apply(self._calculate_domain_priority, axis=1)
        
        # 2. TSS Proximity
        features_df['tss_proximity_score'] = features_df['pos'].apply(self._calculate_tss_score)
        
        # 3. Rarity Score
        # Use gnomad_max_af if available, else 0
        af_col = 'gnomad_max_af' if 'gnomad_max_af' in features_df.columns else 'gnomAD_AF'
        if af_col in features_df.columns:
            features_df['rarity_score'] = features_df[af_col].apply(self._calculate_rarity_score)
        else:
            features_df['rarity_score'] = 1.0 # Assume rare if unknown

        # Combined Regulatory Score
        features_df['regulatory_score'] = (
            PRIORITY_SCORE_WEIGHT * features_df['domain_priority_score'] +
            TSS_DISTANCE_WEIGHT * features_df['tss_proximity_score'] +
            RARITY_WEIGHT * features_df['rarity_score']
        )
        
        # Clip low scores
        features_df.loc[features_df['regulatory_score'] < REGULATORY_CLIP_THRESHOLD, 'regulatory_score'] = 0.0

        logger.info("Computed derived regulatory features")
        return features_df

    def save_regulatory_features(self, features_df: pd.DataFrame) -> bool:
        """Save regulatory features."""
        output_path = self.output_dir / "regulatory_features.parquet"

        try:
            # Select regulatory-specific columns
            reg_cols = [
                'variant_id', 'domain_priority_score', 'tss_proximity_score',
                'rarity_score', 'regulatory_score'
            ]

            available_cols = [col for col in reg_cols if col in features_df.columns]
            reg_df = features_df[available_cols]

            reg_df.to_parquet(output_path, index=False)
            logger.info(f"Saved regulatory features for {len(reg_df)} variants to {output_path}")

            return True
        except Exception as e:
            logger.error(f"Failed to save regulatory features: {e}")
            return False

    def run(self) -> bool:
        """Run the complete regulatory feature computation process."""
        logger.info("Starting regulatory feature computation...")

        # Load annotated variants
        variants_df = self.load_annotated_variants()
        if variants_df is None:
            return False

        # Compute derived features
        variants_df = self.compute_regulatory_features(variants_df)

        # Save results
        if not self.save_regulatory_features(variants_df):
            return False

        logger.info("Regulatory feature computation completed successfully!")
        return True

def main():
    """Main entry point."""
    computer = RegulatoryFeatureComputer()
    success = computer.run()
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
