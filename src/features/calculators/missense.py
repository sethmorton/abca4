#!/usr/bin/env python3
"""
Compute missense variant features from AlphaMissense and ESM.

Joins AlphaMissense pathogenicity scores with ESM embeddings for ABCA4 variants.
"""

import pandas as pd
from pathlib import Path
import sys
from typing import Optional, List
from src.config import logger

CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]

class MissenseFeatureComputer:
    """Compute missense-specific features for variants."""

    def __init__(self, input_dir: Optional[Path] = None, output_dir: Optional[Path] = None):
        processed_root = CAMPAIGN_ROOT / "data_processed"
        self.input_dir = input_dir or (processed_root / "annotations")
        self.raw_dir = CAMPAIGN_ROOT / "data_raw" / "alphamissense"
        self.output_dir = output_dir or (processed_root / "features")
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def load_annotated_variants(self) -> Optional[pd.DataFrame]:
        """Load annotated ABCA4 variants."""
        variants_path = self.input_dir / "abca4_vus_annotated.parquet"

        if not variants_path.exists():
            logger.error(f"Annotated variants not found: {variants_path}")
            return None

        logger.info(f"Loading annotated variants from {variants_path}")
        try:
            df = pd.read_parquet(variants_path)
            logger.info(f"Loaded {len(df)} annotated variants")
            return df
        except Exception as e:
            logger.error(f"Failed to load annotated variants: {e}")
            return None

    def load_alphamissense_data(self) -> Optional[pd.DataFrame]:
        """Load AlphaMissense data for ABCA4."""
        alphamissense_path = self.raw_dir / "alphamissense_abca4.tsv"

        if not alphamissense_path.exists():
            logger.error(f"AlphaMissense data not found: {alphamissense_path}")
            logger.info("Run 'invoke data.download' first to download AlphaMissense data")
            return None

        logger.info(f"Loading AlphaMissense data from {alphamissense_path}")
        try:
            df = pd.read_csv(alphamissense_path, sep='\t')
            logger.info(f"Loaded {len(df)} AlphaMissense entries")
            return df
        except Exception as e:
            logger.error(f"Failed to load AlphaMissense data: {e}")
            return None

    def _add_variant_key(self, df: pd.DataFrame) -> pd.Series:
        return (
            df['chrom'].astype(str).str.removeprefix('chr') + '_' +
            df['pos'].astype(int).astype(str) + '_' +
            df['ref'].astype(str) + '_' +
            df['alt'].astype(str)
        )

    def _build_alphamissense_key(self, alphamissense_df: pd.DataFrame) -> Optional[pd.Series]:
        chrom_col = next((c for c in ('chrom', 'Chromosome', '#CHROM') if c in alphamissense_df.columns), None)
        pos_col = next((c for c in ('pos', 'Position', 'Start') if c in alphamissense_df.columns), None)
        ref_col = next((c for c in ('ref', 'Ref', 'ReferenceAllele') if c in alphamissense_df.columns), None)
        alt_col = next((c for c in ('alt', 'Alt', 'AlternateAllele') if c in alphamissense_df.columns), None)

        if all([chrom_col, pos_col, ref_col, alt_col]):
            return (
                alphamissense_df[chrom_col].astype(str).str.removeprefix('chr') + '_' +
                alphamissense_df[pos_col].astype(int).astype(str) + '_' +
                alphamissense_df[ref_col].astype(str) + '_' +
                alphamissense_df[alt_col].astype(str)
            )
        return None

    def join_alphamissense_scores(self, variants_df: pd.DataFrame,
                                 alphamissense_df: pd.DataFrame) -> pd.DataFrame:
        """Join AlphaMissense scores with variants."""
        logger.info("Joining AlphaMissense scores with variants...")

        variants_with_scores = variants_df.copy()
        variants_with_scores['join_key'] = self._add_variant_key(variants_with_scores)

        alpha_df = alphamissense_df.copy()
        alpha_df['join_key'] = self._build_alphamissense_key(alpha_df)

        if alpha_df['join_key'].isna().all():
            protein_col = next((c for c in ('protein_variant', 'ProteinVariant', 'hgvsp', 'aa_change')
                                if c in alpha_df.columns), None)
            if not protein_col or 'protein_change' not in variants_with_scores.columns:
                logger.warning("Could not align AlphaMissense entries with variants; skipping join")
                return variants_with_scores

            alpha_df['protein_norm'] = alpha_df[protein_col].astype(str).str.replace('p.', '', regex=False)
            variants_with_scores['protein_norm'] = variants_with_scores['protein_change'].astype(str).str.replace('p.', '', regex=False)
            merge_cols = ['protein_norm']
        else:
            merge_cols = ['join_key']

        merge_cols = [col for col in merge_cols if col in alpha_df.columns and col in variants_with_scores.columns]

        merged = variants_with_scores.merge(
            alpha_df,
            on=merge_cols,
            how='left',
            suffixes=('', '_alpha')
        )

        score_col = next((c for c in ('am_pathogenicity', 'alphamissense_score') if c in merged.columns), None)
        class_col = next((c for c in ('am_class', 'alphamissense_class') if c in merged.columns), None)
        esm_col = next((c for c in ('esm1b_score', 'esm1b_log_likelihood_ratio', 'esm_score') if c in merged.columns), None)

        if score_col:
            merged['alphamissense_score'] = merged[score_col]
        if class_col:
            merged['alphamissense_class'] = merged[class_col]
        if esm_col:
            merged['esm_logit_diff'] = merged[esm_col]

        drop_cols = [col for col in ('join_key', 'protein_norm') if col in merged.columns]
        merged = merged.drop(columns=drop_cols)

        return merged

    def add_esm_features(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Ensure ESM-derived features exist (fallback to AlphaMissense metadata)."""
        variants_with_esm = variants_df.copy()
        if 'esm_logit_diff' not in variants_with_esm.columns:
            logger.info("No precomputed ESM scores found; deriving proxy from AlphaMissense")
            if 'alphamissense_score' in variants_with_esm.columns:
                variants_with_esm['esm_logit_diff'] = variants_with_esm['alphamissense_score'] - 0.5
            else:
                variants_with_esm['esm_logit_diff'] = 0.0

        return variants_with_esm

    def compute_missense_features(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Compute derived missense-specific features."""
        logger.info("Computing derived missense features...")

        features_df = variants_df.copy()

        score = features_df['alphamissense_score'].fillna(0.5)
        esm = features_df['esm_logit_diff'].fillna(0.0)

        features_df['alphamissense_pathogenic'] = (score > 0.8).astype(int)
        features_df['alphamissense_benign'] = (score < 0.2).astype(int)
        features_df['missense_combined_score'] = 0.7 * score + 0.3 * esm

        logger.info("Computed derived missense features")
        return features_df

    def save_missense_features(self, features_df: pd.DataFrame) -> bool:
        """Save missense features."""
        output_path = self.output_dir / "missense_features.parquet"

        try:
            # Select missense-specific columns
            missense_cols = [
                'variant_id', 'alphamissense_score', 'alphamissense_class',
                'alphamissense_pathogenic', 'alphamissense_benign',
                'esm_logit_diff', 'esm_embedding_dim', 'missense_combined_score'
            ]

            available_cols = [col for col in missense_cols if col in features_df.columns]
            missense_df = features_df[available_cols]

            missense_df.to_parquet(output_path, index=False)
            logger.info(f"Saved missense features for {len(missense_df)} variants to {output_path}")

            return True
        except Exception as e:
            logger.error(f"Failed to save missense features: {e}")
            return False

    def run(self) -> bool:
        """Run the complete missense feature computation process."""
        logger.info("Starting missense feature computation...")

        # Load annotated variants
        variants_df = self.load_annotated_variants()
        if variants_df is None:
            return False

        # Load AlphaMissense data
        alphamissense_df = self.load_alphamissense_data()
        if alphamissense_df is None:
            return False

        # Join AlphaMissense scores
        variants_df = self.join_alphamissense_scores(variants_df, alphamissense_df)

        # Add ESM features
        variants_df = self.add_esm_features(variants_df)

        # Compute derived features
        variants_df = self.compute_missense_features(variants_df)

        # Save results
        if not self.save_missense_features(variants_df):
            return False

        logger.info("Missense feature computation completed successfully!")
        return True

def main():
    """Main entry point."""
    computer = MissenseFeatureComputer()
    success = computer.run()
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
