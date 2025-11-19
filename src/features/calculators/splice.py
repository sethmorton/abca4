#!/usr/bin/env python3
"""
Compute splice variant features from SpliceAI.

Extracts and processes SpliceAI scores for ABCA4 variants.
"""

import pandas as pd
from pathlib import Path
import sys
from typing import Optional
from src.config import logger

CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]

class SpliceFeatureComputer:
    """Compute splice-specific features for variants."""

    def __init__(self, input_dir: Optional[Path] = None, output_dir: Optional[Path] = None):
        processed_root = CAMPAIGN_ROOT / "data_processed"
        self.input_dir = input_dir or (processed_root / "annotations")
        self.raw_dir = CAMPAIGN_ROOT / "data_raw" / "spliceai"
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

    def load_spliceai_data(self) -> Optional[pd.DataFrame]:
        """Load SpliceAI data for ABCA4 region."""
        tidy_path = self.raw_dir / "spliceai_abca4_scores.parquet"
        legacy_vcf = self.raw_dir / "spliceai_abca4_scores.vcf"

        if tidy_path.exists():
            try:
                df = pd.read_parquet(tidy_path)
                logger.info(f"Loaded {len(df)} SpliceAI entries from {tidy_path}")
                return df
            except Exception as exc:
                logger.error(f"Failed to read {tidy_path}: {exc}")

        if not legacy_vcf.exists():
            logger.error("SpliceAI data not found. Run 'invoke data.download' first.")
            return None

        logger.info(f"Loading SpliceAI data from {legacy_vcf}")
        try:
            spliceai_records = []

            with open(legacy_vcf, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue

                    fields = line.strip().split('\t')
                    if len(fields) < 8:
                        continue

                    chrom, pos, _id, ref, alt, _qual, _filter, info = fields[:8]

                    info_dict = {}
                    for item in info.split(';'):
                        if '=' in item:
                            key, value = item.split('=', 1)
                            info_dict[key] = value

                    spliceai_blob = info_dict.get('SpliceAI')
                    if not spliceai_blob:
                        continue

                    scores = self._parse_spliceai_blob(spliceai_blob, alt)
                    if not scores:
                        continue

                    spliceai_records.append({
                        'chrom': chrom,
                        'pos': int(pos),
                        'ref': ref,
                        'alt': alt,
                        **scores
                    })

            df = pd.DataFrame(spliceai_records)
            logger.info(f"Loaded {len(df)} SpliceAI predictions from legacy VCF")
            return df

        except Exception as e:
            logger.error(f"Failed to load SpliceAI data: {e}")
            return None

    @staticmethod
    def _parse_spliceai_blob(blob: str, alt: str) -> Optional[dict]:
        entries = blob.split(',')
        selected_entry = None
        for entry in entries:
            parts = entry.split('|')
            if len(parts) < 7:
                continue
            allele = parts[0]
            if allele == alt:
                selected_entry = parts
                break
            if selected_entry is None:
                selected_entry = parts  # fallback to first entry

        if not selected_entry or len(selected_entry) < 7:
            return None

        def _to_float(value: str) -> float:
            try:
                return float(value)
            except (TypeError, ValueError):
                return 0.0

        acceptor_gain = _to_float(selected_entry[3])
        acceptor_loss = _to_float(selected_entry[4])
        donor_gain = _to_float(selected_entry[5])
        donor_loss = _to_float(selected_entry[6])

        return {
            'spliceai_acceptor_gain': acceptor_gain,
            'spliceai_acceptor_loss': acceptor_loss,
            'spliceai_donor_gain': donor_gain,
            'spliceai_donor_loss': donor_loss,
            'spliceai_max_score': max(acceptor_gain, acceptor_loss, donor_gain, donor_loss)
        }

    def join_spliceai_scores(self, variants_df: pd.DataFrame,
                           spliceai_df: pd.DataFrame) -> pd.DataFrame:
        """Join SpliceAI scores with variants."""
        logger.info("Joining SpliceAI scores with variants...")

        # Create join key from variant coordinates
        variants_df = variants_df.copy()
        variants_df['join_key'] = (
            variants_df['chrom'].astype(str) + '_' +
            variants_df['pos'].astype(str) + '_' +
            variants_df['ref'] + '_' +
            variants_df['alt']
        )

        spliceai_df = spliceai_df.copy()
        spliceai_df['join_key'] = (
            spliceai_df['chrom'].astype(str) + '_' +
            spliceai_df['pos'].astype(str) + '_' +
            spliceai_df['ref'] + '_' +
            spliceai_df['alt']
        )

        # Left join to keep all variants
        merged_df = variants_df.merge(
            spliceai_df,
            on='join_key',
            how='left',
            suffixes=('', '_spliceai')
        )

        # Fill missing scores with 0
        score_cols = [
            'spliceai_acceptor_gain', 'spliceai_acceptor_loss',
            'spliceai_donor_gain', 'spliceai_donor_loss', 'spliceai_max_score'
        ]

        for col in score_cols:
            if col in merged_df.columns:
                merged_df[col] = merged_df[col].fillna(0)

        # Remove temporary join key
        merged_df = merged_df.drop('join_key', axis=1)

        matched_variants = merged_df['spliceai_max_score'].notna().sum()
        logger.info(f"Matched {matched_variants} variants with SpliceAI scores")

        return merged_df

    def compute_splice_features(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Compute derived splice-specific features."""
        logger.info("Computing derived splice features...")

        features_df = variants_df.copy()

        # SpliceAI pathogenicity thresholds
        features_df['spliceai_pathogenic'] = (
            features_df['spliceai_max_score'] > 0.5
        ).astype(int)

        features_df['spliceai_high_impact'] = (
            features_df['spliceai_max_score'] > 0.8
        ).astype(int)

        # Splice impact categories
        conditions = [
            (features_df['spliceai_max_score'] < 0.1),
            (features_df['spliceai_max_score'] < 0.5),
            (features_df['spliceai_max_score'] < 0.8),
            (features_df['spliceai_max_score'] >= 0.8)
        ]
        choices = ['low', 'moderate', 'high', 'very_high']
        features_df['spliceai_impact_category'] = pd.cut(
            features_df['spliceai_max_score'],
            bins=[0, 0.1, 0.5, 0.8, 1.0],
            labels=['low', 'moderate', 'high', 'very_high'],
            include_lowest=True
        )

        # Delta scores (gain vs loss)
        features_df['spliceai_net_gain'] = (
            features_df['spliceai_acceptor_gain'] + features_df['spliceai_donor_gain'] -
            features_df['spliceai_acceptor_loss'] - features_df['spliceai_donor_loss']
        )

        logger.info("Computed derived splice features")
        return features_df

    def save_splice_features(self, features_df: pd.DataFrame) -> bool:
        """Save splice features."""
        output_path = self.output_dir / "splice_features.parquet"

        try:
            # Select splice-specific columns
            splice_cols = [
                'variant_id', 'spliceai_acceptor_gain', 'spliceai_acceptor_loss',
                'spliceai_donor_gain', 'spliceai_donor_loss', 'spliceai_max_score',
                'spliceai_pathogenic', 'spliceai_high_impact',
                'spliceai_impact_category', 'spliceai_net_gain'
            ]

            available_cols = [col for col in splice_cols if col in features_df.columns]
            splice_df = features_df[available_cols]

            splice_df.to_parquet(output_path, index=False)
            logger.info(f"Saved splice features for {len(splice_df)} variants to {output_path}")

            return True
        except Exception as e:
            logger.error(f"Failed to save splice features: {e}")
            return False

    def run(self) -> bool:
        """Run the complete splice feature computation process."""
        logger.info("Starting splice feature computation...")

        # Load annotated variants
        variants_df = self.load_annotated_variants()
        if variants_df is None:
            return False

        # Load SpliceAI data
        spliceai_df = self.load_spliceai_data()
        if spliceai_df is None:
            return False

        # Join SpliceAI scores
        variants_df = self.join_spliceai_scores(variants_df, spliceai_df)

        # Compute derived features
        variants_df = self.compute_splice_features(variants_df)

        # Save results
        if not self.save_splice_features(variants_df):
            return False

        logger.info("Splice feature computation completed successfully!")
        return True

def main():
    """Main entry point."""
    computer = SpliceFeatureComputer()
    success = computer.run()
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
