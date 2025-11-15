#!/usr/bin/env python3
"""Download SpliceAI scores for the ABCA4 campaign.

Streams precomputed SpliceAI tracks from the public SpliceAI Lookup bucket,
filters to the ABCA4 window, and saves a tidy feature table for downstream
feature builders.
"""

import logging
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional
from urllib.parse import quote

import pandas as pd

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]


class SpliceAIDownloader:
    """Download and normalize SpliceAI features for the ABCA4 window."""

    BASE_GCS_URL = "https://storage.googleapis.com/download/storage/v1/b/tgg-viewer/o"
    DATASET_PREFIX = (
        "ref/GRCh38/spliceai/spliceai_scores.raw.snps_and_indels.hg38.filtered.sorted"
    )

    CHROM = "chr1"
    START = 93500000
    END = 95000000

    def __init__(self, output_dir: Optional[Path] = None,
                 score_threshold: float = 0.2,
                 effects: Optional[List[str]] = None):
        default_dir = CAMPAIGN_ROOT / "data_raw" / "spliceai"
        self.output_dir = output_dir or default_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.score_threshold = score_threshold
        self.score_label = self._format_threshold(score_threshold)
        self.effects = effects or ["gain", "loss"]
        self.region_no_prefix = f"{self.CHROM.lstrip('chr')}:{self.START}-{self.END}"

    @staticmethod
    def _format_threshold(value: float) -> str:
        text = f"{value:.2f}".rstrip('0').rstrip('.')
        return text or "0"

    def _build_track_url(self, effect: str) -> str:
        object_path = (
            f"{self.DATASET_PREFIX}.score_{self.score_label}.splice_{effect}.bed.gz"
        )
        encoded = quote(object_path, safe="")
        return f"{self.BASE_GCS_URL}/{encoded}?alt=media"

    def _stream_effect(self, effect: str) -> List[Dict[str, float]]:
        """Stream a precomputed SpliceAI track for the ABCA4 region."""
        url = self._build_track_url(effect)
        cmd = ["tabix", url, self.region_no_prefix]
        logger.info("Streaming %s track from %s", effect, url)

        try:
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
        except FileNotFoundError:
            logger.error("tabix is required but not installed")
            return []

        records: List[Dict[str, float]] = []
        assert process.stdout is not None
        for line in process.stdout:
            parsed = self._parse_track_line(line.strip(), effect)
            if parsed:
                records.append(parsed)

        stderr = process.stderr.read() if process.stderr else ""
        retcode = process.wait()

        if retcode != 0:
            logger.error("SpliceAI track streaming failed: %s", stderr.strip())
            return []

        logger.info("Fetched %s entries for %s", len(records), effect)
        return records

    @staticmethod
    def _parse_track_line(line: str, fallback_effect: str) -> Optional[Dict[str, float]]:
        if not line:
            return None

        parts = line.split('\t')
        if len(parts) < 4:
            return None

        info_blob = parts[3]
        info: Dict[str, str] = {}
        for chunk in info_blob.split(';'):
            if '=' in chunk:
                key, value = chunk.split('=', 1)
                info[key.strip()] = value.strip()

        allele = info.get('allele') or info.get('allele_with_max_score')
        if not allele:
            return None

        try:
            chrom, pos_str, ref, alt = allele.split('-', 3)
            pos = int(pos_str)
        except ValueError:
            return None

        score_str = info.get('score') or info.get('max_score')
        try:
            score_val = float(score_str) if score_str else 0.0
        except ValueError:
            score_val = 0.0

        effect_label = info.get('effect', f"{fallback_effect}")

        return {
            'chrom': chrom,
            'pos': pos,
            'ref': ref,
            'alt': alt,
            'effect': effect_label,
            'score': score_val,
        }

    def _aggregate_records(self, records: List[Dict[str, float]]) -> pd.DataFrame:
        if not records:
            return pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt'])

        df = pd.DataFrame(records)
        df = df.groupby(['chrom', 'pos', 'ref', 'alt', 'effect'], as_index=False)['score'].max()
        pivot = df.pivot_table(
            index=['chrom', 'pos', 'ref', 'alt'],
            columns='effect',
            values='score',
            aggfunc='max',
        ).reset_index().fillna(0.0)

        rename_map = {
            'acceptor_gain': 'spliceai_acceptor_gain',
            'acceptor_loss': 'spliceai_acceptor_loss',
            'donor_gain': 'spliceai_donor_gain',
            'donor_loss': 'spliceai_donor_loss',
        }

        columns_present = {col: rename_map.get(col, col) for col in pivot.columns}
        pivot = pivot.rename(columns=columns_present)

        for required in rename_map.values():
            if required not in pivot.columns:
                pivot[required] = 0.0

        feature_cols = [
            'spliceai_acceptor_gain',
            'spliceai_acceptor_loss',
            'spliceai_donor_gain',
            'spliceai_donor_loss',
        ]
        pivot['spliceai_max_score'] = pivot[feature_cols].max(axis=1)
        pivot['variant_id'] = (
            pivot['chrom'].astype(str) + '_' +
            pivot['pos'].astype(int).astype(str) + '_' +
            pivot['ref'].astype(str) + '_' +
            pivot['alt'].astype(str)
        )
        pivot['score_threshold'] = self.score_threshold
        return pivot

    def save_features(self, df: pd.DataFrame) -> bool:
        if df.empty:
            logger.warning("No SpliceAI entries were found for the requested region")

        output_path = self.output_dir / "spliceai_abca4_scores.parquet"
        csv_path = self.output_dir / "spliceai_abca4_scores.csv"

        try:
            df.to_parquet(output_path, index=False)
            df.to_csv(csv_path, index=False)
            logger.info("Saved normalized SpliceAI scores to %s", output_path)
            return True
        except Exception as exc:
            logger.error("Failed to save SpliceAI table: %s", exc)
            return False

    def run(self) -> bool:
        logger.info("Starting SpliceAI download and normalization ...")

        all_records: List[Dict[str, float]] = []
        for effect in self.effects:
            all_records.extend(self._stream_effect(effect))

        features_df = self._aggregate_records(all_records)
        return self.save_features(features_df)

def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(description="Download SpliceAI scores for ABCA4")
    parser.add_argument("--output-dir", type=Path,
                        help="Output directory (default: campaigns/abca4/data_raw/spliceai)")
    parser.add_argument("--score-threshold", type=float, default=0.2,
                        help="SpliceAI score threshold for the precomputed tracks (default: 0.2)")

    args = parser.parse_args()

    downloader = SpliceAIDownloader(
        output_dir=args.output_dir,
        score_threshold=args.score_threshold
    )
    success = downloader.run()
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
