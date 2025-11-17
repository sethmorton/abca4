#!/usr/bin/env python3
"""Compute conservation features (phyloP/phastCons) for ABCA4 variants."""

from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import requests
from src.config import logger

CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]


class ConservationFeatureComputer:
    """Fetch UCSC conservation scores for ABCA4 variants."""

    BASE_URL = "https://api.genome.ucsc.edu/getData/track"
    TRACKS = ["phyloP100way", "phastCons100way"]
    CHROM = "chr1"
    CHUNK_SIZE = 200

    def __init__(self,
                 annotations_dir: Optional[Path] = None,
                 output_dir: Optional[Path] = None):
        processed_root = CAMPAIGN_ROOT / "data_processed"
        self.annotations_dir = annotations_dir or (processed_root / "annotations")
        self.output_dir = output_dir or (processed_root / "features")
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def load_variants(self) -> Optional[pd.DataFrame]:
        path = self.annotations_dir / "abca4_vus_annotated.parquet"
        if not path.exists():
            logger.error("Annotated variants not found at %s", path)
            return None
        try:
            df = pd.read_parquet(path)
            logger.info("Loaded %s annotated variants", len(df))
            return df
        except Exception as exc:
            logger.error("Unable to read annotated variants: %s", exc)
            return None

    def _request_track(self, track: str, start: int, end: int) -> Dict[int, float]:
        params = (
            f"genome=hg38;track={track};chrom={self.CHROM};start={start};end={end}"
        )
        url = f"{self.BASE_URL}?{params}"
        try:
            resp = requests.get(url, timeout=30)
            if resp.status_code != 200:
                logger.warning("UCSC request failed (%s): %s", resp.status_code, url)
                return {}
            payload = resp.json()
            entries = payload.get(self.CHROM, [])
        except Exception as exc:
            logger.error("Failed to fetch %s: %s", track, exc)
            return {}

        values: Dict[int, float] = {}
        for entry in entries:
            entry_start = entry.get('start')
            entry_end = entry.get('end', entry_start + 1)
            value = entry.get('value', 0.0)
            if entry_start is None:
                continue
            for pos in range(entry_start, entry_end):
                values[pos] = float(value)
        return values

    def fetch_track(self, track: str, positions: List[int]) -> Dict[int, float]:
        positions = sorted(set(positions))
        scores: Dict[int, float] = {}
        for i in range(0, len(positions), self.CHUNK_SIZE):
            chunk = positions[i:i + self.CHUNK_SIZE]
            if not chunk:
                continue
            window_start = max(chunk[0] - 5, 0)
            window_end = chunk[-1] + 6
            logger.info("Fetching %s for %s-%s", track, window_start, window_end)
            data = self._request_track(track, window_start, window_end)
            for pos in chunk:
                if pos in data:
                    scores[pos] = data[pos]
        return scores

    def compute(self, variants: pd.DataFrame) -> pd.DataFrame:
        variants = variants.copy()
        positions = variants['pos'].astype(int).tolist()

        phyloP = self.fetch_track('phyloP100way', positions)
        phastCons = self.fetch_track('phastCons100way', positions)

        variants['phyloP100way'] = variants['pos'].map(phyloP).fillna(0.0)
        variants['phastCons100way'] = variants['pos'].map(phastCons).fillna(0.0)

        variants['phyloP_high'] = (variants['phyloP100way'] > 2.0).astype(int)
        variants['phastCons_high'] = (variants['phastCons100way'] > 0.8).astype(int)
        variants['conservation_rank'] = variants['phyloP100way'].rank(pct=True)

        for column in ['phyloP100way', 'phastCons100way']:
            mean = variants[column].mean()
            std = variants[column].std(ddof=0) or 1.0
            variants[f"{column}_z"] = (variants[column] - mean) / std

        variants['conservation_score'] = (
            variants['phyloP100way_z'] * 0.6 + variants['phastCons100way_z'] * 0.4
        )

        return variants

    def save(self, df: pd.DataFrame) -> bool:
        output_path = self.output_dir / "conservation_features.parquet"
        try:
            df.to_parquet(output_path, index=False)
            logger.info("Saved %s conservation rows", len(df))
            return True
        except Exception as exc:
            logger.error("Unable to save conservation features: %s", exc)
            return False

    def run(self) -> bool:
        variants = self.load_variants()
        if variants is None:
            return False

        features = self.compute(variants)
        columns = [
            'variant_id', 'phyloP100way', 'phastCons100way',
            'phyloP100way_z', 'phastCons100way_z',
            'phyloP_high', 'phastCons_high', 'conservation_rank',
            'conservation_score'
        ]
        available = [col for col in columns if col in features.columns]
        return self.save(features[available])


def main() -> None:
    computer = ConservationFeatureComputer()
    success = computer.run()
    raise SystemExit(0 if success else 1)


if __name__ == "__main__":
    main()
