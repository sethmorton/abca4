#!/usr/bin/env python3
"""Assemble the full ABCA4 feature matrix from individual feature tables."""

import logging
from pathlib import Path
from typing import Dict, Optional

import pandas as pd


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]


class FeatureAssembler:
    """Combine annotation + feature tables into a single matrix."""

    def __init__(self,
                 annotations_dir: Optional[Path] = None,
                 features_dir: Optional[Path] = None):
        processed_root = CAMPAIGN_ROOT / "data_processed"
        self.annotations_dir = annotations_dir or (processed_root / "annotations")
        self.features_dir = features_dir or (processed_root / "features")
        self.features_dir.mkdir(parents=True, exist_ok=True)

    def load_table(self, path: Path, description: str) -> Optional[pd.DataFrame]:
        if not path.exists():
            logger.warning("Skipping %s (missing %s)", description, path)
            return None
        try:
            df = pd.read_parquet(path)
            logger.info("Loaded %s rows from %s", len(df), path.name)
            return df
        except Exception as exc:
            logger.error("Unable to read %s: %s", path, exc)
            return None

    def run(self) -> bool:
        base = self.load_table(
            self.annotations_dir / "abca4_vus_annotated.parquet",
            "annotated variants",
        )
        if base is None:
            return False

        feature_paths: Dict[str, Path] = {
            'missense': self.features_dir / "missense_features.parquet",
            'splice': self.features_dir / "splice_features.parquet",
            'regulatory': self.features_dir / "regulatory_features.parquet",
            'conservation': self.features_dir / "conservation_features.parquet",
        }

        merged = base.copy()
        merged = merged.set_index('variant_id')

        for name, path in feature_paths.items():
            table = self.load_table(path, f"{name} features")
            if table is None:
                continue
            if 'variant_id' not in table.columns:
                logger.warning("%s table missing variant_id, skipping", name)
                continue
            table = table.drop_duplicates(subset='variant_id')
            merged = merged.join(table.set_index('variant_id'), how='left', rsuffix=f"_{name}")

        merged = merged.reset_index()
        merged = merged.fillna({
            'gnomad_genome_af': 0.0,
            'gnomad_exome_af': 0.0,
            'phyloP100way': 0.0,
            'phastCons100way': 0.0,
            'spliceai_max_score': 0.0,
        })

        output_path = self.features_dir / "abca4_feature_matrix.parquet"
        csv_path = self.features_dir / "abca4_feature_matrix.csv"
        try:
            merged.to_parquet(output_path, index=False)
            merged.to_csv(csv_path, index=False)
            logger.info("Saved unified feature matrix with %s rows", len(merged))
            return True
        except Exception as exc:
            logger.error("Unable to save feature matrix: %s", exc)
            return False


def main() -> None:
    assembler = FeatureAssembler()
    success = assembler.run()
    raise SystemExit(0 if success else 1)


if __name__ == "__main__":
    main()
