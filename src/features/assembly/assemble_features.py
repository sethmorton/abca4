#!/usr/bin/env python3
"""Assemble feature matrix from individual feature tables."""

import logging
import sys
import os
from pathlib import Path
from typing import Dict, Optional

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[1]
CAMPAIGN_ROOT = PROJECT_ROOT
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.config import load_gene_config


class FeatureAssembler:
    """Combine annotation + feature tables into a single matrix."""

    def __init__(self, gene_name: str,
                 annotations_dir: Optional[Path] = None,
                 features_dir: Optional[Path] = None,
                 config: Optional[Dict] = None):
        if not gene_name or not isinstance(gene_name, str) or gene_name.strip() == "":
            raise ValueError("gene_name is required and must be a non-empty string")
        
        self.gene_name = gene_name
        self.config = config or load_gene_config(gene_name)
        self.file_prefix = self.config.get("file_prefix", self.gene_name.lower())
        
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
            self.annotations_dir / f"{self.file_prefix}_vus_annotated.parquet",
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

        output_path = self.features_dir / f"{self.file_prefix}_feature_matrix.parquet"
        csv_path = self.features_dir / f"{self.file_prefix}_feature_matrix.csv"
        try:
            merged.to_parquet(output_path, index=False)
            merged.to_csv(csv_path, index=False)
            logger.info("Saved unified feature matrix with %s rows", len(merged))
            return True
        except Exception as exc:
            logger.error("Unable to save feature matrix: %s", exc)
            return False


def main() -> None:
    import argparse
    
    parser = argparse.ArgumentParser(description="Assemble feature matrix from feature tables")
    parser.add_argument("--gene", type=str, default=os.getenv("GENE_NAME", None),
                       help="Gene symbol (required: pass via --gene or GENE_NAME env var)")
    args = parser.parse_args()
    
    if not args.gene:
        parser.error(
            "❌ ERROR: Gene symbol required but not provided.\n"
            "Please specify one of:\n"
            "  1. Command line: python script.py --gene GENE_NAME\n"
            "  2. Environment: export GENE_NAME=GENE_NAME"
        )
    
    try:
        assembler = FeatureAssembler(args.gene)
        success = assembler.run()
        raise SystemExit(0 if success else 1)
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
