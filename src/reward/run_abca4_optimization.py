#!/usr/bin/env python3
"""Campaign-specific ranking pipeline that mimics a Strand optimization run."""

from __future__ import annotations

import logging
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd
import mlflow

CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT = CAMPAIGN_ROOT.parents[1]
sys.path.insert(0, str(REPO_ROOT))

from strand.core.sequence import Sequence
from strand.rewards.base import BaseRewardBlock, BlockType, RewardBlockMetadata, RewardContext

try:
    from strand.engine.engine import EngineConfig as SDKEngineConfig
    from strand.engine.engine import EngineResults as SDKEngineResults
except Exception:  # pragma: no cover - torch might be missing locally
    SDKEngineConfig = None
    SDKEngineResults = None


@dataclass(slots=True)
class LocalEngineConfig:
    iterations: int = 1
    population_size: int = 0
    seed: int = 0
    timeout_s: float = 0.0
    early_stop_patience: int | None = None
    max_evals: int | None = None
    method: str = "feature-ranking"
    extra: Mapping[str, object] = field(default_factory=dict)
    batching: Mapping[str, object] | None = None
    device: Mapping[str, object] | None = None


@dataclass(slots=True)
class LocalEngineResults:
    best: tuple[Sequence, float] | None = None
    history: list[object] = field(default_factory=list)
    summary: Mapping[str, object] = field(default_factory=dict)


EngineConfig = SDKEngineConfig or LocalEngineConfig
EngineResults = SDKEngineResults or LocalEngineResults


from src.config import logger


class FeatureRewardBlock(BaseRewardBlock):
    """Simple reward block that reads from the sequence metadata."""

    def __init__(self, name: str, feature: str, weight: float = 1.0) -> None:
        metadata = RewardBlockMetadata(
            block_type=BlockType.HEURISTIC,
            description=f"Feature-based reward for {feature}",
        )
        super().__init__(name=name, weight=weight, metadata=metadata)
        self.feature = feature

    def _score(self, sequence: Sequence, context: RewardContext) -> float:
        value = sequence.metadata.get(self.feature, 0.0)
        try:
            return float(value)
        except Exception:
            return 0.0


class OptimizationRunner:
    def __init__(self, feature_path: Path | None = None) -> None:
        features_dir = CAMPAIGN_ROOT / "data_processed" / "features"
        self.feature_matrix_path = feature_path or (features_dir / "abca4_feature_matrix.parquet")
        self.output_dir = features_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def load_feature_matrix(self) -> pd.DataFrame:
        if not self.feature_matrix_path.exists():
            raise FileNotFoundError(
                f"Feature matrix not found: {self.feature_matrix_path}. Run 'invoke features.compute' first."
            )
        df = pd.read_parquet(self.feature_matrix_path)
        logger.info("Loaded %s feature rows", len(df))
        return df

    def build_sequences(self, df: pd.DataFrame) -> list[Sequence]:
        sequences = []
        for _, row in df.iterrows():
            metadata = row.to_dict()
            seq = Sequence(id=row['variant_id'], tokens=row['variant_id'], metadata=metadata)
            sequences.append(seq)
        return sequences

    def compute_scores(self, sequences: list[Sequence]) -> pd.DataFrame:
        blocks = [
            FeatureRewardBlock("regulatory", "regulatory_score", weight=0.35),
            FeatureRewardBlock("splice", "spliceai_max_score", weight=0.25),
            FeatureRewardBlock("missense", "missense_combined_score", weight=0.25),
            FeatureRewardBlock("conservation", "conservation_score", weight=0.15),
        ]

        scores = []
        for seq in sequences:
            total = 0.0
            block_values = {}
            for block in blocks:
                value = block.score(seq, context=RewardContext())
                block_values[f"score_{block.name}"] = value
                total += value
            scores.append({
                'variant_id': seq.id,
                'reward': total,
                **block_values,
                **seq.metadata,
            })

        ranked = pd.DataFrame(scores).sort_values('reward', ascending=False).reset_index(drop=True)
        return ranked

    def log_mlflow(self, ranked: pd.DataFrame, top_k: int) -> None:
        mlflow.set_experiment("abca4-optimization")
        with mlflow.start_run(run_name=f"abca4_rank_{datetime.utcnow().isoformat(timespec='seconds')}"):
            config = EngineConfig(
                iterations=1,
                population_size=len(ranked),
                method="feature-ranking",
                extra={"top_k": top_k},
            )
            mlflow.log_params({
                "iterations": config.iterations,
                "population_size": config.population_size,
                "method": config.method,
                "top_k": top_k,
            })

            summary = {
                "top_reward": float(ranked.loc[0, 'reward']),
                "mean_top_reward": float(ranked.head(top_k)['reward'].mean()),
            }
            mlflow.log_metrics(summary)

            artifact_path = self.output_dir / "abca4_top_variants.json"
            ranked.head(top_k).to_json(artifact_path, orient='records', indent=2)
            mlflow.log_artifact(str(artifact_path))

    def save_outputs(self, ranked: pd.DataFrame, top_k: int) -> Path:
        output_path = self.output_dir / "abca4_ranked_variants.parquet"
        ranked.to_parquet(output_path, index=False)
        ranked.head(top_k).to_csv(self.output_dir / "abca4_top_variants.csv", index=False)
        logger.info("Saved ranked variants to %s", output_path)
        return output_path

    def run(self, top_k: int = 100) -> Path:
        df = self.load_feature_matrix()
        sequences = self.build_sequences(df)
        ranked = self.compute_scores(sequences)
        self.log_mlflow(ranked, top_k)
        return self.save_outputs(ranked, top_k)


def main() -> None:
    runner = OptimizationRunner()
    runner.run(top_k=100)


if __name__ == "__main__":
    main()
