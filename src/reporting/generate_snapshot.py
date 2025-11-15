#!/usr/bin/env python3
"""Build lightweight ABCA4 campaign snapshot for reporting."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(name)s | %(levelname)s | %(message)s",
)
LOGGER = logging.getLogger(__name__)

CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]
FEATURE_MATRIX = CAMPAIGN_ROOT / "data_processed" / "features" / "abca4_feature_matrix.parquet"
RANKED_VARIANTS = CAMPAIGN_ROOT / "data_processed" / "features" / "abca4_ranked_variants.parquet"
REPORT_DIR = CAMPAIGN_ROOT / "data_processed" / "reports"
REPORT_DIR.mkdir(parents=True, exist_ok=True)


def _load_feature_matrix() -> pd.DataFrame:
    if not FEATURE_MATRIX.exists():
        raise FileNotFoundError(
            f"Missing feature matrix: {FEATURE_MATRIX}. Run `invoke features.compute` first."
        )
    df = pd.read_parquet(FEATURE_MATRIX)
    LOGGER.info("Loaded %s feature rows", len(df))
    return df


def _load_ranked_variants(base_df: pd.DataFrame) -> pd.DataFrame:
    if RANKED_VARIANTS.exists():
        LOGGER.info("Using ranked variants from %s", RANKED_VARIANTS)
        return pd.read_parquet(RANKED_VARIANTS)

    LOGGER.warning("Ranked variants not found. Falling back to feature matrix ordering.")
    return base_df.assign(composite_score=0.0)


def _build_summary(df: pd.DataFrame, ranked: pd.DataFrame) -> dict:
    summary = {
        "total_variants": int(len(df)),
        "clin_sig_breakdown": df["clinical_significance"].value_counts().to_dict(),
        "coding_breakdown": df["coding_impact"].value_counts().to_dict(),
        "regulatory_breakdown": df["regulatory_type"].value_counts().to_dict(),
        "mean_gnomad_af": float(df.get("gnomad_max_af", pd.Series(dtype=float)).fillna(0).mean()),
        "median_spliceai": float(df.get("spliceai_max_score", pd.Series(dtype=float)).fillna(0).median()),
        "median_missense": float(df.get("missense_combined_score", pd.Series(dtype=float)).fillna(0).median()),
    }

    top_variants = ranked.head(50)
    summary["top_k"] = [
        {
            "variant_id": row["variant_id"],
            "clinical_significance": row.get("clinical_significance"),
            "regulatory_region": row.get("regulatory_region"),
            "gnomad_max_af": row.get("gnomad_max_af"),
            "spliceai_max_score": row.get("spliceai_max_score"),
            "missense_combined_score": row.get("missense_combined_score"),
            "composite_score": row.get("composite_score"),
        }
        for _, row in top_variants.iterrows()
    ]
    return summary


def _write_markdown(summary: dict) -> Path:
    md_path = REPORT_DIR / "abca4_snapshot.md"
    with open(md_path, "w") as handle:
        handle.write("# ABCA4 Campaign Snapshot\n\n")
        handle.write(f"**Total variants:** {summary['total_variants']:,}\n\n")
        handle.write("## Clinical significance\n")
        for label, count in summary["clin_sig_breakdown"].items():
            handle.write(f"- {label}: {count}\n")
        handle.write("\n## Coding impact\n")
        for label, count in summary["coding_breakdown"].items():
            handle.write(f"- {label}: {count}\n")
        handle.write("\n## Top variants\n")
        for idx, variant in enumerate(summary["top_k"], start=1):
            handle.write(
                f"{idx}. {variant['variant_id']} | clinsig={variant['clinical_significance']} | "
                f"reg={variant['regulatory_region']} | splice={variant['spliceai_max_score']:.3f} | "
                f"missense={variant['missense_combined_score']:.3f} | AF={variant['gnomad_max_af']:.4g}\n"
            )
    LOGGER.info("Wrote markdown snapshot to %s", md_path)
    return md_path


def main() -> None:
    df = _load_feature_matrix()
    ranked = _load_ranked_variants(df)
    summary = _build_summary(df, ranked)

    json_path = REPORT_DIR / "abca4_snapshot.json"
    with open(json_path, "w") as handle:
        json.dump(summary, handle, indent=2)
    LOGGER.info("Wrote summary JSON to %s", json_path)

    _write_markdown(summary)


if __name__ == "__main__":
    main()
