#!/usr/bin/env python3
"""
Add proper Step 5 clustering to variants_scored.parquet.

Ensures every variant has:
- cluster_id: Domain or consequence-based cluster assignment
- cluster_target: τⱼ coverage threshold per cluster
- coverage_by_cluster: cov_j(S) = max impact_score in cluster
"""

import sys
from pathlib import Path
from typing import Optional

import pandas as pd
from src.config import logger

CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]


class ClusteringProcessor:
    """Add Step 5 clustering to scored variants."""

    def __init__(self, features_dir: Optional[Path] = None):
        processed_root = CAMPAIGN_ROOT / "data_processed"
        self.features_dir = features_dir or (processed_root / "features")
        self.features_dir.mkdir(parents=True, exist_ok=True)

    def load_scored_variants(self) -> Optional[pd.DataFrame]:
        """Load variants_scored.parquet."""
        scored_path = self.features_dir / "variants_scored.parquet"
        if not scored_path.exists():
            logger.error(f"Missing variants_scored.parquet at {scored_path}")
            return None
        
        try:
            df = pd.read_parquet(scored_path)
            logger.info(f"Loaded {len(df)} scored variants")
            return df
        except Exception as e:
            logger.error(f"Failed to load scored variants: {e}")
            return None

    def assign_clusters_consequence(self, df: pd.DataFrame) -> pd.DataFrame:
        """Assign clusters based on variant consequence (primary mechanism)."""
        df_clustered = df.copy()
        
        def consequence_to_cluster(consequence_str):
            """Map VEP consequence to cluster category."""
            if not isinstance(consequence_str, str) or pd.isna(consequence_str):
                return "other"
            
            consequence_lower = consequence_str.lower()
            
            # Check for LoF (high priority)
            if any(x in consequence_lower for x in ["frameshift", "stop_gained", "stop_lost"]):
                if "frameshift" in consequence_lower:
                    return "LoF_frameshift"
                return "LoF_stop"
            
            # Check for canonical splice variants
            if any(x in consequence_lower for x in ["splice_acceptor", "splice_donor"]):
                if "acceptor" in consequence_lower:
                    return "LoF_splice_acceptor"
                return "LoF_splice_donor"
            
            # Check for missense
            if "missense" in consequence_lower:
                return "missense"
            
            # Check for splice region/proxy (moderate impact)
            if any(x in consequence_lower for x in ["splice_region", "splice_polypyrimidine"]):
                return "splice_region"
            
            # Check for indels
            if any(x in consequence_lower for x in ["inframe_deletion", "inframe_insertion", "disruptive"]):
                return "inframe_indel"
            
            # Check for synonymous
            if "synonymous" in consequence_lower:
                return "synonymous"
            
            # Check for intron
            if "intron" in consequence_lower:
                return "intron"
            
            # Check for UTR
            if "utr" in consequence_lower:
                return "utr"
            
            return "other"
        
        if "vep_consequence" in df_clustered.columns:
            df_clustered["cluster_id"] = df_clustered["vep_consequence"].apply(consequence_to_cluster)
            logger.info("Using vep_consequence for mechanism-based clustering")
        else:
            df_clustered["cluster_id"] = "unknown"
            logger.warning("vep_consequence not found; using 'unknown' for all")
        
        return df_clustered

    def compute_cluster_targets(self, df: pd.DataFrame, threshold_factor: float = 0.8) -> pd.DataFrame:
        """Compute τⱼ = threshold_factor × max(model_score) per cluster."""
        df_with_targets = df.copy()
        
        if "model_score" not in df_with_targets.columns:
            logger.warning("model_score not found; using random fallback")
            df_with_targets["model_score"] = 0.5
        
        # Ensure model_score is numeric and bounded [0, 1]
        df_with_targets["model_score"] = pd.to_numeric(
            df_with_targets["model_score"], errors="coerce"
        ).fillna(0.5).clip(0, 1)
        
        cluster_targets = {}
        for cluster_id, group in df_with_targets.groupby("cluster_id"):
            max_score = group["model_score"].max()
            tau_j = max_score * threshold_factor
            cluster_targets[cluster_id] = tau_j
            logger.info(f"  {cluster_id}: max_score={max_score:.4f}, τⱼ={tau_j:.4f}")
        
        df_with_targets["cluster_target"] = df_with_targets["cluster_id"].map(cluster_targets)
        
        logger.info(f"Computed cluster targets for {len(cluster_targets)} clusters (threshold={threshold_factor})")
        return df_with_targets

    def compute_cluster_coverage(self, df: pd.DataFrame) -> pd.DataFrame:
        """Compute cov_j(S) = max(model_score) in each cluster."""
        df_with_coverage = df.copy()
        
        cluster_coverage = {}
        for cluster_id, group in df_with_coverage.groupby("cluster_id"):
            max_score = group["model_score"].max()
            cluster_coverage[cluster_id] = max_score
            logger.info(f"  {cluster_id}: cov_j(S)={max_score:.4f}")
        
        df_with_coverage["coverage_by_cluster"] = df_with_coverage["cluster_id"].map(cluster_coverage)
        
        logger.info(f"Computed coverage metrics for {len(cluster_coverage)} clusters")
        return df_with_coverage

    def save_clustered_variants(self, df: pd.DataFrame) -> bool:
        """Save variants with clustering info to variants_scored.parquet."""
        output_path = self.features_dir / "variants_scored.parquet"
        
        try:
            df.to_parquet(output_path, index=False)
            logger.info(f"Saved {len(df)} clustered variants to {output_path}")
            
            # Log summary
            if "cluster_id" in df.columns:
                cluster_counts = df["cluster_id"].value_counts()
                logger.info(f"Cluster distribution:\n{cluster_counts}")
            
            return True
        except Exception as e:
            logger.error(f"Failed to save clustered variants: {e}")
            return False

    def run(self, threshold_factor: float = 0.8) -> bool:
        """Execute full clustering pipeline."""
        logger.info("Starting clustering processor...")
        
        # Load
        df = self.load_scored_variants()
        if df is None:
            return False
        
        # Assign clusters based on consequence/mechanism
        logger.info("Assigning clusters...")
        df = self.assign_clusters_consequence(df)
        
        # Compute targets
        logger.info("Computing cluster targets...")
        df = self.compute_cluster_targets(df, threshold_factor)
        
        # Compute coverage
        logger.info("Computing cluster coverage...")
        df = self.compute_cluster_coverage(df)
        
        # Save
        logger.info("Saving results...")
        return self.save_clustered_variants(df)


def main():
    """Main entry point."""
    processor = ClusteringProcessor()
    success = processor.run(threshold_factor=0.8)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
